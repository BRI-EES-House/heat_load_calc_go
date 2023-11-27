package heat_load_calc

// 応答係数の初項、指数項別応答係数、公比の計算

import (
	"math"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

const nRoot = 12

type ResponseFactor struct {
	rft0 float64   // 貫流応答係数の初項
	rfa0 float64   // 吸熱応答係数の初項
	rft1 []float64 // 貫流応答係数
	rfa1 []float64 // 吸熱応答係数
	row  []float64 // 公比
}

func NewResponseFactor(
	rft0 float64,
	rfa0 float64,
	rft1 []float64,
	rfa1 []float64,
	row []float64,
) *ResponseFactor {
	return &ResponseFactor{
		rft0: rft0,
		rfa0: rfa0,
		rft1: rft1,
		rfa1: rfa1,
		row:  row,
	}
}

/*
	Args:
		u_w: 熱貫流率, W/m2K
		r_i: 室内側（総合）熱伝達抵抗, m2K/W

	Returns:
		応答係数
*/
func create_for_steady(u_w float64, r_i float64) *ResponseFactor {
	// 開口部の室内表面から屋外までの熱コンダクタンス, W/m2K
	u_so := 1.0 / (1.0/u_w - r_i)

	return NewResponseFactor(
		1.0,
		1.0/u_so,
		make([]float64, nRoot),
		make([]float64, nRoot),
		make([]float64, nRoot),
	)
}

/*
	応答係数を作成する（地盤以外に用いる）

	裏面に、熱抵抗をもち、熱容量は 0.0 の層を追加する。
	Args:
		cs: 単位面積あたりの熱容量, kJ/m2K, [layer数]
		rs: 熱抵抗, m2K/W, [layer数]
		r_o: 室外側熱伝達抵抗, m2K/W

	Returns:
		応答係数
*/
func create_for_unsteady_not_ground(cs []float64, rs []float64, r_o float64) *ResponseFactor {
	// 裏面に熱容量 0.0 、熱抵抗 r_o の層を加える。
	cs = append(cs[:], 0.0)
	rs = append(rs[:], r_o)

	// 単位変換 kJ/m2K -> J/m2K
	floats.Scale(1000.0, cs)

	// 応答係数
	rft0, rfa0, rft1, rfa1, row := calc_response_factor_non_residential(cs, rs)

	return NewResponseFactor(rft0, rfa0, rft1, rfa1, row)
}

/*
	応答係数を作成する（地盤用）

	地盤層を追加する。
	Args:
		cs: 単位面積あたりの熱容量, kJ/m2K, [layer数]
		rs: 熱抵抗, m2K/W, [layer数]

	Returns:
		応答係数
*/
func create_for_unsteady_ground(cs []float64, rs []float64) *ResponseFactor {
	// 裏面に地盤の層を加える。
	cs = append(cs[:], 3300.0*3.0)
	rs = append(rs[:], 3.0/1.0)

	// 単位変換 kJ/m2K -> J/m2K
	floats.Scale(1000.0, cs)

	// 応答係数
	rft0, rfa0, rft1, rfa1, row := calc_response_factor(true, cs, rs)

	// 貫流応答係数の上書
	// 土壌の計算は吸熱応答のみで計算するため、畳み込み積分に必要な指数項別応答係数はすべて０にする
	// 貫流応答の初項は年平均気温掛かる係数であることから１とし、計算された貫流応答係数をすべて上書きする
	rft0 = 1.0
	rft1 = make([]float64, nRoot)

	return NewResponseFactor(rft0, rfa0, rft1, rfa1, row)
}

// ラプラス変数の設定
// alp: 固定根
func get_laps(alp []float64) []float64 {
	// 与えるラプラス変数の個数
	n := len(alp) * 2

	// ラプラス変数の配列
	laps := make([]float64, n)

	for i := 1; i <= n; i++ {
		if i%2 == 0 {
			// 偶数番目はαをそのまま入力
			laps[i-1] = alp[int((i-1)/2)]
		} else if i == 1 {
			// 最初はα1/√(α2/α1）とする
			laps[i-1] = alp[0] / math.Sqrt(alp[1]/alp[0])
		} else {
			// それ以外は等比数列で補間
			lngL := int(math.Ceil(float64(i-1) / 2))
			laps[i-1] = alp[lngL] / math.Sqrt(alp[lngL]/alp[lngL-1])
		}
	}

	return laps
}

/*
   固定根を取得する。

   Args:
       is_ground: 熱容量を持つ外皮が地盤かどうか。（地盤でない場合は、一般部位または間仕切り）

   Returns:
       固定根
*/
func get_alpha_m(is_ground bool) []float64 {
	if is_ground {
		// 地盤の場合
		return []float64{
			1.05699306612549e-08,
			3.27447457677204e-08,
			1.01440436059147e-07,
			3.14253839100315e-07,
			9.73531652917036e-07,
			3.01591822058485e-06,
			9.34305801562961e-06,
			2.89439987091205e-05,
			8.96660450863221e-05,
			2.77777777777778e-04,
		}
	} else {
		// 地盤以外の場合
		return []float64{
			2.0000e-06,
			7.9000e-06,
			3.1205e-05,
			1.2325e-04,
			4.8687e-04,
			1.9231e-03,
			7.5964e-03,
			3.0006e-02,
		}
	}
}

func _sum(a []float64) float64 {
	var sum float64
	for _, v := range a {
		sum += v
	}
	return sum
}

// 壁体の単位応答の計算
// Args:
// 		C_i_k_p: 壁体の熱容量, kJ/m2K, [layer数]
// 		R_i_k_p: 壁体の熱抵抗, m2K/W, [layer数]
// 		laps : ラプラス変数
// 		alp : 固定根
// 		M : 積分の上限
func get_step_reps_of_wall(C_i_k_p, R_i_k_p, laps, alp []float64, M int) (float64, float64, []float64, []float64) {
	// 吸熱、貫流の各伝達関数ベクトルの初期化
	nlaps := len(laps)
	matGA := mat.NewVecDense(nlaps, nil)
	matGT := mat.NewVecDense(nlaps, nil)

	// 吸熱、貫流の各伝達関数ベクトルの初期化
	dblAT0 := 1.0
	dblAA0 := _sum(R_i_k_p)

	// GA(0), GT(0)
	dblGA0 := dblAA0
	dblGT0 := dblAT0

	// 壁体の熱容量が0（定常）の場合
	// 定常部位であっても、そのまま処理を継続する（計算上は問題ないため）
	// if abs(dblCtotal) < 0.001:
	//    pass //　暫定処理（VBAではここで処理を抜ける）

	// 吸熱、貫流の各伝達関数ベクトルの作成
	for lngI := 0; lngI < len(laps); lngI++ {
		// 伝達関数の計算
		GA, GT := calc_transfer_function(C_i_k_p, R_i_k_p, laps[lngI])

		// 吸熱、貫流の各伝達関数ベクトルの作成
		matGA.SetVec(lngI, GA-dblGA0)
		matGT.SetVec(lngI, GT-dblGT0)
	}

	// 伝達関数の係数を求めるための左辺行列を作成
	nroot := len(alp)
	matF := mat.NewDense(nlaps, nroot, nil)
	for lngI, lap := range laps {
		for lngJ, root := range alp {
			matF.Set(lngI, lngJ, lap/(lap+root))
		}
	}

	// 最小二乗法のための係数行列を作成
	var matU mat.Dense
	matU.Mul(matF.T(), matF)

	// 最小二乗法のための定数項行列を作成
	var matFGA, matFGT mat.Dense
	matFGA.Apply(func(i, j int, v float64) float64 { return v * matGA.AtVec(i) }, matF)
	matFGT.Apply(func(i, j int, v float64) float64 { return v * matGT.AtVec(i) }, matF)
	CA := make([]float64, nroot)
	CT := make([]float64, nroot)
	for lngI := 0; lngI < nroot; lngI++ {
		matFGA.RawRowView(lngI)
		for lngJ := 0; lngJ < nlaps; lngJ++ {
			CA[lngI] += matFGA.At(lngJ, lngI)
			CT[lngI] += matFGT.At(lngJ, lngI)
		}
	}
	matCA := mat.NewVecDense(nroot, CA)
	matCT := mat.NewVecDense(nroot, CT)

	// 伝達関数の係数を計算
	var matU_LU mat.LU
	var matAA, matAT mat.VecDense
	matU_LU.Factorize(&matU)
	err1 := matU_LU.SolveVecTo(&matAA, false, matCA)
	if err1 != nil {
		panic(err1)
	}
	err2 := matU_LU.SolveVecTo(&matAT, false, matCT)
	if err2 != nil {
		panic(err2)
	}

	// 伝達関数の係数を一次元配列に変換
	dblAT := matAT.RawVector()
	dblAA := matAA.RawVector()

	if dblAT.Inc != 1 {
		panic("matAT.RawVector().Inc != 1")
	}
	if dblAA.Inc != 1 {
		panic("matAA.RawVector().Inc != 1")
	}

	return dblAT0, dblAA0, dblAT.Data, dblAA.Data
}

/*
	伝達関数の計算
    Args:
        C_i_k_p: 層の熱容量[J/m2K]
        R_i_k_p: 層の熱抵抗[m2K/W]
        laps: ラプラス変数[1/s]

    Returns:
        貫流伝達関数[K]
*/
func calc_transfer_function(C_i_k_p []float64, R_i_k_p []float64, laps float64) (float64, float64) {
	// 四端子行列の初期化
	matFt := mat.NewDense(2, 2, []float64{
		1.0, 0.0,
		0.0, 1.0,
	})
	var matFi *mat.Dense
	for lngK, R_k := range R_i_k_p {
		C_k := C_i_k_p[lngK]

		// ---- 四端子基本行列 matFi ----
		if math.Abs(C_k) < 0.001 {
			// 定常部位（空気層等）の場合
			matFi = mat.NewDense(2, 2, []float64{
				1.0, R_k,
				0.0, 1.0,
			})
		} else {
			// 非定常部位の場合
			dblTemp := math.Sqrt(R_k * C_k * laps)
			dblCosh := math.Cosh(dblTemp)
			dblSinh := math.Sinh(dblTemp)

			matFi = mat.NewDense(2, 2, []float64{
				dblCosh, R_k / dblTemp * dblSinh,
				dblTemp / R_k * dblSinh, dblCosh,
			})
		}

		// ---- 四端子行列 matFt ----
		matFt.Mul(matFt, matFi)
	}

	// 吸熱、貫流の各伝達関数ベクトルの作成
	GA := matFt.At(0, 1) / matFt.At(1, 1)
	GT := 1.0 / matFt.At(1, 1)

	return GA, GT
}

// 壁体の単位応答の計算（非住宅向け重み付き最小二乗法適用）
// Args:
//	C_i_k_p: 熱容量の配列
//	R_i_k_p: 熱抵抗の配列
//	laps: 時刻の配列
//	alp: 温度時定数の配列
//	M: 応答係数で作成する項数
func get_step_reps_of_wall_weighted(C_i_k_p, R_i_k_p []float64, laps []float64, alp []float64, weight float64) (float64, float64, []float64, []float64) {
	// 四端子基本行列の初期化
	matFi := make([]*mat.Dense, len(C_i_k_p))

	// 吸熱、貫流の各伝達関数ベクトルの初期化
	nlaps := len(laps)
	matGA := make([]float64, nlaps)
	matGT := make([]float64, nlaps)

	// 単位貫流応答、単位吸熱応答の初期化
	dblAT0 := 1.0
	dblAA0 := _sum(R_i_k_p)

	if nlaps == 0.0 {
		return dblAT0, dblAA0, matGA, matGT
	}

	// GA(0), GT(0)
	dblGA0 := dblAA0
	dblGT0 := dblAT0

	// 吸熱、貫流の各伝達関数ベクトルの作成
	matFt := mat.NewDense(2, 2, nil)
	for lngI := 0; lngI < len(laps); lngI++ {
		for lngK, R_k := range R_i_k_p {
			C_k := C_i_k_p[lngK]

			// ---- 四端子基本行列 matFi ----
			var matFiElem *mat.Dense
			if math.Abs(C_k) < 0.001 {
				// 定常部位（空気層等）の場合
				matFiElem = mat.NewDense(2, 2, []float64{
					1.0, R_k,
					0.0, 1.0,
				})
			} else {
				// 非定常部位の場合
				dblTemp := math.Sqrt(R_k * C_k * laps[lngI])
				dblCosh := math.Cosh(dblTemp)
				dblSinh := math.Sinh(dblTemp)

				matFiElem = mat.NewDense(2, 2, []float64{
					dblCosh, R_k / dblTemp * dblSinh,
					dblTemp / R_k * dblSinh, dblCosh,
				})
			}

			matFi[lngK] = matFiElem

			// ---- 四端子行列 matFt ----
			if lngK == 0 {
				// 室内側1層目の場合は、四端子行列に四端子基本行列をコピーする
				matFt.CloneFrom(matFiElem)
			} else {
				// 室内側2層目以降は、四端子基本行列を乗算
				matFt.Mul(matFt, matFiElem)
			}
		}

		// 吸熱、貫流の各伝達関数ベクトルの作成
		matGA[lngI] = matFt.At(0, 1)/matFt.At(1, 1) - dblGA0
		matGT[lngI] = 1.0/matFt.At(1, 1) - dblGT0
	}

	// 伝達関数の係数を求めるための左辺行列を作成
	nroot := len(alp)
	matF := mat.NewDense(nlaps, nroot, nil)
	for lngI, lap := range laps {
		for lngJ, root := range alp {
			matF.Set(lngI, lngJ, lap/(lap+root))
		}
	}

	// 最小二乗法のための係数行列を作成
	matU := mat.NewDense(nroot, nroot, nil)
	for lngK := 0; lngK < nroot; lngK++ {
		for lngJ := 0; lngJ < nroot; lngJ++ {
			for lngI := 0; lngI < nlaps; lngI++ {
				matU.Set(lngK, lngJ, matU.At(lngK, lngJ)+math.Pow(laps[lngI], weight)*matF.At(lngI, lngK)*matF.At(lngI, lngJ))
			}
		}
	}

	// 最小二乗法のための定数項行列を作成
	matCA := mat.NewVecDense(nroot, nil)
	matCT := mat.NewVecDense(nroot, nil)
	for lngK := 0; lngK < nroot; lngK++ {
		for lngI := 0; lngI < nlaps; lngI++ {
			laps_w_F := math.Pow(laps[lngI], weight) * matF.At(lngI, lngK)
			matCA.SetVec(lngK, matCA.At(lngK, 0)+laps_w_F*matGA[lngI])
			matCT.SetVec(lngK, matCT.At(lngK, 0)+laps_w_F*matGT[lngI])
		}
	}

	// 伝達関数の係数を計算
	var matU_LU mat.LU
	var matAA, matAT mat.VecDense
	matU_LU.Factorize(matU)
	err1 := matU_LU.SolveVecTo(&matAA, false, matCA)
	if err1 != nil {
		panic(err1)
	}
	err2 := matU_LU.SolveVecTo(&matAT, false, matCT)
	if err2 != nil {
		panic(err2)
	}

	// 伝達関数の係数を一次元配列に変換
	dblAT := matAT.RawVector()
	dblAA := matAA.RawVector()

	if dblAT.Inc != 1 {
		panic("matAT.RawVector().Inc != 1")
	}
	if dblAA.Inc != 1 {
		panic("matAA.RawVector().Inc != 1")
	}

	return dblAT0, dblAA0, dblAT.Data, dblAA.Data
}

// 二等辺三角波励振の応答係数、指数項別応答係数、公比の計算
func get_RFTRI(alp []float64, AT0 float64, AA0 float64, AT []float64, AA []float64, M int) ([]float64, []float64, []float64, []float64, []float64) {
	// 二等辺三角波励振の応答係数の配列を初期化
	dblRFT := make([]float64, M)
	dblRFA := make([]float64, M)

	// 二等辺三角波励振の応答係数の初項を計算
	dblTemp := make([]float64, len(alp))
	for i, val := range alp {
		dblTemp[i] = val * 900
	}

	// 二等辺三角波励振の応答係数の二項目以降を計算
	dblE1 := make([]float64, len(alp))
	for i, val := range dblTemp {
		dblE1[i] = (1.0 - math.Exp(-val)) / val
	}
	dblRFT[0] = AT0 + floats.Dot(dblE1, AT)
	dblRFA[0] = AA0 + floats.Dot(dblE1, AA)

	// 二等辺三角波励振の応答係数の二項目以降を計算
	for lngJ := 1; lngJ < M; lngJ++ {
		for i, val := range dblTemp {
			temp := (1.0 - math.Exp(-val))
			dblE1[i] = temp * temp * math.Exp(-(float64(lngJ)-1.0)*val) / val
		}
		dblRFT[lngJ] = -floats.Dot(dblE1, AT)
		dblRFA[lngJ] = -floats.Dot(dblE1, AA)
	}

	// 指数項別応答係数、公比を計算
	for i, val := range dblTemp {
		temp := (1.0 - math.Exp(-val))
		dblE1[i] = temp * temp / val
	}
	dblRFT1 := make([]float64, len(AT))
	dblRFA1 := make([]float64, len(AA))
	for i := range AT {
		dblRFT1[i] = -AT[i] * dblE1[i]
		dblRFA1[i] = -AA[i] * dblE1[i]
	}
	dblRow := make([]float64, len(alp))
	for i, val := range dblTemp {
		dblRow[i] = math.Exp(-val)
	}

	return dblRFT, dblRFA, dblRFT1, dblRFA1, dblRow
}

// 応答係数
func calc_response_factor(is_ground bool, cs, rs []float64) (float64, float64, []float64, []float64, []float64) {
	NcalTime := 50                  // 応答係数を作成する時間数[h]
	M := int(NcalTime*3600/900) + 1 //応答係数で作成する項数

	//  固定根, 一般部位の場合[8], 地盤の場合[10]
	alpha_m := get_alpha_m(is_ground)

	// ラプラス変数の設定
	laps := get_laps(alpha_m)

	// 単位応答の計算
	AT0, AA0, AT, AA := get_step_reps_of_wall(cs, rs, laps, alpha_m, M)

	// 二等辺三角波励振の応答係数、指数項別応答係数、公比の計算
	RFT, RFA, RFT1, RFA1, Row := get_RFTRI(alpha_m, AT0, AA0, AT, AA, M)

	RFT0 := RFT[0] // 貫流応答係数の初項
	RFA0 := RFA[0] // 吸熱応答係数の初項
	//Nroot := len(alpha_m) // 根の数

	RFT1_12 := make([]float64, 12)
	RFA1_12 := make([]float64, 12)
	Row_12 := make([]float64, 12)
	copy(RFT1_12, RFT1)
	copy(RFA1_12, RFA1)
	copy(Row_12, Row)

	return RFT0, RFA0, RFT1_12, RFA1_12, Row_12
}

func logspace(start, stop float64, num int, base float64) []float64 {
	if num <= 0 {
		return []float64{}
	}

	step := (stop - start) / float64(num-1)
	logspaceSlice := make([]float64, num)

	for i := 0; i < num; i++ {
		logspaceSlice[i] = math.Pow(base, start+float64(i)*step)
	}

	return logspaceSlice
}

// 応答係数（非住宅用　住宅との相違は固定根と重み付き最小二乗法を使用する点）
func calc_response_factor_non_residential(C_i_k_p, R_i_k_p []float64) (float64, float64, []float64, []float64, []float64) {
	NcalTime := 50                           // 応答係数を作成する時間数[h]
	M := int(float64(NcalTime*3600)/900) + 1 // 応答係数で作成する項数

	// 固定根, 初項 1/(86400*365)、終項 1/600、項数 10
	first_term := math.Log10(1.0 / (86400.0 * 365.0))
	last_term := math.Log10(1.0 / 900.0)
	alpha_m := logspace(first_term, last_term, 10, 10)
	alpha_m_temp := make([]float64, len(alpha_m))
	copy(alpha_m_temp, alpha_m)

	nroot := len(alpha_m) // 根の数
	// 実際に応答係数計算に使用する固定根を選定する
	GA := make([]float64, nroot)
	GT := make([]float64, nroot)
	// 固定根をラプラスパラメータとして伝達関数を計算
	for i, lap := range alpha_m_temp {
		GA[i], GT[i] = calc_transfer_function(C_i_k_p, R_i_k_p, lap)
	}
	GT2 := make([]float64, nroot+2)
	// 配列0に定常の伝達関数を入力
	GT2[0] = 1.0
	// 配列の最後にs=∞の伝達関数を入力
	GT2[nroot+2-1] = 0.0
	// それ以外に計算した伝達関数を代入
	for i := 0; i < nroot; i++ {
		GT2[i+1] = GT[i]
	}

	// 採用する固定根の場合1
	is_adopts := make([]float64, nroot)
	for i := 0; i <= nroot+1; i++ {
		for j := i + 1; j < nroot+1; j++ {
			// 伝達関数が3%以上変化した根だけ採用する
			if math.Abs(GT2[j]-GT2[i]) > 0.03 {
				is_adopts[j-1] = 1.0
				i = j - 1
				break
			}
		}
	}

	//採用する根を数える
	adopt_nroot := floats.Sum(is_adopts)
	if adopt_nroot > 0.0 {
		// 不採用の固定根を削除
		floats.Reverse(is_adopts)
		for i, adopts := range is_adopts {
			if adopts == 0.0 {
				alpha_m_temp = append(alpha_m_temp[:nroot-i-1], alpha_m_temp[nroot-i:]...)
			}
		}
	} else {
		alpha_m_temp = []float64{}
	}

	// 採用する根を数える
	//adopt_nroot := floats.Sum(is_adopts)

	// ラプラス変数の設定
	laps := get_laps(alpha_m_temp)

	// 単位応答の計算
	AT0, AA0, AT, AA := get_step_reps_of_wall_weighted(C_i_k_p, R_i_k_p, laps, alpha_m_temp, 0.0)

	// 二等辺三角波励振の応答係数、指数項別応答係数、公比の計算
	RFT, RFA, RFT1, RFA1, Row := get_RFTRI(alpha_m_temp, AT0, AA0, AT, AA, M)

	RFT0 := RFT[0] // 貫流応答係数の初項
	RFA0 := RFA[0] // 吸熱応答係数の初項
	//Nroot := len(alpha_m) // 根の数

	RFT1_12 := make([]float64, 12)
	RFA1_12 := make([]float64, 12)
	Row_12 := make([]float64, 12)

	for i, alpha := range alpha_m {
		for j, alpha_temp := range alpha_m_temp {
			if alpha == alpha_temp {
				RFT1_12[i] = RFT1[j]
				RFA1_12[i] = RFA1[j]
				Row_12[i] = Row[j]
			}
		}
	}

	return RFT0, RFA0, RFT1_12, RFA1_12, Row_12
}
