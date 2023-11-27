package heat_load_calc

import (
	"gonum.org/v1/gonum/mat"
)

type Equipments struct {
	is_radiative_heating_is []bool
	is_radiative_cooling_is []bool
	q_rs_h_max_is           []float64
	q_rs_c_max_is           []float64
	beta_h_is               []float64
	beta_c_is               []float64
	has_flr                 bool
	f_flr_h_js_is           mat.Matrix // 放射成分(床暖房が無い場合はゼロ行列)
	f_flr_c_js_is           mat.Matrix // 放射成分(床暖房が無い場合はゼロ行列)→常のゼロでは？
	_hes                    []interface{}
	_ces                    []interface{}
	_n_rm                   int
	_n_b                    int
}

type EquipmentBase struct {
	id      int    // ID
	name    string //名前
	room_id int    //暖房・冷房する室のID
}

type EquipmentRAC struct {
	q_min float64 //最小暖房・冷房能力, W
	q_max float64 //最大暖房・冷房能力, W
	v_min float64 //最小風量, m3/min
	v_max float64 //最大風量, m3/min
	bf    float64 //バイパスファクター
}

type HeatingEquipmentRAC struct {
	EquipmentBase
	EquipmentRAC
}

type EquipmentFloor struct {
	boundary_id      int     // 放射暖房・冷房が設置されている境界の番号
	max_capacity     float64 // 面積当たりの放熱能力, W/m2
	area             float64 // 面積, m2
	convection_ratio float64 // 対流成分比率, -
}

type HeatingEquipmentFloorHeating struct {
	EquipmentBase
	EquipmentFloor
}

type CoolingEquipmentRAC struct {
	EquipmentBase
	EquipmentRAC
}

type CoolingEquipmentFloorCooling struct {
	EquipmentBase
	EquipmentFloor
}

type EquipmentInfo struct {
	HeatingEquipments []EquipmentData `json:"heating_equipments"`
	CoolingEquipments []EquipmentData `json:"cooling_equipments"`
}

// NOTE: 構造体のフィールドに対応するJSONキーが存在しない場合、そのフィールドにはゼロ値が設定されます。
type EquipmentData struct {
	EquipmentType string `json:"equipment_type"`
	Id            int    `json:"id"`
	Name          string `json:"name"`
	Property      struct {
		SpaceId         int     `json:"space_id"`
		QMin            float64 `json:"q_min"`
		QMax            float64 `json:"q_max"`
		VMin            float64 `json:"v_min"`
		VMax            float64 `json:"v_max"`
		Bf              float64 `json:"bf"`
		BoundaryId      int     `json:"boundary_id"`
		MaxCapacity     float64 `json:"max_capacity"`
		Area            float64 `json:"area"`
		ConvectionRatio float64 `json:"convection_ratio"`
	} `json:"property"`
}

/*
設備に関する情報を辞書形式で受け取り、データクラスに変換して保持する。
暖房・冷房それぞれにおいて、
辞書の中の "equipment_type" の種類に応じて対応するデータクラスを生成する。

Args:

	dict_equipments: 設備の情報が記された辞書
	n_rm: 部屋の数
	n_b: 境界の数
	bs: Boundariesクラス

Notes:

	ここで Boundaries クラスは、境界IDと室IDとの対応関係を見ることだけに使用される。
	放射暖冷房に関する設備情報には対応する境界IDしか記されていない。
	一方で、放射暖冷房においても、beta, f_flr の係数を計算する際には、
	その放射暖冷房がどの室に属しているのかの情報が必要になるため、
	Equipments を initialize する際に、あらかじめ放射暖冷房にも room_id を付与しておくこととする。
*/
func NewEquipments(eq *EquipmentsJson, n_rm int, n_b int, bs *Boundaries) *Equipments {
	heating_eq := eq.HeatingEquipments
	hes := make([]interface{}, len(heating_eq))
	for i, he := range heating_eq {
		hes[i] = _create_heating_equipment(&he, bs)
	}

	cooling_eq := eq.CoolingEquipments
	ces := make([]interface{}, len(cooling_eq))
	for i, ce := range cooling_eq {
		ces[i] = _create_cooling_equipment(&ce, bs)
	}

	f_flr_h_js_is := _get_f_flr_js_is(hes, n_rm, n_b)
	f_flr_c_js_is := _get_f_flr_js_is(ces, n_rm, n_b)
	has_flr := false
	r, c := f_flr_h_js_is.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			if f_flr_h_js_is.At(i, j) != 0.0 {
				has_flr = true
				break
			}
			if f_flr_c_js_is.At(i, j) != 0.0 {
				has_flr = true
				break
			}
		}
		if has_flr {
			break
		}
	}
	return &Equipments{
		_hes:                    hes,
		_ces:                    ces,
		is_radiative_heating_is: _get_is_radiative_is(hes, n_rm),
		is_radiative_cooling_is: _get_is_radiative_is(ces, n_rm),
		q_rs_h_max_is:           _get_q_rs_max_is(hes, n_rm),
		q_rs_c_max_is:           _get_q_rs_max_is(ces, n_rm),
		beta_h_is:               _get_beta_is(hes, n_rm),
		beta_c_is:               _get_beta_is(ces, n_rm),
		has_flr:                 has_flr,
		f_flr_h_js_is:           f_flr_h_js_is,
		f_flr_c_js_is:           f_flr_c_js_is,
		_n_rm:                   n_rm,
		_n_b:                    n_b,
	}
}

func _create_heating_equipment(eq *EquipmentJson, bs *Boundaries) interface{} {
	heType := eq.EquipmentType
	id := eq.Id
	name := eq.Name
	prop := eq.Property

	if heType == "rac" {
		return &HeatingEquipmentRAC{
			EquipmentBase: EquipmentBase{
				id:      id,
				name:    name,
				room_id: prop.SpaceId,
			},
			EquipmentRAC: EquipmentRAC{
				q_min: prop.Qmin,
				q_max: prop.Q_max,
				v_min: prop.V_min,
				v_max: prop.V_max,
				bf:    prop.Bf,
			},
		}
	} else if heType == "floor_heating" {
		roomID := bs.get_room_id_by_boundary_id(prop.BoundaryId)
		return &HeatingEquipmentFloorHeating{
			EquipmentBase: EquipmentBase{
				id:      id,
				name:    name,
				room_id: roomID,
			},
			EquipmentFloor: EquipmentFloor{
				boundary_id:      prop.BoundaryId,
				max_capacity:     prop.MaxCapacity,
				area:             prop.Area,
				convection_ratio: prop.ConvectionRatio,
			},
		}
	} else {
		panic("Invalid heating equipment type")
	}
}

func _create_cooling_equipment(eq *EquipmentJson, bs *Boundaries) interface{} {
	ceType := eq.EquipmentType
	id := eq.Id
	name := eq.Name
	prop := eq.Property

	if ceType == "rac" {
		return CoolingEquipmentRAC{
			EquipmentBase: EquipmentBase{
				id:      id,
				name:    name,
				room_id: prop.SpaceId,
			},
			EquipmentRAC: EquipmentRAC{
				q_min: prop.Qmin,
				q_max: prop.Q_max,
				v_min: prop.V_min,
				v_max: prop.V_max,
				bf:    prop.Bf,
			},
		}
	} else if ceType == "floor_cooling" {
		room_id := bs.get_room_id_by_boundary_id(prop.BoundaryId)
		return CoolingEquipmentFloorCooling{
			EquipmentBase: EquipmentBase{
				id:      id,
				name:    name,
				room_id: room_id,
			},
			EquipmentFloor: EquipmentFloor{
				boundary_id:      prop.BoundaryId,
				max_capacity:     prop.MaxCapacity,
				area:             prop.Area,
				convection_ratio: prop.ConvectionRatio,
			},
		}
	} else {
		panic("Invalid cooling equipment type")
	}
}

func _get_q_rs_max_is(es []interface{}, n_rm int) []float64 {
	q_rs_max_is := make([]float64, n_rm)

	for i := 0; i < len(es); i++ {
		switch es[i].(type) {
		case HeatingEquipmentFloorHeating, CoolingEquipmentFloorCooling:
			b := es[i].(EquipmentBase)
			f := es[i].(EquipmentFloor)
			q_rs_max_is[b.room_id] += f.max_capacity * f.area
		}
	}

	return q_rs_max_is
}

/*
室に放射暖冷房があるか否かを判定する。

Returns:

	放射暖冷房の有無, [i, 1]
*/
func _get_is_radiative_is(es []interface{}, n_rm int) []bool {
	is_radiative_is := make([]bool, n_rm)
	for i := 0; i < len(es); i++ {
		switch es[i].(type) {
		case HeatingEquipmentFloorHeating:
			is_radiative_is[es[i].(HeatingEquipmentFloorHeating).room_id] = true
		case CoolingEquipmentFloorCooling:
			is_radiative_is[es[i].(CoolingEquipmentFloorCooling).room_id] = true
		}
	}
	return is_radiative_is
}

// 室ごとの空調機器による対流成分の比率(-) [I]
func _get_beta_is(es []interface{}, n_rm int) []float64 {
	if len(es) == 0 {
		return make([]float64, n_rm)
	}

	// 室・空調機器ごとの対流成分比率(-) [I,Es]
	f_beta_eqp_ks_is := _get_f_beta_eqp_ks_is(es, n_rm)

	// 空調機器の最大放熱能力が設置される部屋に属する全ての機器の最大放熱能力の合計に占める割合(-) [I,Es]
	r_max_ks_is := _get_r_max_ks_is(es, n_rm)

	// 空調機器の最大放熱能力の対流成分が設置される部屋に属する全ての機器の最大放熱能力の合計に占める割合(-) [I,Es]
	var beta mat.Dense
	beta.MulElem(f_beta_eqp_ks_is, r_max_ks_is)

	// 室ごとの空調機器による対流成分の比率(-) [I]
	ret := make([]float64, n_rm)
	for i := 0; i < n_rm; i++ {
		for k := 0; k < len(es); k++ {
			ret[i] += beta.At(k, i)
		}
	}

	return ret
}

// 室と床暖房の設置をあわらす。設置された室・空調機のペアの要素にのみ1が設定されている。
// サイズ [I, Es] 疎行列である。
func _get_p_ks_is(es []interface{}, n_rm int) *mat.Dense {
	if len(es) == 0 {
		return mat.NewDense(1, n_rm, nil)
	}

	p_ks_is := mat.NewDense(len(es), n_rm, nil)

	for k, e := range es {
		switch e.(type) {
		case HeatingEquipmentFloorHeating, CoolingEquipmentFloorCooling:
			eb := e.(EquipmentBase)
			p_ks_is.Set(k, eb.room_id, 1.0)
		}
	}

	return p_ks_is
}

// f_beta_eqp[I, Es] (-): 室・機器ごとの対流成分比率
func _get_f_beta_eqp_ks_is(es []interface{}, n_rm int) *mat.Dense {
	if len(es) == 0 {
		return mat.NewDense(1, n_rm, nil)
	}
	f_beta_eqp_ks_is := mat.NewDense(len(es), n_rm, nil)
	for k, e := range es {
		switch e.(type) {
		case *HeatingEquipmentFloorHeating:
			f_beta_eqp_ks_is.Set(k, e.(*HeatingEquipmentFloorHeating).room_id, e.(*HeatingEquipmentFloorHeating).convection_ratio)
		case *CoolingEquipmentFloorCooling:
			f_beta_eqp_ks_is.Set(k, e.(*CoolingEquipmentFloorCooling).room_id, e.(*CoolingEquipmentFloorCooling).convection_ratio)
		}
	}
	return f_beta_eqp_ks_is
}

// r_max [I,Es] : 機器の最大放熱能力が設置される部屋に属する全ての機器の最大放熱能力の合計に占める割合[-]
func _get_r_max_ks_is(es []interface{}, n_rm int) *mat.Dense {
	if len(es) == 0 {
		return mat.NewDense(1, n_rm, nil)
	}
	q_max_ks_is := mat.NewDense(len(es), n_rm, nil)
	for k, e := range es {
		switch e.(type) {
		case *HeatingEquipmentFloorHeating:
			q_max_ks_is.Set(k, e.(*HeatingEquipmentFloorHeating).room_id, e.(*HeatingEquipmentFloorHeating).max_capacity*e.(*HeatingEquipmentFloorHeating).area)
		case *CoolingEquipmentFloorCooling:
			q_max_ks_is.Set(k, e.(*CoolingEquipmentFloorCooling).room_id, e.(*CoolingEquipmentFloorCooling).max_capacity*e.(*CoolingEquipmentFloorCooling).area)
		}
	}

	// 室ごとに最大放熱能力を加算
	sum_of_q_max_is := make([]float64, n_rm)
	for i := 0; i < len(es); i++ {
		for j := 0; j < n_rm; j++ {
			sum_of_q_max_is[j] += q_max_ks_is.At(i, j)
		}
	}
	for i := 0; i < len(es); i++ {
		// ゼロ除算の回避
		if sum_of_q_max_is[i] == 0.0 {
			sum_of_q_max_is[i] = 1.0
		}
	}

	// 機器の最大放熱能力が設置される部屋に属する全ての機器の最大補熱能力の合計に占める割合
	for i := 0; i < len(es); i++ {
		for j := 0; j < n_rm; j++ {
			q_max_ks_is.Set(i, j, q_max_ks_is.At(i, j)/sum_of_q_max_is[j])
		}
	}
	return q_max_ks_is
}

// 放射成分
func _get_f_flr_js_is(es []interface{}, n_rm int, n_b int) *mat.Dense {
	if len(es) == 0 {
		return mat.NewDense(n_b, n_rm, nil)
	}

	f_flr_eqp_js_ks := _get_f_flr_eqp_js_ks(es, n_b)    // 境界j,機器esごとの対流成分比率 (疎)
	f_beta_eqp_ks_is := _get_f_beta_eqp_ks_is(es, n_rm) // 部屋i,機器esごとの対流成分比率（疎）

	// 機器の最大放熱能力が設置される部屋に属する全ての機器の最大補熱能力の合計に占める割合 [-]
	r_max_ks_is := _get_r_max_ks_is(es, n_rm)

	// 室ごとの空調機器による対流成分の比率(-) [I]
	beta_is := _get_beta_is(es, n_rm)

	// 室への床暖房の設置状況 (1 or 0) [I,Es]
	p_ks_is := _get_p_ks_is(es, n_rm)

	// 部屋i,機器esごとの放射成分比率（疎）
	// 放射成分比率 := 1 - 対流成分比率
	var subResult mat.Dense
	subResult.Sub(p_ks_is, f_beta_eqp_ks_is)

	// 部屋i,機器esごとの最大放射能力 [W]（疎）
	var elementMulResult mat.Dense
	elementMulResult.MulElem(&subResult, r_max_ks_is)

	// Compute np.dot(f_flr_eqp_js_ks, elementMulResult)
	var dotResult mat.Dense
	dotResult.Mul(f_flr_eqp_js_ks, &elementMulResult)

	// Create a diagonal matrix with diagonal elements as 1/(1 - beta_is)
	betaInverse := mat.NewDense(len(beta_is), len(beta_is), nil)
	for i := 0; i < len(beta_is); i++ {
		betaInverse.Set(i, i, 1.0/(1.0-beta_is[i]))
	}

	// Multiply dotResult with betaInverse
	var finalResult mat.Dense
	finalResult.Mul(&dotResult, betaInverse)

	return &finalResult
}

// f_flr_eqp[J, E] [-]: 対流成分比率比率
func _get_f_flr_eqp_js_ks(es []interface{}, n_b int) *mat.Dense {
	if len(es) == 0 {
		return nil
	}

	f_flr_eqp_js_ks := mat.NewDense(n_b, len(es), nil)

	for k, e := range es {
		switch e.(type) {
		case HeatingEquipmentFloorHeating, CoolingEquipmentFloorCooling:
			hefh := e.(*HeatingEquipmentFloorHeating)
			f_flr_eqp_js_ks.Set(hefh.boundary_id, k, hefh.convection_ratio)
		}
	}

	return f_flr_eqp_js_ks
}

var __f_l_cl_wgt_is_is_n mat.Dense
var __f_l_cl_cst_is_n mat.VecDense

func (eq *Equipments) make_get_f_l_cl_funcs() func(*mat.VecDense, *mat.VecDense, []float64) (*mat.VecDense, *mat.Dense) {
	/*
		Args:
			l_cs_is_n: ステップ n からステップ n+1 における室 i の暖冷房設備の顕熱処理量（暖房を正・冷房を負とする）, W, [i, 1]
			theta_r_is_n_pls: ステップ n+1 における室 i の温度, degree C, [i, 1]
			x_r_ntr_is_n_pls: ステップ n+1 における室 i の加湿・除湿を行わない場合の絶対湿度, kg/kg(DA), [i, 1]

		Returns:
			タプル
				ステップ n　からステップ n+1 における係数 f_l_cl_wgt, kg/s(kg/kg(DA)), [i, i]
				ステップ n　からステップ n+1 における係数 f_l_cl_cst, kg/s, [i, 1]
	*/
	get_f_l_cl := func(
		l_cs_is_n *mat.VecDense,
		theta_r_is_n_pls *mat.VecDense,
		x_r_ntr_is_n_pls []float64,
	) (*mat.VecDense, *mat.Dense) {
		n_rm := eq._n_rm
		if __f_l_cl_wgt_is_is_n.IsEmpty() {
			__f_l_cl_wgt_is_is_n = *mat.NewDense(n_rm, n_rm, nil)
		}
		if __f_l_cl_cst_is_n.IsEmpty() {
			__f_l_cl_cst_is_n = *mat.NewVecDense(n_rm, nil)
		}

		f_l_cl_wgt_is_is_n := &__f_l_cl_wgt_is_is_n
		f_l_cl_cst_is_n := &__f_l_cl_cst_is_n

		for _, ce := range eq._ces {
			ls_a, ls_b := eq._get_ls_a_ls_b(
				l_cs_is_n,
				theta_r_is_n_pls,
				x_r_ntr_is_n_pls,
				ce.(interface{}),
			)

			// 係数 la 及び lb それぞれ合計する。
			// la [i,i] kg/s(kg/kg(DA))
			// lb [i,1] kg/kg(DA)
			// TODO: La は正負が仕様書と逆になっている
			f_l_cl_wgt_is_is_n.Sub(f_l_cl_wgt_is_is_n, ls_a)
			f_l_cl_cst_is_n.AddVec(f_l_cl_cst_is_n, ls_b)
		}

		return f_l_cl_cst_is_n, f_l_cl_wgt_is_is_n
	}

	return get_f_l_cl
}

func (eq *Equipments) _get_ls_a_ls_b(
	l_cs_is_n mat.Vector,
	theta_r_is_n_pls mat.Vector,
	x_r_ntr_is_n_pls []float64,
	ce interface{},
) (mat.Matrix, mat.Vector) {

	switch ce.(type) {
	case CoolingEquipmentRAC:
		return eq._func_rac(
			l_cs_is_n,
			theta_r_is_n_pls,
			x_r_ntr_is_n_pls,
			ce.(CoolingEquipmentRAC),
		)
	case CoolingEquipmentFloorCooling:
		panic("NotImplementedError")
	default:
		panic("NotImplementedError")
	}
}

var __brmx_is_is *mat.Dense
var __brxc_is *mat.VecDense

func (eq *Equipments) _func_rac(
	l_cs_is_n mat.Vector,
	theta_r_is_n_pls mat.Vector,
	x_r_ntr_is_n_pls []float64,
	ce CoolingEquipmentRAC,
) (mat.Matrix, mat.Vector) {

	// 室の数
	n_rm := eq._n_rm

	// Lcsは加熱が正で表される。
	// 加熱時は除湿しない。
	// 以下の取り扱いを簡単にするため（冷房負荷を正とするため）、正負を反転させる
	q_s_i_n := -l_cs_is_n.AtVec(ce.room_id)
	theta_r_i_n_pls := theta_r_is_n_pls.AtVec(ce.room_id)
	x_r_ntr_i_n_pls := x_r_ntr_is_n_pls[ce.room_id]

	v_rac_i_n := _get_vac_rac_i_n(
		ce.q_max,
		ce.q_min,
		q_s_i_n,
		ce.v_max/60,
		ce.v_min/60,
	)

	theta_rac_ex_srf_i_n_pls := _get_theta_rac_ex_srf_i_n_pls(
		ce.bf,
		q_s_i_n,
		theta_r_i_n_pls,
		v_rac_i_n,
	)

	x_rac_ex_srf_i_n_pls := _get_x_rac_ex_srf_i_n_pls(theta_rac_ex_srf_i_n_pls)

	var brmx_rac_is float64
	if x_r_ntr_i_n_pls > x_rac_ex_srf_i_n_pls && q_s_i_n > 0.0 {
		brmx_rac_is = rho_a * v_rac_i_n * (1 - ce.bf)
	} else {
		brmx_rac_is = 0.0
	}

	var brcx_rac_is float64
	if x_r_ntr_i_n_pls > x_rac_ex_srf_i_n_pls && q_s_i_n > 0.0 {
		brcx_rac_is = rho_a * v_rac_i_n * (1 - ce.bf) * x_rac_ex_srf_i_n_pls
	} else {
		brcx_rac_is = 0.0
	}

	if __brmx_is_is == nil {
		__brmx_is_is = mat.NewDense(n_rm, n_rm, nil)
	}
	if __brxc_is == nil {
		__brxc_is = mat.NewVecDense(n_rm, nil)
	}

	brmx_is_is := __brmx_is_is
	brxc_is := __brxc_is

	brmx_is_is.Set(ce.room_id, ce.room_id, brmx_rac_is)
	brxc_is.SetVec(ce.room_id, brcx_rac_is)

	return brmx_is_is, brxc_is
}

/*
ルームエアコンディショナーの室内機の熱交換器表面の絶対湿度を求める。
Args:

	theta_rac_ex_srf_i_n_pls: ステップ n+1 における室 i に設置されたルームエアコンディショナーの室内機の熱交換器表面温度,degree C

Returns:

	ステップ n+1 における室 i に設置されたルームエアコンディショナーの室内機の熱交換器表面の絶対湿度, kg/kg(DA)

Notes:

	繰り返し計算（温度と湿度） eq.12
*/
func _get_x_rac_ex_srf_i_n_pls(theta_rac_ex_srf_i_n_pls float64) float64 {
	return get_x(get_p_vs(theta_rac_ex_srf_i_n_pls))
}

/*
ステップ n+1 における室 i に設置されたルームエアコンディショナーの室内機の熱交換器表面温度を計算する。
Args:

	bf_rac_i: 室 i に設置されたルームエアコンディショナーの室内機の熱交換器のバイパスファクター, -
	q_s_i_n: ステップ n から n+1 における室 i の顕熱負荷, W
	theta_r_i_n_pls: ステップ n+1 における室 i の温度, degree C
	v_rac_i_n: ステップ n から n+1 における室 i に設置されたルームエアコンディショナーの吹き出し風量, m3/s

Returns:

	ステップ n+1 における室 i に設置されたルームエアコンディショナーの室内機の熱交換器表面温度, degree C

Notes:

	繰り返し計算（温度と湿度） eq.14
*/
func _get_theta_rac_ex_srf_i_n_pls(
	bf_rac_i float64,
	q_s_i_n float64,
	theta_r_i_n_pls float64,
	v_rac_i_n float64,
) float64 {
	return theta_r_i_n_pls - q_s_i_n/(c_a*rho_a*v_rac_i_n*(1.0-bf_rac_i))
}

/*
ルームエアコンディショナーの吹き出し風量を顕熱負荷に応じて計算する。

Args:

	q_rac_max_i: 室 i に設置されたルームエアコンディショナーの最大能力, W
	q_rac_min_i: 室 i に設置されたルームエアコンディショナーの最小能力, W
	q_s_i_n:　ステップ n からステップ n+1 における室 i の顕熱負荷, W
	v_rac_max_i: 室 i に設置されたルームエアコンディショナーの最小能力時における風量, m3/s
	v_rac_min_i: 室 i に設置されたルームエアコンディショナーの最大能力時における風量, m3/s

Returns:

	室iに設置されたルームエアコンディショナーの吹き出し風量, m3/s

Notes:

	繰り返し計算（湿度と潜熱） eq.14
*/
func _get_vac_rac_i_n(
	q_rac_max_i float64,
	q_rac_min_i float64,
	q_s_i_n float64,
	v_rac_max_i float64,
	v_rac_min_i float64,
) float64 {

	// 吹き出し風量（仮）, m3/s
	v := v_rac_min_i*(q_rac_max_i-q_s_i_n)/(q_rac_max_i-q_rac_min_i) +
		v_rac_max_i*(q_rac_min_i-q_s_i_n)/(q_rac_min_i-q_rac_max_i)

	// 下限値・上限値でクリップ後の吹き出し風量, m3/s
	var v_rac_i_n float64
	if v < v_rac_min_i {
		v_rac_i_n = v_rac_min_i
	} else if v > v_rac_max_i {
		v_rac_i_n = v_rac_max_i
	} else {
		v_rac_i_n = v
	}

	return v_rac_i_n
}
