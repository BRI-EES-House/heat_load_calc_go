package heat_load_calc

import (
	"gonum.org/v1/gonum/mat"
)

var __next_temp_and_load__kc *mat.Dense
var __next_temp_and_load__nt *mat.VecDense
var __next_temp_and_load__theta_set *mat.VecDense
var __next_temp_and_load__c *mat.VecDense
var __next_temp_and_load__lr_set *mat.VecDense
var __next_temp_and_load__lc_set *mat.VecDense
var __next_temp_and_load__r *mat.VecDense

func get_next_temp_and_load(
	ac_demand_is_ns *ScheduleData,
	brc_ot_is_n []float64,
	brm_ot_is_is_n *mat.Dense,
	brl_ot_is_is_n *mat.Dense,
	theta_lower_target_is_n *mat.VecDense,
	theta_upper_target_is_n *mat.VecDense,
	operation_mode_is_n []OperationMode,
	is_radiative_heating_is []bool,
	is_radiative_cooling_is []bool,
	lr_h_max_cap_is []float64,
	lr_cs_max_cap_is []float64,
	theta_natural_is_n []float64,
	n int,
) (*mat.VecDense, *mat.VecDense, *mat.VecDense) {

	var nn int
	if n < 0 {
		c := ac_demand_is_ns.Len()
		nn = n + c
	} else {
		nn = n
	}
	ac_demand_is_n := ac_demand_is_ns.Get(nn)

	roomShape := len(operation_mode_is_n)

	// 室の数
	n_room := roomShape

	// 係数　kt, W / K, [i, i], float型
	kt := brm_ot_is_is_n

	// 係数 kc, [i, i], float型
	if __next_temp_and_load__kc == nil {
		__next_temp_and_load__kc = mat.NewDense(n_room, n_room, nil)
		for i := 0; i < n_room; i++ {
			__next_temp_and_load__kc.Set(i, i, 1)
		}
	}
	kc := __next_temp_and_load__kc

	// 係数 kr, [i, i], float型
	kr := brl_ot_is_is_n

	// 係数 k, W, [i, 1], float型
	k := brc_ot_is_n

	// 室温指定を表す係数 theta_set, [i, 1], int型
	// 指定する = 0, 指定しない = 1
	// 室温を指定しない場合は、 operation_mode が STOP_CLOSE or STOP_OPEN の場合である。
	// 後で再計算する際に、負荷が機器容量を超えている場合は、最大暖房／冷房負荷で処理されることになるため、
	// 室温を指定しない場合は、この限りではない。
	// ---
	// nt:
	// nt = 0 （室温を指定する） に対応する要素に、ターゲットとなるOTを代入する。
	// nt = 1 （室温を指定しない）場合は、theta_set は 0 にしなければならない。
	// ---
	// 対流空調指定を表す係数 c , [i, 1], int型
	// 指定する = 0, 指定しない = 1
	// 対流空調を指定しない場合は、対流空調をしている場合に相当するので、
	//   operation_mode が　HEATING でかつ、 is_radiative_heating_is が false の場合か、
	//   operation_mode が COOLING でかつ、 is_radiative_cooling_is が false の場合
	// のどちらかである。
	// ---
	// 放射空調指定を表す係数 r, [i, 1], int型
	// 指定する = 0, 指定しない = 1
	// 放射空調を指定しない場合は、放射空調をしている場合に相当するので、
	//   operation_mode が　HEATING でかつ、 is_radiative_heating_is が true の場合か、
	//   operation_mode が COOLING でかつ、 is_radiative_cooling_is が true の場合
	// のどちらかである。
	if __next_temp_and_load__nt == nil {
		__next_temp_and_load__nt = mat.NewVecDense(roomShape, nil)
	}
	if __next_temp_and_load__theta_set == nil {
		__next_temp_and_load__theta_set = mat.NewVecDense(roomShape, nil)
	}
	if __next_temp_and_load__c == nil {
		__next_temp_and_load__c = mat.NewVecDense(roomShape, nil)
	}
	if __next_temp_and_load__r == nil {
		__next_temp_and_load__r = mat.NewVecDense(roomShape, nil)
	}
	r := __next_temp_and_load__r
	c := __next_temp_and_load__c
	nt := __next_temp_and_load__nt
	theta_set := __next_temp_and_load__theta_set
	for i := 0; i < roomShape; i++ {
		if operation_mode_is_n[i] == HEATING &&
			theta_natural_is_n[i] < theta_lower_target_is_n.AtVec(i) {
			// 室温指定係数
			nt.SetVec(i, 0)

			// 室温指定
			val := theta_lower_target_is_n.AtVec(i)*ac_demand_is_n.AtVec(i) + theta_natural_is_n[i]*(1.0-ac_demand_is_n.AtVec(i))
			theta_set.SetVec(i, val)

			// 対流空調指定 c / 放射空調指定 r
			if is_radiative_heating_is[i] {
				c.SetVec(i, 0)
				r.SetVec(i, 1)
			} else {
				c.SetVec(i, 1)
				r.SetVec(i, 0)
			}
		} else if operation_mode_is_n[i] == COOLING &&
			theta_upper_target_is_n.AtVec(i) < theta_natural_is_n[i] {

			// 室温指定係数
			nt.SetVec(i, 0)

			// 室温指定
			val := theta_upper_target_is_n.AtVec(i)*ac_demand_is_n.AtVec(i) + theta_natural_is_n[i]*(1.0-ac_demand_is_n.AtVec(i))
			theta_set.SetVec(i, val)

			// 対流空調指定 c / 放射空調指定 r
			if is_radiative_cooling_is[i] {
				c.SetVec(i, 0)
				r.SetVec(i, 1)
			} else {
				c.SetVec(i, 1)
				r.SetVec(i, 0)
			}
		} else {
			// 室温指定係数
			nt.SetVec(i, 1) // 室温を指定しない
			// 室温指定
			theta_set.SetVec(i, 0.0)
			// 対流空調指定 c
			c.SetVec(i, 0)
			// 放射空調指定
			r.SetVec(i, 0)
		}
	}

	if __next_temp_and_load__lc_set == nil {
		__next_temp_and_load__lc_set = mat.NewVecDense(roomShape, nil)
	}
	lc_set := __next_temp_and_load__lc_set
	lc_set.Zero()

	// r = 0 （放射空調を指定する）に対応する要素に、0.0 を代入する。
	// 放射空調を指定する場合は空調をしていないことに相当するため。ただし、後述する、最大能力で動く場合は、その値を代入することになる。
	// 放射空調を指定しない場合は、 lr_set には 0.0 を入れなければならない。
	// 結果として、ここでは、あらゆるケースで 0.0 が代入される。
	if __next_temp_and_load__lr_set == nil {
		__next_temp_and_load__lr_set = mat.NewVecDense(roomShape, nil)
	}
	lr_set := __next_temp_and_load__lr_set
	lr_set.Zero()

	// theta 温度, degree C, [i, 1]
	// lc 対流空調負荷, W, [i, 1]
	// lr 放射空調負荷, W, [i, 1]
	theta, lc, lr := get_load_and_temp(kt, kc, kr, k, nt, theta_set, c, lc_set, r, lr_set)

	// 計算された放射空調負荷が最大放熱量を上回る場合は、放熱量を最大放熱量に固定して、対流空調負荷を未知数として再計算する。
	over_lr := make([]bool, roomShape)
	for i := range over_lr {
		over_lr[i] = lr.AtVec(i) > lr_h_max_cap_is[i]
	}

	for i, v := range over_lr {
		if v {
			// 対流負荷を未知数とする。
			c.SetVec(i, 1)

			// 放射負荷を最大放熱量に指定する。
			r.SetVec(i, 0)
			lr_set.SetVec(i, lr_h_max_cap_is[i])
		}
	}

	// 計算された放射空調負荷が最大放熱量を下回る場合は、放熱量を最大放熱量に固定して、対流空調負荷を未知数として再計算する。
	// 注意：冷房の最大放熱量は正の値で指定される。一方、計算される負荷（lr）は、冷房の場合、負の値で指定される。
	under_lr := make([]bool, roomShape)
	for i := range under_lr {
		under_lr[i] = lr.AtVec(i) < -lr_cs_max_cap_is[i]
	}

	for i, v := range under_lr {
		if v {
			// 対流負荷を未知数とする。
			c.SetVec(i, 1)

			// 放射負荷を最大放熱量に指定する。
			r.SetVec(i, 0)
			lr_set.SetVec(i, -lr_cs_max_cap_is[i])
		}
	}

	// 放射暖房を最大放熱量に指定して再計算する。
	theta, lc, lr = get_load_and_temp(kt, kc, kr, k, nt, theta_set, c, lc_set, r, lr_set)

	return theta, lc, lr
}

var __m_nt *mat.Dense
var __m_c *mat.Dense
var __m_r *mat.Dense
var __load_and_temp__x1 mat.Dense
var __load_and_temp__x1_term1 mat.Dense
var __load_and_temp__x1_term2 mat.Dense
var __load_and_temp__x1_term3 mat.Dense
var __load_and_temp__x2 mat.VecDense
var __load_and_temp__x2_term1 mat.VecDense
var __load_and_temp__x2_term2 mat.VecDense
var __load_and_temp__x2_term3 mat.VecDense
var __load_and_temp__v mat.VecDense
var __load_and_temp__theta_rq mat.VecDense
var __load_and_temp__lc_rq mat.VecDense
var __load_and_temp__lr_rq mat.VecDense

/*
Args:

	kt: 係数 kt, W/K, [i, i], float型
	kc: 係数 kc, [i, i], float型
	kr: 係数 kr, [i, i], float型
	k: 係数 k, W, [i, 1], float型
	nt: 室温指定を表す係数, [i, 1], int型
	theta_set: 室温を指定する場合の室温, degree C, [i, 1], float型
	c: 対流空調指定を表す係数, [i, 1], int型
	lc_set: 対流空調を指定する場合の放熱量, W, [i, 1], float型
	r: 放射空調指定を表す係数, [i, 1], int型
	lr_set: 放射空調を指定する場合の放熱量, W, [i, 1], float型

Returns:

	温度, degree C, [i, 1]
	対流空調負荷, W, [i, 1]
	放射空調負荷, W, [i, 1]

Notes:

	各係数によって、
	kt theta = kc Lc + kr Lr + k
	が維持される。

	室温・対流空調・放射空調 を指定する場合は、指定しない = 1, 指定する = 0 とする。
	theta_set, lc_set, lr_set について、指定しない場合は必ず 0.0 とする。
*/
func get_load_and_temp(
	kt mat.Matrix,
	kc mat.Matrix,
	kr *mat.Dense,
	k []float64,
	nt *mat.VecDense,
	theta_set mat.MutableVector,
	c *mat.VecDense,
	lc_set mat.MutableVector,
	r *mat.VecDense,
	lr_set mat.MutableVector,
) (*mat.VecDense, *mat.VecDense, *mat.VecDense) {

	n := nt.Len()

	// Modify theta_set where nt == 1
	for i := 0; i < n; i++ {
		// 室温を指定しない場合 nt = 1
		// 室温を指定しない場合 theta_set は 0.0 とする。
		if nt.AtVec(i) == 1 {
			theta_set.SetVec(i, 0.0)
		}

		// 対流空調負荷を指定しない場合 c = 1
		// 対流空調負荷を指定しない場合 lc_set = 0.0 とする。
		if c.AtVec(i) == 1 {
			lc_set.SetVec(i, 0.0)
		}

		// 放射空調負荷を指定しない場合 r = 1
		// 放射空調負荷を指定しない場合 lr_set = 0.0 とする。
		if r.AtVec(i) == 1 {
			lr_set.SetVec(i, 0.0)
		}
	}

	// nt, c, r matrix
	if __m_nt == nil || __m_c == nil || __m_r == nil {
		__m_nt, __m_c, __m_r = mat.NewDense(n, n, nil), mat.NewDense(n, n, nil), mat.NewDense(n, n, nil)
	}
	m_nt, m_c, m_r := __m_nt, __m_c, __m_r
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			m_nt.Set(i, j, nt.AtVec(j))
			m_c.Set(i, j, c.AtVec(j))
			m_r.Set(i, j, r.AtVec(j))
		}
	}

	// Calculate x1 and x2

	// kt * nt.T - kc * c.T - kr * r.T
	x1, x1_term1, x1_term2, x1_term3 := &__load_and_temp__x1, &__load_and_temp__x1_term1, &__load_and_temp__x1_term2, &__load_and_temp__x1_term3
	x1_term1.MulElem(kt, m_nt)
	x1_term2.MulElem(kc, m_c)
	x1.Sub(x1_term1, x1_term2)
	if kr != nil {
		x1_term3.MulElem(kr, m_r)
		x1.Sub(x1, x1_term3)
	}

	// -np.dot(kt, theta_set) + np.dot(kc, lc_set) + np.dot(kr, lr_set) + k
	x2, x2_term1, x2_term2, x2_term3 := &__load_and_temp__x2, &__load_and_temp__x2_term1, &__load_and_temp__x2_term2, &__load_and_temp__x2_term3
	x2_term1.MulVec(kt, theta_set)
	x2_term2.MulVec(kc, lc_set)
	x2.SubVec(x2_term2, x2_term1)
	if kr != nil {
		x2_term3.MulVec(kr, lr_set)
		x2.AddVec(x2, x2_term3)
	}
	x2.AddVec(x2, mat.NewVecDense(len(k), k))

	v := &__load_and_temp__v
	v.SolveVec(x1, x2)

	// 求めるべき数値
	// nt, c, r それぞれ、1の場合（値を指定しない場合）は、vで表される値が入る。
	// 反対に、 0 の場合（値を指定する場合）、は、それぞれ、theta_set, lc_set, lr_set の値が入る。
	theta_rq, lc_rq, lr_rq := &__load_and_temp__theta_rq, &__load_and_temp__lc_rq, &__load_and_temp__lr_rq
	__combineVectorTo(theta_rq, n, v, theta_set, nt)
	__combineVectorTo(lc_rq, n, v, lc_set, c)
	__combineVectorTo(lr_rq, n, v, lr_set, r)

	return theta_rq, lc_rq, lr_rq
}

// 長さnのベクトルA, ベクトルBをベクトルrの比率で按分する。rは0-1の値を取り、結果に含まれるAの比率を示す
// 混合結果 X := A * r + B * (1-r) である。
func __combineVectorTo(x *mat.VecDense, n int, a mat.Vector, b mat.Vector, r mat.Vector) {
	// X := A * r + B * (1-r)  を式変形する。
	//   := A * r + B - B * r
	//   := A * r - B * r + B
	//   := (A - B) * r  + B

	x.SubVec(a, b)     // A-B
	x.MulElemVec(x, r) //(A-B)*r
	x.AddVec(x, b)     // (A-B)*r + B
}
