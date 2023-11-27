package heat_load_calc

// **** 日射吸収量 ****
// See: https://github.com/HC-Energy/hc-energy-doc-ja/blob/main/docs/source/contents/03_26_eval_solar_absorption.rst

import (
	"gonum.org/v1/gonum/mat"
)

/*
	設置された家具による透過日射吸収熱量の時間平均値を求める。

	Args:
		q_trs_sor_is_ns:

	Returns:
		ステップ n からステップ n+1 における室 i に設置された家具による透過日射吸収熱量時間平均値, W, [i, n]
	Notes:
		式(1)
*/
func get_q_sol_frt_is_ns(q_trs_sor_is_ns mat.Matrix) *mat.Dense {
	var q_sol_frt_is_ns mat.Dense

	q_sol_frt_is_ns.Scale(r_sol_frt, q_trs_sor_is_ns)

	return &q_sol_frt_is_ns
}

/*
	境界jの透過日射吸収熱量を求める。

	Args:
		p_is_js: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [i, j]
		a_s_js: 境界 j の面積, m2, [j, 1]
		p_s_sol_abs_js: 境界 j において透過日射を吸収するか否かを表す係数（吸収する場合は 1 とし、吸収しない場合は 0 とする。）, -, [j, 1]
		p_js_is: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [j, i]
		q_trs_sol_is_ns:　ステップ n における室 i の透過日射熱量, W, [i, 1]

	Returns:
       ステップ n における境界 j の透過日射吸収熱量, W/m2, [j, n]

	Notes:
		式(2)
*/
func get_q_s_sol_js_ns(
	p_is_js mat.Matrix,
	a_s_js mat.Vector,
	p_s_sol_abs_js mat.Vector,
	p_js_is mat.Matrix,
	q_trs_sol_is_ns mat.Matrix,
) *mat.Dense {
	// 室iにおける日射が吸収される境界の面積の合計, m2, [i, 1]
	var a_s_js_times_p_s_sol_abs_js mat.VecDense
	a_s_js_times_p_s_sol_abs_js.MulElemVec(a_s_js, p_s_sol_abs_js)
	var a_srf_abs_is mat.VecDense
	a_srf_abs_is.MulVec(p_is_js, &a_s_js_times_p_s_sol_abs_js)

	// ステップnの境界jにおける透過日射吸収熱量, W/m2, [j, n]

	// q_trs_sol_is_ns / a_srf_abs_is
	var q_trs_sol_div_a_srf_abs_is mat.Dense
	q_trs_sol_div_a_srf_abs_is.Apply(func(i, j int, v float64) float64 {
		return v / a_srf_abs_is.AtVec(i)
	}, q_trs_sol_is_ns)

	// np.dot(p_js_is, q_trs_sol_is_ns / a_srf_abs_is) * p_s_sol_abs_js * (1.0 - _get_r_sol_frt())
	var temp_result mat.Dense
	temp_result.Mul(p_js_is, &q_trs_sol_div_a_srf_abs_is)
	temp_result.Apply(func(i, j int, v float64) float64 {
		return v * p_s_sol_abs_js.AtVec(i)
	}, &temp_result)
	temp_result.Scale(1.0-r_sol_frt, &temp_result)

	return &temp_result
}

/*
	室の日射が吸収される境界の面積の合計を取得する。

	Args:
        p_is_js: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [i, j]
        a_s_js: 境界 j の面積, m2, [j, 1]
        p_s_sol_abs_js: 境界 j において透過日射を吸収するか否かを表す係数（吸収する場合は 1 とし、吸収しない場合は 0 とする。）, -, [j, 1]
    Returns:
        室iの日射が吸収される境界の面積の合計, m2, [i, 1]
    Note:
        式(3)
*/
func _get_a_s_abs_is(p_is_js mat.Matrix, a_s_js mat.Matrix, p_s_sol_abs_js mat.Matrix) mat.Matrix {
	var tmp mat.Dense
	tmp.MulElem(p_is_js, a_s_js)
	var res mat.Dense
	res.Mul(p_is_js, &tmp)
	return &res
}

// 室内侵入日射のうち家具に吸収される割合, -
const r_sol_frt = 0.5
