package heat_load_calc

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/floats"
)

/*
PMVを計算する。
    人体周りの熱伝達率を着衣温度との収束計算をした上でPMVを計算する。

    Args:
        p_a_is_n:　ステップnにおける室iの水蒸気圧, Pa, [i, 1]
        theta_r_is_n: ステップ n における室 i の空気温度, degree C, [i, 1]
        theta_mrt_is_n: ステップ n における室 i の在室者の平均放射温度, degree C, [i, 1]
        clo_is_n: ステップ n における室 i のClo値, [i, 1]
        v_hum_is_n: ステップ n における室 i の在室者周りの風速, m/s, [i, 1]
        met_is: 室 i の在室者のMet値, [i, 1]
        method: PMV計算時の熱伝達率計算に収束計算を行うか固定値を使用するか
    Returns:
        ステップ n における室 i の在室者のPMV, [i, 1]
*/
func get_pmv_is_n(
	p_a_is_n []float64,
	theta_r_is_n []float64,
	theta_mrt_is_n []float64,
	clo_is_n []float64,
	v_hum_is_n []float64,
	met_is []float64,
	method string,
) []float64 {
	// (1) ステップ n における室 i の在室者周りの対流熱伝達率, W/m2K, [i, 1]
	// (2) ステップ n における室 i の在室者周りの放射熱伝達率, W/m2K, [i, 1]
	// (3) ステップ n における室 i の在室者周りの総合熱伝達率, W/m2K, [i, 1]
	h_hum_c_is_n, h_hum_r_is_n, h_hum_is_n := get_h_hum(
		theta_mrt_is_n,
		theta_r_is_n,
		clo_is_n,
		v_hum_is_n,
		method,
		met_is,
	)

	c := len(theta_r_is_n)
	pmv_is_n := make([]float64, c)

	for i := 0; i < c; i++ {
		// ステップnにおける室iの在室者の作用温度, degree C, [i, 1]
		theta_ot := (h_hum_r_is_n[i]*theta_mrt_is_n[i] + h_hum_c_is_n[i]*theta_r_is_n[i]) / h_hum_is_n[i]

		// ステップ n における室 i の在室者の着衣抵抗, m2K/W, [i, 1]
		i_cl := _get_i_cl_is_n(clo_is_n[i])

		// ステップ n における室 i の在室者の着衣面積率, [i, 1]
		f_cl := _get_f_cl_is_n(i_cl)

		// ステップ n における室 i の在室者の代謝量, W/m2, [i, 1]
		m_is := _get_m_is(met_is[i])

		// ステップnにおける室iの在室者の厚着時のPMV, [i, 1]
		pmv_is_n[i] = _get_pmv_is_n(theta_r_is_n[i], p_a_is_n[i], h_hum_is_n[i], theta_ot, i_cl, f_cl, m_is)
	}

	return pmv_is_n
}

/*
PPDを計算する。

    Args:
        pmv_is_n: ステップ n における室 i の在室者のPMV, [i, 1]

    Returns:
        ステップ n における室 i の在室者のPPD, [i, 1]
*/
func get_ppd_is_n(pmv_is_n []float64) []float64 {
	ppd_is_n := make([]float64, len(pmv_is_n))
	for i, pmv := range pmv_is_n {
		pmv2 := pmv * pmv
		pmv4 := pmv2 * pmv2
		ppd_is_n[i] = 100.0 - 95.0*math.Exp(-0.03353*pmv4-0.2179*pmv2)
	}
	return ppd_is_n
}

/*
指定したPMVを満たすOTを計算する

Args:
	clo_is_n: ステップnにおける室iの在室者のClo値, [i, 1]
	p_a_is_n: ステップnにおける室iの水蒸気圧, Pa, [i, 1]
	h_hum_is_n: ステップnにおける室iの在室者周りの総合熱伝達率, W/m2K, [i, 1]
	met_is: 室 i における在室者のMet, [i, 1]
	pmv_target_is_n: ステップnにおける室iの在室者の目標PMV, [i, 1]

Returns:
	ステップnにおける室iの在室者の目標OT, [i, 1]

Notes:
	ISOで定める計算方法ではなく、前の時刻に求めた人体周りの熱伝達率、着衣温度を使用して収束計算が生じないようにしている。
*/
func get_theta_ot_target(
	clo_is_n float64,
	p_a_is_n float64,
	h_hum_is_n float64,
	met_is float64,
	pmv_target_is_n float64,
) float64 {
	// ステップ n における室 i の在室者の着衣抵抗, m2K/W, [i, 1]
	i_cl_is_n := _get_i_cl_is_n(clo_is_n)

	// 室 i の在室者の代謝量（人体内部発熱量）, W/m2
	m_is := _get_m_is(met_is)

	// ステップ n における室 i の着衣面積率, [i, 1]
	f_cl_is_n := _get_f_cl_is_n(i_cl_is_n)

	return _get_theta_ot_target_is_n(p_a_is_n, h_hum_is_n, pmv_target_is_n, i_cl_is_n, m_is, f_cl_is_n)
}

/*
在室者周りの熱伝達率を計算する。

Args:
	theta_mrt_is_n: ステップnにおける室iの在室者の平均放射温度, degree C, [i, 1]
	theta_r_is_n: ステップnにおける室iの空気温度, degree C, [i, 1]
	clo_is_n: CLO値, [i, 1]
	v_hum_is_n: ステップnにおける室iの在室者周りの風速, m/s, [i, 1]
	method: 在室者周りの熱伝達率を求める際に収束計算を行うかどうか
	met_is: 室 i の在室者のMet値, [i, 1]
Returns:
	(1) ステップ n における室 i の在室者周りの対流熱伝達率, W/m2K, [i, 1]
	(2) ステップ n における室 i の在室者周りの放射熱伝達率, W/m2K, [i, 1]
	(3) ステップ n における室 i の在室者周りの総合熱伝達率, W/m2K, [i, 1]
*/
func get_h_hum(
	theta_mrt_is_n []float64,
	theta_r_is_n []float64,
	clo_is_n []float64,
	v_hum_is_n []float64,
	method string,
	met_is []float64,
) ([]float64, []float64, []float64) {
	h_hum_c_is_n, h_hum_r_is_n := _get_h_hum_c_is_n_and_h_hum_r_is_n(
		theta_mrt_is_n,
		theta_r_is_n,
		clo_is_n,
		v_hum_is_n,
		method,
		met_is,
	)

	h_hum_is_n := _get_h_hum_is_n(h_hum_c_is_n, h_hum_r_is_n)

	return h_hum_c_is_n, h_hum_r_is_n, h_hum_is_n
}

/*
在室者周りの総合熱伝達率を計算する。
calculate the integrated heat transfere coefficient around the ocupant.
Args:
	h_hum_c_is_n: ステップ n における室 i の在室者周りの対流熱伝達率, W/m2K, [i, 1]
	h_hum_r_is_n: ステップ n における室 i の在室者周りの放射熱伝達率, W/m2K, [i, 1]
Returns:
	ステップ n における室 i の在室者周りの総合熱伝達率, W/m2K, [i, 1]
Note:
	eq.(4)
*/
func _get_h_hum_is_n(h_hum_c_is_n []float64, h_hum_r_is_n []float64) []float64 {
	h_hum_is_n := make([]float64, len(h_hum_c_is_n))
	floats.AddTo(h_hum_is_n, h_hum_c_is_n, h_hum_r_is_n)
	return h_hum_is_n
}

/*
在室者周りの熱伝達率を計算する。

    Args:
        theta_mrt_is_n: ステップnにおける室iの在室者の平均放射温度, degree C, [i, 1]
        theta_r_is_n: ステップnにおける室iの空気温度, degree C, [i, 1]
        clo_is_n: CLO値, [i, 1]
        v_hum_is_n: ステップnにおける室iの在室者周りの風速, m/s, [i, 1]
        method: 在室者周りの熱伝達率を求める際に収束計算を行うかどうか
        met_is: 室 i の在室者のMet値, [i, 1]
    Returns:
        (1) ステップ n における室 i の在室者周りの対流熱伝達率, W/m2K, [i, 1]
        (2) ステップ n における室 i の在室者周りの放射熱伝達率, W/m2K, [i, 1]
        (3) ステップ n における室 i の在室者周りの総合熱伝達率, W/m2K, [i, 1]
*/
func _get_h_hum_c_is_n_and_h_hum_r_is_n(
	theta_mrt_is_n []float64,
	theta_r_is_n []float64,
	clo_is_n []float64,
	v_hum_is_n []float64,
	method string,
	met_is []float64,
) ([]float64, []float64) {
	//NOTE: fsolveはFortranのMINPACK相当だが、Goでは相当の実装が見当たらない
	//ref: https://github.com/fortran-lang/minpack

	l := len(theta_r_is_n)
	h_hum_c_is_n := make([]float64, l)
	h_hum_r_is_n := make([]float64, l)

	if method == "convergence" {
		// 収束計算による方法

		var theta_r, v_hum, theta_mrt, clo, m float64

		f := func(t float64) float64 {
			// ステップnにおける室iの在室者周りの対流熱伝達率, W/m2K, [i, 1]
			h_hum_c := _get_h_hum_c_is_n_convergence(theta_r, t, v_hum)

			// ステップnにおける室iの在室者周りの放射熱伝達率, W/m2K, [i, 1]
			h_hum_r := _get_h_hum_r_is_n_convergence(t, theta_mrt)

			// ステップnにおける室iの在室者の作用温度, degree C, [i, 1]
			theta_ot := _get_theta_ot_is_n(h_hum_r, theta_mrt, h_hum_c, theta_r)

			theta_cl := _get_theta_cl_is_n(clo, theta_ot, m, h_hum_r, h_hum_c)

			return theta_cl - t
		}

		// Bisection method
		findRoot := func(a float64, b float64, tol float64, maxIter int) (float64, error) {
			if f(a)*f(b) >= 0 {
				return 0, fmt.Errorf("no root found in the interval [%f, %f] (f(%f)=%f, f(%f)=%f)", a, b, a, f(a), b, f(b))
			}

			var c float64
			for i := 0; i < maxIter; i++ {
				c = (a + b) / 2

				if f(c) == 0 || (b-a)/2 < tol {
					return c, nil
				}

				if f(c)*f(a) < 0 {
					b = c
				} else {
					a = c
				}
			}
			return 0, fmt.Errorf("failed to find root within %d iterations", maxIter)
		}

		a := 0.0
		b := 100.0
		tol := 1e-6
		maxIter := 100

		for i := 0; i < l; i++ {
			theta_r = theta_r_is_n[i]
			v_hum = v_hum_is_n[i]
			theta_mrt = theta_mrt_is_n[i]
			clo = clo_is_n[i]

			// 室 i の在室者の代謝量（人体内部発熱量）, W/m2
			m = met_is[i]

			// ステップnにおける室iの在室者の着衣温度, degree C, [i, 1]
			theta_cl, err := findRoot(a, b, tol, maxIter)
			if err != nil {
				panic(err)
			}

			// ステップnにおける室iの在室者周りの対流熱伝達率, W/m2K, [i, 1]
			h_hum_c_is_n[i] = _get_h_hum_c_is_n_convergence(theta_r, theta_cl, v_hum)

			// ステップnにおける室iの在室者周りの放射熱伝達率, W/m2K, [i, 1]
			h_hum_r_is_n[i] = _get_h_hum_r_is_n_convergence(theta_cl, theta_mrt)
		}

	} else if method == "constant" {
		// 収束計算によらない方法
		for i := 0; i < l; i++ {
			// ステップnにおける室iの在室者の着衣温度, degree C, [i, 1]
			theta_r := theta_r_is_n[i]

			// ステップnにおける室iの在室者周りの対流熱伝達率, W/m2K, [i, 1]
			h_hum_c_is_n[i] = _get_h_hum_c_is_n_constant(theta_r)

			// ステップnにおける室iの在室者周りの放射熱伝達率, W/m2K, [i, 1]
			h_hum_r_is_n[i] = _get_h_hum_r_is_n_constant(theta_r)
		}
	} else {
		panic("Invalid method")
	}

	return h_hum_c_is_n, h_hum_r_is_n
}

/*
PMVを計算する
Calculate the PMV of a occupant.
Args:
	theta_r_is_n: ステップ n における室 i の空気温度, degree C, [i, 1]
	p_a_is_n: ステップ n における室 i の水蒸気圧, Pa, [i, 1]
	h_hum_is_n: ステップ n における室 i の在室者周りの総合熱伝達率, W/m2K, [i, 1]
	theta_ot_is_n: ステップ n における室 i の在室者の作用温度, degree C, [i, 1]
	i_cl_is_n: ステップ n における室 i の在室者の着衣抵抗, m2K/W, [i, 1]
	f_cl_is_n: ステップ n における室 i の在室者の着衣面積率, [i, 1]
	m_is: 室 i の在室者の代謝量, W/m2, [i, 1]
Returns:
	ステップ n における室 i の在室者のPMV, [i, 1]
Notes:
	eq.(1)
*/
func _get_pmv_is_n(
	theta_r_is_n float64,
	p_a_is_n float64,
	h_hum_is_n float64,
	theta_ot_is_n float64,
	i_cl_is_n float64,
	f_cl_is_n float64,
	m_is float64,
) float64 {
	return (0.303*math.Exp(-0.036*m_is) + 0.028) * (m_is - // 活動量, W/m2
		3.05e-3*(5733.0-6.99*m_is-p_a_is_n) - // 皮膚からの潜熱損失, W/m2
		math.Max(0.42*(m_is-58.15), 0.0) - // 発汗熱損失, W/m2
		1.7e-5*m_is*(5867.0-p_a_is_n) - // 呼吸に伴う潜熱損失, W/m2
		0.0014*m_is*(34.0-theta_r_is_n) - // 呼吸に伴う顕熱損失, W/m2 ( = 呼吸量, (g/s)/m2 ✕ (34.0 - 室温)
		f_cl_is_n*h_hum_is_n*(35.7-0.028*m_is-theta_ot_is_n)/(1+i_cl_is_n*f_cl_is_n*h_hum_is_n)) // 着衣からの熱損失
}

/*
目標作用温度を計算する。
Calculate the target operative temperature.
Args:
	p_a_is_n: ステップ n における室 i の水蒸気圧, Pa, [i, 1]
	h_hum_is_n: ステップ n における室 i の在室者周りの総合熱伝達率, W/m2K, [i, 1]
	pmv_target_is_n: ステップ n における室 i の在室者の目標PMV, [i, 1]
	i_cl_is_n: ステップ n における室 i の在室者の着衣抵抗, m2K/W, [i, 1]
	m_is: 室 i の在室者の代謝量, W/m2, [i, 1]
	f_cl_is_n: ステップ n における室 i の在室者の着衣面積率, [i, 1]
Returns:
	ステップ n における室 i の在室者の目標作用温度, degree C, [i, 1]
NOte:
	eq.(3)
*/
func _get_theta_ot_target_is_n(p_a_is_n float64, h_hum_is_n float64, pmv_target_is_n float64, i_cl_is_n float64, m_is, f_cl_is_n float64) float64 {
	return (pmv_target_is_n/(0.303*math.Exp(-0.036*m_is)+0.028) - m_is +
		3.05e-3*(5733.0-6.99*m_is-p_a_is_n) +
		math.Max(0.42*(m_is-58.15), 0.0) +
		1.7e-5*m_is*(5867.0-p_a_is_n) +
		0.0014*m_is*34.0 +
		f_cl_is_n*h_hum_is_n*(35.7-0.028*m_is)/(1+i_cl_is_n*f_cl_is_n*h_hum_is_n)) /
		(0.0014*m_is + f_cl_is_n*h_hum_is_n/(1+i_cl_is_n*f_cl_is_n*h_hum_is_n))
}

/*
在室者の作用温度を計算する。
Calculate the operative temperature of the occupant.
Args:
	h_hum_r_is_n: ステップ n における室 i の在室者周りの放射熱伝達率, W/m2K, [i, 1]
	theta_mrt_is_n: ステップ n における室 i の在室者の平均放射温度, degree C, [i, 1]
	h_hum_c_is_n: ステップ n における室 i の在室者周りの対流熱伝達率, W/m2K, [i, 1]
	theta_r_is_n: ステップ n における室 i の空気温度, degree C, [i, 1]
Returns:
	ステップ n における室 i の在室者の作用温度, degree C, [i, 1]
Notes:
	eq.(5)
*/
func _get_theta_ot_is_n(h_hum_r_is_n float64, theta_mrt_is_n float64, h_hum_c_is_n float64, theta_r_is_n float64) float64 {
	return (h_hum_r_is_n*theta_mrt_is_n + h_hum_c_is_n*theta_r_is_n) / (h_hum_r_is_n + h_hum_c_is_n)
}

/*
人体周りの対流熱伝達率を計算する。（収束計算による方法の場合）
Calculate the convective heat transfer coefficient of the occupant. (Convergence Method)
Args:
	theta_r_is_n: ステップ n における室 i の空気温度, degree C, [i, 1]
	theta_cl_is_n: ステップ n における室 i の在室者の着衣温度, degree C, [i, 1]
	v_hum_is_n: ステップ n における室 i の在室者周りの風速, m/s, [i, 1]
Returns:
	ステップ n の室 i における在室者周りの対流熱伝達率, W/m2K, [i, 1]
Note:
	eq.(5a)
*/
func _get_h_hum_c_is_n_convergence(theta_r_is_n float64, theta_cl_is_n float64, v_hum_is_n float64) float64 {
	return math.Max(12.1*math.Sqrt(v_hum_is_n), 2.38*math.Pow(math.Abs(theta_cl_is_n-theta_r_is_n), 0.25))
}

/*
人体周りの対流熱伝達率を計算する。（収束計算によらず定数を用いる場合）
Calculate the convective heat transfer coefficient of the occupant. (Constant Method)
Args:
	theta_r_is_n: ステップ n における室 i の空気温度, degree C, [i, 1]
Returns:
	ステップ n の室 i における在室者周りの対流熱伝達率, W/m2K, [i, 1]
Note:
	eq.(5b)
*/
func _get_h_hum_c_is_n_constant(theta_r_is_n float64) float64 {
	return 4.0
}

/*
在室者周りの放射熱伝達率を計算する。（収束計算による方法の場合）
Calculate the radiative heat transfer coefficient of the occupnat. (Convergence Method)
Args:
	theta_cl_is_n: ステップ n における室 i の在室者の着衣温度, degree C, [i, 1]
	theta_mrt_is_n: ステップ n における室 i の在室者の平均放射温度, degree C, [i, 1]
Returns:
	ステップ n における室 i の在室者周りの放射熱伝達率, W/m2K, [i, 1]
Note:
	eq.(6a), eq.(6b), eq.(6c)
*/
func _get_h_hum_r_is_n_convergence(theta_cl_is_n float64, theta_mrt_is_n float64) float64 {
	// ステップ n における室 i の在室者の着衣温度, K, [i, 1]
	t_cl_is_n := theta_cl_is_n + 273.0

	// ステップ n における室 i の在室者の平均放射温度, K, [i, 1]
	t_mrt_is_n := theta_mrt_is_n + 273.0

	t_cl_2 := t_cl_is_n * t_cl_is_n
	t_cl_3 := t_cl_2 * t_cl_is_n

	t_mrt_2 := t_mrt_is_n * t_mrt_is_n
	t_mrt_3 := t_mrt_2 * t_mrt_is_n

	return 3.96e-8 * (t_cl_3 + t_cl_2*t_mrt_is_n + t_cl_is_n*t_mrt_2 + t_mrt_3)
}

/*
	人体周りの放射熱伝達率を計算する。（周桑計算によらず定数を用いる場合）
    Calculate the radiative heat transfer coefficient of the occupant. (Constant Method)
    Args:
        theta_r_is_n: ステップ n における室 i の空気温度, degree C, [i, 1]
    Returns:
        ステップ n の室 i における在室者周りの放射熱伝達率, W/m2K, [i, 1]
    Note:
        eq.(6d)
*/
func _get_h_hum_r_is_n_constant(theta_r_is_n float64) float64 {
	const K3 = (20.0 + 273.15) * (20.0 + 273.15) * (20.0 + 273.15)
	return 4 * 3.96e-8 * K3
}

/*
	着衣温度を計算する。
	Calculate the temperature on the clothes of the cuuupant.
	Args:
		clo_is_n: ステップ n における室 i の在室者のClo値, [i, 1]
		theta_ot_is_n: ステップ n における室 i の在室者の作用温度, degree C, [i, 1]
		met_is: 室 i の在室者のMet値, [i, 1]
		h_hum_r_is_n: ステップ n における室 i の在室者周りの放射熱伝達率, W/m2K, [i, 1]
		h_hum_c_is_n: ステップ n における室 i の在室者周りの対流熱伝達率, W/m2K, [i, 1]
	Returns:
		ステップnにおける室iの着衣温度, degree C, [i, 1]
*/
func _get_theta_cl_is_n(
	clo_is_n float64,
	theta_ot_is_n float64,
	met_is float64,
	h_hum_r_is_n float64,
	h_hum_c_is_n float64,
) float64 {
	// 室 i の在室者の代謝量（人体内部発熱量）, W/m2, [i, 1]
	m_is := _get_m_is(met_is)

	// ステップnにおける室iの在室者の着衣抵抗, m2K/W, [i, 1]
	i_cl_is_n := _get_i_cl_is_n(clo_is_n)

	// ステップnにおける室iの在室者の着衣面積率, [i]
	f_cl_is_n := _get_f_cl_is_n(i_cl_is_n)

	// ステップnにおける室iの在室者の着衣温度, degree C
	t_cl_i_n := (35.7-0.028*m_is-theta_ot_is_n)/(1+i_cl_is_n*f_cl_is_n*(h_hum_r_is_n+h_hum_c_is_n)) + theta_ot_is_n

	return t_cl_i_n
}

/*
代謝量を得る。

    Args:
        met_is: 居室 i の在室者のMet値, [i, 1]

    Returns:
        居室 i の在室者の代謝量, W/m2, [i, 1]

    Notes:
        代謝量は1.0 met に固定とする。
        1.0 met は、ISOにおける、Resting - Seated, quiet に相当
        1 met = 58.15 W/m2
*/
func _get_m_is(met float64) float64 {
	return met * 58.15
}

/*
着衣面積率を計算する。

    Args:
        i_cl_is_n: ステップ n における室 i の在室者の着衣抵抗, m2K/W, [i, 1]

    Returns:
        ステップ n における室 i の在室者の着衣面積率, [i, 1]
*/
func _get_f_cl_is_n(i_cl float64) float64 {
	if i_cl <= 0.078 {
		return 1.00 + 1.290*i_cl
	} else {
		return 1.05 + 0.645*i_cl
	}
}

/*
Clo値から着衣抵抗を計算する。

    Args:
        clo_is_n: ステップ n における室 i の在室者のClo値, [i, 1]

    Returns:
        ステップ n における室 i の在室者の着衣抵抗, m2K/W, [i, 1]

    Notes:
        1 clo = 0.155 m2K/W
*/
func _get_i_cl_is_n(clo float64) float64 {
	return clo * 0.155
}
