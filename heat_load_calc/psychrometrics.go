package heat_load_calc

import (
	"math"
)

/*
相対湿度を計算する。

    Args:
        p_v: 水蒸気圧, Pa
        p_vs: 飽和水蒸気圧, Pa

    Returns:
        相対湿度, %

    Notes:
        省エネ基準第11章「その他」第5節「湿り空気」式(1)
*/
func get_h(p_v, p_vs float64) float64 {
	return p_v / p_vs * 100.0
}

/*
水蒸気圧から絶対湿度を計算する。

    Args:
        p_v: 水蒸気圧, Pa

    Returns:
        絶対湿度, kg/kgDA

    Notes:
        省エネ基準第11章「その他」第5節「湿り空気」式(2)
        ただし、省エネ基準の式は一部間違っていると考えられる。(2020/1/5時点)
*/
func get_x(p_v float64) float64 {
	f := _get_f()

	return 0.622 * p_v / (f - p_v)
}

/*
絶対湿度から水蒸気圧を求める。

    Args:
        x_r_is_n: ステップnにおける室iの絶対湿度, kg/kgDA, [i]
    Returns:
        ステップnにおける室iの水蒸気圧, Pa, [i]

    Notes:
        省エネ基準第11章「その他」第5節「湿り空気」式(4)
        ただし、省エネ基準の式は絶対湿度の単位として(g/kg(DA))を使用しているが、
        ここでは、kg/kg(DA)に統一した。
*/
// 大気圧, Pa
func get_p_v_r_is_n(x_r_is_n []float64) []float64 {
	f := _get_f()

	p_v_r_is_n := make([]float64, len(x_r_is_n))
	for i := 0; i < len(x_r_is_n); i++ {
		x := x_r_is_n[i]
		p_v_r_is_n[i] = f * x / (x + 0.62198)
	}

	return p_v_r_is_n
}

/*
飽和水蒸気圧を計算する。

    Args:
        theta: 空気温度, degree C

    Returns:
        飽和水蒸気圧, Pa

    Notes:
        省エネ基準
*/
func get_p_vs(theta float64) float64 {
	// 絶対温度の計算
	t := theta + 273.15

	const a1 = -6096.9385
	const a2 = 21.2409642
	const a3 = -0.02711193
	const a4 = 0.00001673952
	const a5 = 2.433502
	const b1 = -6024.5282
	const b2 = 29.32707
	const b3 = 0.010613863
	const b4 = -0.000013198825
	const b5 = -0.49382577

	var p_vs_is float64
	if theta >= 0.0 {
		p_vs_is = math.Exp(a1/t + a2 + a3*t + a4*t*t + a5*math.Log(t))
	} else {
		p_vs_is = math.Exp(b1/t + b2 + b3*t + b4*t*t + b5*math.Log(t))
	}

	return p_vs_is
}

/*
大気圧を求める。

    Returns:
        大気圧, Pa
*/
func _get_f() float64 {
	return 101325.0
}
