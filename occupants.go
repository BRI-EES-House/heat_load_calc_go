package main

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
	1人あたりの人体発湿を計算する。

    Args:
        theta_r_is_n: ステップnの室iにおける室温, degree C, [i, 1]

    Returns:
        ステップnの室iにおける1人あたりの人体発湿, kg/s, [i, 1]
*/
func get_x_hum_psn_is_n(theta_r_is_n mat.Vector) mat.Vector {
	q_hum_psn_is_n := get_q_hum_psn_is_n(theta_r_is_n)

	x_hum_psn_is_n := mat.NewVecDense(theta_r_is_n.Len(), nil)
	for i := 0; i < theta_r_is_n.Len(); i++ {
		x_hum_psn_is_n.SetVec(i, (119.0-q_hum_psn_is_n.AtVec(i))/l_wtr)
	}

	return x_hum_psn_is_n
}

/*
	1人あたりの人体発湿を計算する。

    Args:
        theta_r_is_n: ステップnの室iにおける室温, degree C, [i, 1]

    Returns:
        ステップnの室iにおける1人あたりの人体発熱, W, [i, 1]
*/
func get_q_hum_psn_is_n(theta_r_is_n mat.Vector) mat.Vector {
	q_hum_psn_is_n := mat.NewVecDense(theta_r_is_n.Len(), nil)
	for i := 0; i < theta_r_is_n.Len(); i++ {
		q_hum_psn_is_n.SetVec(i, math.Max(63.0-4.0*(theta_r_is_n.AtVec(i)-24.0), 119.0))
	}
	return q_hum_psn_is_n
}

func get_clo_heavy() float64 {
	return 1.1
}

func get_clo_middle() float64 {
	return 0.7
}

func get_clo_light() float64 {
	return 0.3
}
