package heat_load_calc

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
	1人あたりの人体発湿を計算する。

	Args:
		q_hum_psn_is_n:  ステップnからステップn+1における室iの1人あたりの人体発熱, W, [i]

	Returns:
		ステップnからステップn+1における室iの1人あたりの人体発湿, kg/s, [i, 1]
*/
func get_x_hum_psn_is_n(q_hum_psn_is_n *mat.VecDense) *mat.VecDense {
	x_hum_psn_is_n := mat.NewVecDense(q_hum_psn_is_n.Len(), nil)
	for i := 0; i < x_hum_psn_is_n.Len(); i++ {
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
func get_q_hum_psn_is_n(theta_r_is_n *mat.VecDense) *mat.VecDense {
	q_hum_psn_is_n := mat.NewVecDense(theta_r_is_n.Len(), nil)
	for i := 0; i < q_hum_psn_is_n.Len(); i++ {
		q_hum_psn_is_n.SetVec(i, math.Min(63.0-4.0*(theta_r_is_n.AtVec(i)-24.0), 119.0))
	}
	return q_hum_psn_is_n
}

const clo_heavy = 1.1

const clo_middle = 0.7

const clo_light = 0.3
