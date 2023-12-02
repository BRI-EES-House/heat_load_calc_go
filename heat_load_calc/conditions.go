package heat_load_calc

import "gonum.org/v1/gonum/mat"

type Conditions struct {
	operation_mode_is_n     []OperationMode // ステップnにおける室iの運転状態, [i, 1]
	theta_r_is_n            *mat.VecDense   //ステップnにおける室iの空気温度, degree C, [i, 1]
	theta_mrt_hum_is_n      *mat.VecDense   // ステップnにおける室iの在室者の平均放射温度, degree C, [i, 1]
	x_r_is_n                *mat.VecDense   // ステップnにおける室iの絶対湿度, kg/kgDA, [i, 1]
	theta_dsh_srf_a_js_ms_n *mat.Dense      // ステップnの境界jにおける項別公比法の指数項mの吸熱応答の項別成分, degree C, [j, m] (m=12)
	theta_dsh_srf_t_js_ms_n *mat.Dense      // ステップnの境界jにおける項別公比法の指数項mの貫流応答の項別成分, degree C, [j, m] (m=12)
	q_s_js_n                *mat.VecDense   // ステップnの境界jにおける表面熱流（壁体吸熱を正とする）, W/m2, [j, 1]
	theta_frt_is_n          *mat.VecDense   // ステップnの室iにおける家具の温度, degree C, [i, 1]
	x_frt_is_n              *mat.VecDense   // ステップnの室iにおける家具の絶対湿度, kg/kgDA, [i, 1]
	theta_ei_js_n           *mat.VecDense   // [i, 1]
}

func NewConditions(
	operation_mode_is_n []OperationMode,
	theta_r_is_n *mat.VecDense,
	theta_mrt_hum_is_n *mat.VecDense,
	x_r_is_n *mat.VecDense,
	theta_dsh_srf_a_js_ms_n *mat.Dense,
	theta_dsh_srf_t_js_ms_n *mat.Dense,
	q_s_js_n *mat.VecDense,
	theta_frt_is_n *mat.VecDense,
	x_frt_is_n *mat.VecDense,
	theta_ei_js_n *mat.VecDense,
) *Conditions {
	return &Conditions{
		operation_mode_is_n:     operation_mode_is_n,
		theta_r_is_n:            theta_r_is_n,
		theta_mrt_hum_is_n:      theta_mrt_hum_is_n,
		x_r_is_n:                x_r_is_n,
		theta_dsh_srf_a_js_ms_n: theta_dsh_srf_a_js_ms_n,
		theta_dsh_srf_t_js_ms_n: theta_dsh_srf_t_js_ms_n,
		q_s_js_n:                q_s_js_n,
		theta_frt_is_n:          theta_frt_is_n,
		x_frt_is_n:              x_frt_is_n,
		theta_ei_js_n:           theta_ei_js_n,
	}
}

type GroundConditions struct {
	theta_dsh_srf_a_js_ms_n *mat.Dense // ステップnの境界jにおける項別公比法の指数項mの吸熱応答の項別成分, degree C, [j, m] (m=12)
	theta_dsh_srf_t_js_ms_n *mat.Dense // ステップnの境界jにおける項別公比法の指数項mの貫流応答の項別成分, degree C, [j, m] (m=12)
	q_srf_js_n              []float64  // ステップnの境界jにおける表面熱流（壁体吸熱を正とする）, W/m2, [j, 1]
}

func initialize_conditions(n_spaces, n_bdries int) *Conditions {
	// 空間iの数
	total_number_of_spaces := n_spaces

	// 統合された境界j*の数
	total_number_of_bdry := n_bdries

	// ステップnにおける室iの運転状態, [i, 1]
	// 初期値を暖房・冷房停止で窓「閉」とする。
	operation_mode_is_n := make([]OperationMode, total_number_of_spaces)
	for i := range operation_mode_is_n {
		operation_mode_is_n[i] = STOP_CLOSE
	}

	// ステップnにおける室iの空気温度, degree C, [i, 1]
	// 初期値を15℃とする。
	theta_r_is_n := mat.NewVecDense(total_number_of_spaces, nil)
	for i := 0; i < total_number_of_spaces; i++ {
		theta_r_is_n.SetVec(i, 15.0)
	}

	// ステップnにおける室iの在室者の平均放射温度, degree C, [i, 1]
	// 初期値を15℃と設定する。
	theta_mrt_hum_is_n := mat.NewVecDense(total_number_of_spaces, nil)
	for i := 0; i < total_number_of_spaces; i++ {
		theta_mrt_hum_is_n.SetVec(i, 15.0)
	}

	// ステップnにおける室iの絶対湿度, kg/kgDA, [i, 1]
	// 初期値を空気温度20℃相対湿度40%の時の値とする。
	x_r_is_n := mat.NewVecDense(total_number_of_spaces, nil)
	initial_xr_in := get_x(get_p_vs(20.0) * 0.4)
	for i := 0; i < total_number_of_spaces; i++ {
		x_r_is_n.SetVec(i, initial_xr_in)
	}

	// ステップnの統合された境界j*における指数項mの吸熱応答の項別成分, degree C, [j*, 12]
	// 初期値を0.0℃とする。
	theta_dsh_srf_a_js_ms_n0 := mat.NewDense(total_number_of_bdry, 12, nil)

	// ステップnの統合された境界j*における指数項mの貫流応答の項別成分, degree C, [j*, 12]
	// 初期値を0.0℃とする。
	theta_dsh_srf_t_js_ms_n0 := mat.NewDense(total_number_of_bdry, 12, nil)

	// ステップnの境界jにおける表面熱流（壁体吸熱を正とする）, W/m2, [j, 1]
	// 初期値を0.0W/m2とする。
	q_srf_jstrs_n := mat.NewVecDense(total_number_of_bdry, nil)

	// ステップnの室iにおける家具の温度, degree C, [i]
	// 初期値を15℃とする。
	theta_frt_is_n := mat.NewVecDense(total_number_of_spaces, nil)
	for i := 0; i < total_number_of_spaces; i++ {
		theta_frt_is_n.SetVec(i, 15.0)
	}

	// ステップnの室iにおける家具の絶対湿度, kg/kgDA, [i, 1]
	// 初期値を空気温度20℃相対湿度40%の時の値とする。
	x_frt_is_n := mat.NewVecDense(total_number_of_spaces, nil)
	initial_x_frt_in := get_x(get_p_vs(20.0) * 0.4)
	for i := 0; i < total_number_of_spaces; i++ {
		x_frt_is_n.SetVec(i, initial_x_frt_in)
	}

	// theta_ei_js_n
	theta_ei_js_n := mat.NewVecDense(total_number_of_bdry, nil)
	for i := 0; i < total_number_of_bdry; i++ {
		theta_ei_js_n.SetVec(i, 15.0)
	}

	return NewConditions(
		operation_mode_is_n,
		theta_r_is_n,
		theta_mrt_hum_is_n,
		x_r_is_n,
		theta_dsh_srf_a_js_ms_n0,
		theta_dsh_srf_t_js_ms_n0,
		q_srf_jstrs_n,
		theta_frt_is_n,
		x_frt_is_n,
		theta_ei_js_n,
	)
}

func initialize_ground_conditions(n_grounds int) *GroundConditions {
	if n_grounds == 0 {
		return nil
	}

	// ステップnの統合された境界j*における指数項mの吸熱応答の項別成分, degree C, [j*, 12]
	// 初期値を0.0℃とする。
	theta_dsh_srf_a_js_ms_n0 := mat.NewDense(n_grounds, 12, nil)

	// ステップnの統合された境界j*における指数項mの貫流応答の項別成分, degree C, [j*, 12]
	// 初期値を0.0℃とする。
	theta_dsh_srf_t_js_ms_n0 := mat.NewDense(n_grounds, 12, nil)

	// ステップnの境界jにおける表面熱流（壁体吸熱を正とする）, W/m2, [j, 1]
	// 初期値を0.0W/m2とする。
	q_srf_js_n0 := make([]float64, n_grounds)

	return &GroundConditions{
		theta_dsh_srf_a_js_ms_n: theta_dsh_srf_a_js_ms_n0,
		theta_dsh_srf_t_js_ms_n: theta_dsh_srf_t_js_ms_n0,
		q_srf_js_n:              q_srf_js_n0,
	}
}

func update_conditions_by_ground_conditions(is_ground []bool, c *Conditions, gc *GroundConditions) {
	if gc == nil {
		return
	}
	gidx := 0
	_, l := gc.theta_dsh_srf_a_js_ms_n.Dims()
	for i, ground := range is_ground {
		if ground {
			for j := 0; j < l; j++ {
				c.theta_dsh_srf_a_js_ms_n.Set(i, j, gc.theta_dsh_srf_a_js_ms_n.At(gidx, j))
			}
			c.q_s_js_n.SetVec(i, gc.q_srf_js_n[gidx])
			gidx++
		}
	}
}
