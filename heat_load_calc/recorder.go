package heat_load_calc

import (
	"fmt"
	"strings"
	"time"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

type Recorder struct {
	YEAR int

	_itv                 Interval          // インターバル
	_id_rm_is            []int             // 室のID
	_id_bdry_js          []int             // 境界のID
	_n_step_i            int               // 瞬時値の行数
	_n_step_a            int               // 平均・積算値の行数
	theta_o_ns           []float64         // ステップ n における外気温度, degree C, [n+1], 出力名："out_temp"
	x_o_ns               []float64         // ステップ n における外気絶対湿度, kg/kg(DA), [n+1], 出力名："out_abs_humid"
	theta_r_is_ns        *mat.Dense        // ステップ n における室 i の室温, degree C, [i, n+1], 出力名："rm[i]_t_r"
	rh_r_is_ns           [][]float64       // ステップ n における室 i の相対湿度, %, [i, n+1], 出力名："rm[i]_rh_r"
	x_r_is_ns            *mat.Dense        // ステップ n における室 i の絶対湿度, kg/kgDA, [i, n+1], 出力名："rm[i]_x_r"
	theta_mrt_hum_is_ns  [][]float64       // ステップ n における室 i の平均放射温度, degree C, [i, n+1], 出力名："rm[i]_mrt"
	theta_ot             [][]float64       // ステップ n における室 i の作用温度, degree C, [i, n+1], 出力名："rm[i]_ot"
	q_trs_sol_is_ns      [][]float64       // ステップ n における室 i の窓の透過日射熱取得, W, [i, n+1], 出力名："rm[i]_q_sol_t"
	theta_frt_is_ns      *mat.Dense        // ステップ n の室 i における家具の温度, degree C, [i, n+1], 出力名："rm[i]_t_fun"
	q_sol_frt_is_ns      [][]float64       // ステップ n の室 i における家具吸収日射熱量, W, [i, n+1], 出力名："rm[i]_q_s_sol_fun"
	x_frt_is_ns          *mat.Dense        // ステップ n の室 i における家具の絶対湿度, kg/kgDA, [i, n+1], 出力名："rm[i]_x_fun"
	pmv_is_ns            [][]float64       // ステップ n の室 i におけるPMV実現値, [i, n+1], 出力名："rm[i]_pmv"
	ppd_is_ns            [][]float64       // ステップ n の室 i におけるPPD実現値, [i, n+1], 出力名："rm[i]_ppd"
	theta_s_js_ns        *mat.Dense        // ステップ n の境界 j の室内側表面温度, degree C, [j, n+1], 出力名:"rm[i]_b[j]_t_s
	theta_ei_js_ns       [][]float64       // ステップ n の境界 j の等価温度, degree C, [j, n+1], 出力名:"rm[i]_b[j]_t_e
	theta_rear_js_ns     [][]float64       // ステップ n の境界 j の裏面温度, degree C, [j, n+1], 出力名:"rm[i]_b[j]_t_b
	h_s_r_js_ns          [][]float64       // ステップ n の境界 j の表面放射熱伝達率, W/m2K, [j, n+1], 出力名:"rm[i]_b[j]_hir_s
	q_r_js_ns            *mat.Dense        // ステップ n の境界 j の表面放射熱流, W, [j, n+1], 出力名:"rm[i]_b[j]_qir_s
	h_s_c_js_ns          [][]float64       // ステップ n の境界 j の表面対流熱伝達率, W/m2K, [j, n+1], 出力名:"rm[i]_b[j]_hic_s
	q_c_js_ns            *mat.Dense        // ステップ n の境界 j の表面対流熱流, W, [j, n+1], 出力名:"rm[i]_b[j]_qic_s
	q_i_sol_s_ns_js      [][]float64       // ステップ n の境界 j の表面日射熱流, W, [j, n+1], 出力名:"rm[i]_b[j]_qisol_s
	q_s_js_ns            [][]float64       // ステップ n の境界 j の表面日射熱流, W, [j, n+1], 出力名:"rm[i]_b[j]_qiall_s
	operation_mode_is_ns [][]OperationMode // ステップ n における室 i の運転状態（平均値）, [i, n], 出力名："rm[i]_ac_operate"
	ac_demand_is_ns      [][]float64       // ステップ n における室 i の空調需要（平均値）, [i, n], 出力名："rm[i]_occupancy"
	h_hum_c_is_ns        [][]float64       // ステップ n における室 i の人体周辺対流熱伝達率（平均値）, W/m2K, [i, n], 出力名："rm[i]_hc_hum"
	h_hum_r_is_ns        [][]float64       // ステップ n における室 i の人体放射熱伝達率（平均値）, W/m2K, [i, n], 出力名："rm[i]_hr_hum"
	q_gen_is_ns          [][]float64       // ステップ n の室 i における人体発熱を除く内部発熱, W, [i, n], 出力名："rm[i]_q_s_except_hum"
	x_gen_is_ns          [][]float64       // ステップ n の室 i における人体発湿を除く内部発湿, kg/s, [i, n], 出力名："rm[i]_q_l_except_hum"
	q_hum_is_ns          [][]float64       // ステップ n の室 i における人体発熱, W, [i, n], 出力名："rm[i]_q_hum_s"
	x_hum_is_ns          [][]float64       // ステップ n の室 i における人体発湿, kg/s, [i, n], 出力名："rm[i]_q_hum_l"
	l_cs_is_ns           [][]float64       // ステップ n の室 i における対流空調顕熱負荷, W, [i, n], 出力名："rm[i]_l_s_c"
	l_rs_is_ns           [][]float64       // ステップ n の室 i における放射空調顕熱負荷, W, [i, n], 出力名："rm[i]_l_s_r"
	l_cl_is_ns           [][]float64       // ステップ n の室 i における対流空調潜熱負荷（加湿側を正とする）, W, [i, n], 出力名："rm[i]_l_l_c"
	q_frt_is_ns          mat.Matrix        // ステップ n の室 i における家具取得熱量, W, [i, n], 出力名："rm[i]_q_s_fun"
	q_l_frt_is_ns        mat.Matrix        // ステップ n の室 i における家具取得水蒸気量, kg/s, [i, n], 出力名："rm[i]_q_l_fun"
	v_reak_is_ns         [][]float64       // ステップ n の室 i におけるすきま風量, m3/s, [i, n], 出力名："rm[i]_v_reak"
	v_ntrl_is_ns         [][]float64       // ステップ n の室 i における自然換気量, m3/s, [i, n], 出力名："rm[i]_v_ntrl"
	v_hum_is_ns          [][]float64       // ステップ　n　の室　i　における人体廻りの風速, m/s, [i, n], 出力名："rm[i]_v_hum"
	clo_is_ns            [][]float64       // ステップ n の室 i におけるClo値, [i, n], 出力名："rm[i]_clo"
	_output_list_room_a  []struct {
		column_name string
		field_name  string
	}
	_output_list_room_i []struct {
		column_name string
		field_name  string
	}
	_output_list_boudary_i []struct {
		column_name string
		field_name  string
	}
}

func make2dim(i int, j int) [][]float64 {
	ar := make([][]float64, i)
	mem := make([]float64, i*j)
	for _i := range ar {
		ar[_i], mem = mem[:j], mem[j:]
	}
	return ar
}

func make2dim_OperationMode(i int, j int) [][]OperationMode {
	ar := make([][]OperationMode, i)
	mem := make([]OperationMode, i*j)
	for _i := range ar {
		ar[_i], mem = mem[:j], mem[j:]
	}
	return ar
}

func NewRecorder(n_step_main int, id_rm_is []int, id_bdry_js []int, itv Interval) *Recorder {
	var r Recorder

	r.YEAR = 1989
	r._itv = itv
	n_rm := len(id_rm_is)
	n_boundries := len(id_bdry_js)
	r._id_rm_is = id_rm_is
	r._id_bdry_js = id_bdry_js
	r._n_step_i = n_step_main + 1
	r._n_step_a = n_step_main

	// --------------- 室に関するもの ---------------

	// ステップ n における外気温度, degree C, [n+1], 出力名："out_temp"
	r.theta_o_ns = make([]float64, r._n_step_i)

	// ステップ n における外気絶対湿度, kg/kg(DA), [n+1], 出力名："out_abs_humid"
	r.x_o_ns = make([]float64, r._n_step_i)

	// ステップ　n　における室　i　の室温, degree C, [i, n+1], 出力名："rm[i]_t_r"
	r.theta_r_is_ns = mat.NewDense(n_rm, r._n_step_i, nil)

	// ステップ n における室 i の相対湿度, %, [i, n+1], 出力名："rm[i]_rh_r"
	r.rh_r_is_ns = make2dim(n_rm, r._n_step_i)

	// ステップ n における室 i の絶対湿度, kg/kgDA, [i, n+1], 出力名："rm[i]_x_r"
	r.x_r_is_ns = mat.NewDense(n_rm, r._n_step_i, nil)

	// ステップ n における室 i の平均放射温度, degree C, [i, n+1], 出力名："rm[i]_mrt"
	r.theta_mrt_hum_is_ns = make2dim(n_rm, r._n_step_i)

	// ステップ n における室 i の作用温度, degree C, [i, n+1], 出力名："rm[i]_ot"
	r.theta_ot = make2dim(n_rm, r._n_step_i)

	// ステップ n における室 i の窓の透過日射熱取得, W, [i, n+1], 出力名："rm[i]_q_sol_t"
	r.q_trs_sol_is_ns = make2dim(n_rm, r._n_step_i)

	// ステップ n の室 i における家具の温度, degree C, [i, n+1], 出力名："rm[i]_t_fun"
	r.theta_frt_is_ns = mat.NewDense(n_rm, r._n_step_i, nil)

	// ステップ n の室 i における家具吸収日射熱量, W, [i, n+1], 出力名："rm[i]_q_s_sol_fun"
	r.q_sol_frt_is_ns = make2dim(n_rm, r._n_step_i)

	// ステップ n の室 i における家具の絶対湿度, kg/kgDA, [i, n+1], 出力名："rm[i]_x_fun"
	r.x_frt_is_ns = mat.NewDense(n_rm, r._n_step_i, nil)

	// ステップ n の室 i におけるPMV実現値, [i, n+1], 出力名："rm[i]_pmv"
	r.pmv_is_ns = make2dim(n_rm, r._n_step_i)

	// ステップ n の室 i におけるPPD実現値, [i, n+1], 出力名："rm[i]_ppd"
	r.ppd_is_ns = make2dim(n_rm, r._n_step_i)

	// --------------- 境界に関するもの ---------------

	// ステップ n の境界 j の室内側表面温度, degree C, [j, n+1], 出力名:"rm[i]_b[j]_t_s
	r.theta_s_js_ns = mat.NewDense(n_boundries, r._n_step_i, nil)

	// ステップ n の境界 j の等価温度, degree C, [j, n+1], 出力名:"rm[i]_b[j]_t_e
	r.theta_ei_js_ns = make2dim(n_boundries, r._n_step_i)

	// ステップ n の境界 j の裏面温度, degree C, [j, n+1], 出力名:"rm[i]_b[j]_t_b
	r.theta_rear_js_ns = make2dim(n_boundries, r._n_step_i)

	// ステップ n の境界 j の表面放射熱流, W, [j, n+1], 出力名:"rm[i]_b[j]_qir_s
	r.h_s_r_js_ns = make2dim(n_boundries, r._n_step_i)

	// ステップ n の境界 j の表面対流熱伝達率, W/m2K, [j, n+1], 出力名:"rm[i]_b[j]_hic_s
	r.q_r_js_ns = mat.NewDense(n_boundries, r._n_step_i, nil)

	// ステップ n の境界 j の表面対流熱流, W, [j, n+1], 出力名:"rm[i]_b[j]_qic_s
	r.h_s_c_js_ns = make2dim(n_boundries, r._n_step_i)

	// ステップ n の境界 j の表面日射熱流, W, [j, n+1], 出力名:"rm[i]_b[j]_qisol_s
	r.q_c_js_ns = mat.NewDense(n_boundries, r._n_step_i, nil)

	// ステップ n の境界 j の表面日射熱流, W, [j, n+1], 出力名:"rm[i]_b[j]_qiall_s
	r.q_i_sol_s_ns_js = make2dim(n_boundries, r._n_step_i)

	// ステップ n の境界 j の係数cvl, degree C, [j, n+1], 出力名:"rm[i]_b[j]_cvl
	r.q_s_js_ns = make2dim(n_boundries, r._n_step_i)

	// --------------- 積算値 ---------------

	// ステップ n における室 i の運転状態（平均値）, [i, n], 出力名："rm[i]_ac_operate"
	r.operation_mode_is_ns = make2dim_OperationMode(n_rm, r._n_step_a)

	// ステップ n における室 i の空調需要（平均値）, [i, n], 出力名："rm[i]_occupancy"
	r.ac_demand_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n における室 i の人体周辺対流熱伝達率（平均値）, W/m2K, [i, n], 出力名："rm[i]_hc_hum"
	r.h_hum_c_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n における室 i の人体放射熱伝達率（平均値）, W/m2K, [i, n], 出力名："rm[i]_hr_hum"
	r.h_hum_r_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における人体発熱を除く内部発熱, W, [i, n], 出力名："rm[i]_q_s_except_hum"
	r.q_gen_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における人体発湿を除く内部発湿, kg/s, [i, n], 出力名："rm[i]_q_l_except_hum"
	r.x_gen_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における人体発熱, W, [i, n], 出力名："rm[i]_q_hum_s"
	r.q_hum_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における人体発湿, kg/s, [i, n], 出力名："rm[i]_q_hum_l"
	r.x_hum_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における対流空調顕熱負荷, W, [i, n], 出力名："rm[i]_l_s_c"
	r.l_cs_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における放射空調顕熱負荷, W, [i, n], 出力名："rm[i]_l_s_r"
	r.l_rs_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における対流空調潜熱負荷（加湿側を正とする）, W, [i, n], 出力名："rm[i]_l_l_c"
	r.l_cl_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における家具取得熱量, W, [i, n], 出力名："rm[i]_q_s_fun"
	r.q_frt_is_ns = mat.NewDense(n_rm, r._n_step_a, nil)

	// ステップ n の室 i における家具取得水蒸気量, kg/s, [i, n], 出力名："rm[i]_q_l_fun"
	r.q_l_frt_is_ns = mat.NewDense(n_rm, r._n_step_a, nil)

	// ステップ n の室 i におけるすきま風量, m3/s, [i, n], 出力名："rm[i]_v_reak"
	r.v_reak_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i における自然換気量, m3/s, [i, n], 出力名："rm[i]_v_ntrl"
	r.v_ntrl_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ　n　の室　i　における人体廻りの風速, m/s, [i, n], 出力名："rm[i]_v_hum"
	r.v_hum_is_ns = make2dim(n_rm, r._n_step_a)

	// ステップ n の室 i におけるClo値, [i, n], 出力名："rm[i]_clo"
	r.clo_is_ns = make2dim(n_rm, r._n_step_a)

	r._output_list_room_a = []struct{ column_name, field_name string }{
		{"operation_mode_is_ns", "ac_operate"},
		{"ac_demand_is_ns", "occupancy"},
		{"h_hum_c_is_ns", "hc_hum"},
		{"h_hum_r_is_ns", "hr_hum"},
		{"q_gen_is_ns", "q_s_except_hum"},
		{"x_gen_is_ns", "q_l_except_hum"},
		{"q_hum_is_ns", "q_hum_s"},
		{"x_hum_is_ns", "q_hum_l"},
		{"l_cs_is_ns", "l_s_c"},
		{"l_rs_is_ns", "l_s_r"},
		{"l_cl_is_ns", "l_l_c"},
		{"q_frt_is_ns", "q_s_fun"},
		{"q_l_frt_is_ns", "q_l_fun"},
		{"v_reak_is_ns", "v_reak"},
		{"v_ntrl_is_ns", "v_ntrl"},
		{"v_hum_is_ns", "v_hum"},
		{"clo_is_ns", "clo"},
	}
	r._output_list_room_i = []struct{ column_name, field_name string }{
		{"theta_r_is_ns", "t_r"},
		{"rh_r_is_ns", "rh_r"},
		{"x_r_is_ns", "x_r"},
		{"theta_mrt_hum_is_ns", "mrt"},
		{"theta_ot", "ot"},
		{"q_trs_sol_is_ns", "q_sol_t"},
		{"theta_frt_is_ns", "t_fun"},
		{"q_sol_frt_is_ns", "q_s_sol_fun"},
		{"x_frt_is_ns", "x_fun"},
		{"pmv_is_ns", "pmv"},
		{"ppd_is_ns", "ppd"},
	}

	r._output_list_boudary_i = []struct{ column_name, field_name string }{
		{"theta_s_js_ns", "t_s"},
		{"theta_ei_js_ns", "t_e"},
		{"theta_rear_js_ns", "t_b"},
		{"h_s_r_js_ns", "hir_s"},
		{"q_r_js_ns", "qir_s"},
		{"h_s_c_js_ns", "hic_s"},
		{"q_c_js_ns", "qic_s"},
		{"q_i_sol_s_ns_js", "qisol_s"},
		{"q_s_js_ns", "qiall_s"},
		{"f_cvl", "f_cvl"},
	}

	return &r
}

func (r *Recorder) pre_recording(
	weather *Weather,
	scd *Schedule,
	bs *Boundaries,
	q_sol_frt_is_ns *mat.Dense,
	q_s_sol_js_ns *mat.Dense,
) {
	// 注意：用意された1年分のデータと実行期間が異なる場合があるためデータスライスする必要がある。

	// ---瞬時値---

	// ステップ n における外気温度, ℃, [n+1]
	copy(r.theta_o_ns, weather.theta_o_ns_plus[0:r._n_step_i])

	// ステップ n における外気絶対湿度, kg/kg(DA), [n+1]
	copy(r.x_o_ns, weather.x_o_ns_plus.RawVector().Data[0:r._n_step_i])

	// ステップ n における室 i の窓の透過日射熱取得, W, [i, n+1]
	for i, _ := range r._id_rm_is {
		// ステップ n における室 i の窓の透過日射熱取得, W, [i, n+1]
		copy(r.q_trs_sol_is_ns[i], bs.q_trs_sol_is_ns.RawRowView(i)[0:r._n_step_i])

		// ステップ n における室 i に設置された備品等による透過日射吸収熱量, W, [i, n+1]
		copy(r.q_sol_frt_is_ns[i], q_sol_frt_is_ns.RawRowView(i)[0:r._n_step_i])
	}

	for i, _ := range r._id_bdry_js {
		// ステップ n の境界 j の表面日射熱流, W, [j, n+1]
		copy(r.q_i_sol_s_ns_js[i], q_s_sol_js_ns.RawRowView(i)[0:r._n_step_i])
		floats.Scale(bs.a_s_js.AtVec(i), r.q_i_sol_s_ns_js[i]) // 面積を掛けておく
	}

	for i := range r._id_bdry_js {
		h_s_c, h_s_r := bs.h_s_c_js.AtVec(i), bs.h_s_r_js.AtVec(i)
		for n := 0; n < r._n_step_i; n++ {
			// ステップ n の境界 j の表面対流熱伝達率, W/m2K, [j, n+1]
			r.h_s_c_js_ns[i][n] = h_s_c

			// ステップ n の境界 j の表面放射熱伝達率, W/m2K, [j, n+1]
			r.h_s_r_js_ns[i][n] = h_s_r
		}
	}

	// ---平均値・積算値---

	off := 0
	for n := 0; n < len(r.ac_demand_is_ns[0]); n++ {
		for i := range r._id_rm_is {
			// ステップ n の室 i における当該時刻の空調需要, [i, n]
			r.ac_demand_is_ns[i][n] = scd.ac_demand_is_ns.Data[off]

			// ステップnの室iにおける人体発熱を除く内部発熱, W, [i, 8760*4]
			r.q_gen_is_ns[i][n] = scd.q_gen_is_ns.Data[off]

			// ステップ n の室 i における人体発湿を除く内部発湿, kg/s, [i, n]
			r.x_gen_is_ns[i][n] = scd.x_gen_is_ns.Data[off]
		}
		off++
	}
}

func (r *Recorder) post_recording(rms *Rooms, bs *Boundaries, f_mrt_is_js *mat.Dense, es *Equipments) {
	// ---瞬時値---

	// 相対湿度を計算する
	p_v_is_ns := make([][]float64, len(r._id_rm_is))
	for i := range r._id_rm_is {
		// ステップ n における室 i の水蒸気圧, Pa, [i, n+1]
		p_v_is_ns[i] = get_p_v_r_is_n(r.x_r_is_ns.RawRowView(i))

		for n := 0; n < r._n_step_i; n++ {
			// ステップ n の室 i における飽和水蒸気圧, Pa, [i, n+1]
			p_vs_is_ns := get_p_vs(r.theta_r_is_ns.At(i, n))

			// ステップnの室iにおける相対湿度, %, [i, n+1]
			r.rh_r_is_ns[i][n] = get_h(p_v_is_ns[i][n], p_vs_is_ns)
		}
	}

	// ステップnの境界jにおける表面熱流（壁体吸熱を正とする）のうち放射成分, W, [j, n]
	var temp1, temp2 mat.Dense
	temp1.Mul(bs.p_js_is, f_mrt_is_js)
	temp2.Mul(&temp1, r.theta_s_js_ns)
	temp2.Sub(&temp2, r.theta_s_js_ns)
	__ScaleRows(&temp2, bs.h_s_r_js)
	__ScaleRows(&temp2, bs.a_s_js)
	r.q_r_js_ns = &temp2
	//ss.h_s_r_js * ss.a_s_js * (np.dot(np.dot(ss.p_js_is, ss.f_mrt_is_js), r.theta_s_js_ns) - r.theta_s_js_ns)

	// ステップnの境界jにおける表面熱流（壁体吸熱を正とする）のうち対流成分, W, [j, n+1]
	var temp3 mat.Dense
	temp3.Mul(bs.p_js_is, r.theta_r_is_ns)
	temp3.Sub(&temp3, r.theta_s_js_ns)
	__ScaleRows(&temp3, bs.h_s_c_js)
	__ScaleRows(&temp3, bs.a_s_js)
	r.q_c_js_ns = &temp3
	// ss.h_s_c_js * ss.a_s_js * (np.dot(ss.p_js_is, r.theta_r_is_ns) - r.theta_s_js_ns)

	// ---平均値・瞬時値---

	// ステップnの室iにおける家具取得熱量, W, [i, n]
	// ステップ n+1 の温度を用いてステップ n からステップ n+1 の平均的な熱流を求めている（後退差分）
	//r.q_frt_is_ns = np.delete(rms.g_sh_frt_is*(r.theta_r_is_ns-r.theta_frt_is_ns), 0, 1)
	var temp5 mat.Dense
	temp5.Sub(r.theta_r_is_ns, r.theta_frt_is_ns)
	__ScaleRows(&temp5, rms.g_sh_frt_is)
	r.q_frt_is_ns = temp5.Slice(0, len(r._id_rm_is), 1, r._n_step_i)

	// ステップ n の室 i の家具等から空気への水分流, kg/s, [i, n]
	// ステップ n+1 の湿度を用いてステップ n からステップ n+1 の平均的な水分流を求めている（後退差分）
	//r.q_l_frt_is_ns = np.delete(rms.g_lh_frt_is*(r.x_r_is_ns-r.x_frt_is_ns), 0, 1)
	var temp6 mat.Dense
	temp6.Sub(r.x_r_is_ns, r.x_frt_is_ns)
	__ScaleRows(&temp6, rms.g_lh_frt_is)
	r.q_l_frt_is_ns = temp6.Slice(0, len(r._id_rm_is), 1, r._n_step_i)

	for i := range r._id_rm_is {
		r.clo_is_ns[i] = get_clo_is_ns(r.operation_mode_is_ns[i])
	}

	// ステップ n+1 のPMVを計算するのに、ステップ n からステップ n+1 のClo値を用いる。
	// 現在、Clo値の配列数が1つ多いバグがあるため、適切な長さになるようにスライスしている。
	// TODO: 本来であれば、助走期間における、n=-1 の時の値を用いないといけないが、とりあえず、配列最後の値を先頭に持ってきて代用している。
	clo_pls := make2dim(len(r._id_rm_is), r._n_step_i)
	for i := range r._id_rm_is {
		copy(clo_pls[i][1:], r.clo_is_ns[i][:r._n_step_i-1])
		clo_pls[i][0] = r.clo_is_ns[i][r._n_step_i-2]
	}

	for i := range r._id_rm_is {
		r.v_hum_is_ns[i] = _get_v_hum_is_n(
			r.operation_mode_is_ns[i],
			es.is_radiative_cooling_is[i],
			es.is_radiative_heating_is[i],
		)
	}

	// ステップ n+1 のPMVを計算するのに、ステップ n からステップ n+1 の人体周りの風速を用いる。
	// TODO: 本来であれば、助走期間における、n=-1 の時の値を用いないといけないが、とりあえず、配列最後の値を先頭に持ってきて代用している。
	v_hum_pls := make2dim(len(r._id_rm_is), r._n_step_i)
	for i := range r._id_rm_is {
		copy(v_hum_pls[i][1:], r.v_hum_is_ns[i][:r._n_step_i-1])
		v_hum_pls[i][0] = r.v_hum_is_ns[i][r._n_step_i-2]
	}

	// ---瞬時値---

	// ステップ n の室 i におけるPMV実現値, [i, n+1]
	var met_is = make([]float64, len(p_v_is_ns[0]))
	for i := range r._id_rm_is {
		for n := 0; n < len(met_is); n++ {
			met_is[n] = rms.met_is[i]
		}
		r.pmv_is_ns[i] = get_pmv_is_n(
			p_v_is_ns[i],
			r.theta_r_is_ns.RawRowView(i),
			r.theta_mrt_hum_is_ns[i],
			clo_pls[i],
			v_hum_pls[i],
			met_is,
			"convergence",
		)
	}

	// ステップ n の室 i におけるPPD実現値, [i, n+1]
	for i := range r._id_rm_is {
		r.ppd_is_ns[i] = get_ppd_is_n(r.pmv_is_ns[i])
	}
}

func (self *Recorder) recording(
	n int,
	theta_r_is_n_pls *mat.VecDense,
	theta_mrt_hum_is_n_pls *mat.VecDense,
	x_r_is_n_pls *mat.VecDense,
	theta_frt_is_n_pls *mat.VecDense,
	x_frt_is_n_pls *mat.VecDense,
	theta_ei_js_n_pls *mat.VecDense,
	q_s_js_n_pls *mat.VecDense,
	theta_ot_is_n_pls *mat.VecDense,
	theta_s_js_n_pls *mat.VecDense,
	theta_rear_js_n_pls *mat.VecDense,
	f_cvl_js_n_pls *mat.VecDense,
	operation_mode_is_n []OperationMode,
	l_cs_is_n *mat.VecDense,
	l_rs_is_n *mat.VecDense,
	l_cl_is_n *mat.VecDense,
	h_hum_c_is_n []float64,
	h_hum_r_is_n []float64,
	q_hum_is_n *mat.VecDense,
	x_hum_is_n *mat.VecDense,
	v_leak_is_n *mat.VecDense,
	v_vent_ntr_is_n []float64,
) {
	// 瞬時値の書き込み

	if n >= -1 {

		// 瞬時値出力のステップ番号
		n_i := n + 1

		r, _ := theta_r_is_n_pls.Dims()
		for i := 0; i < r; i++ {
			// 次の時刻に引き渡す値
			self.theta_r_is_ns.SetCol(n_i, theta_r_is_n_pls.RawVector().Data)
			for i := 0; i < r; i++ {
				self.theta_mrt_hum_is_ns[i][n_i] = theta_mrt_hum_is_n_pls.AtVec(i)
			}
			self.x_r_is_ns.SetCol(n_i, x_r_is_n_pls.RawVector().Data)
			self.theta_frt_is_ns.SetCol(n_i, theta_frt_is_n_pls.RawVector().Data)
			self.x_frt_is_ns.SetCol(n_i, x_frt_is_n_pls.RawVector().Data)
			for i := 0; i < r; i++ {
				self.theta_ei_js_ns[i][n_i] = theta_ei_js_n_pls.AtVec(i)
				self.q_s_js_ns[i][n_i] = q_s_js_n_pls.AtVec(i)
			}

			// 次の時刻に引き渡さない値
			for i := 0; i < r; i++ {
				self.theta_ot[i][n_i] = theta_ot_is_n_pls.AtVec(i)
			}
			self.theta_s_js_ns.SetCol(n_i, theta_s_js_n_pls.RawVector().Data)
			for i := 0; i < r; i++ {
				self.theta_rear_js_ns[i][n_i] = theta_rear_js_n_pls.AtVec(i)
			}
		}
	}

	// 平均値・積算値の書き込み

	if n >= 0 {

		// 平均値出力のステップ番号
		n_a := n

		r, _ := theta_r_is_n_pls.Dims()
		for i := 0; i < r; i++ {
			// 次の時刻に引き渡す値
			self.operation_mode_is_ns[i][n_a] = operation_mode_is_n[i]

			// 次の時刻に引き渡さない値
			// 積算値
			self.l_cs_is_ns[i][n_a] = l_cs_is_n.AtVec(i)
			self.l_rs_is_ns[i][n_a] = l_rs_is_n.AtVec(i)
			self.l_cl_is_ns[i][n_a] = l_cl_is_n.AtVec(i)
			// 平均値
			self.h_hum_c_is_ns[i][n_a] = h_hum_c_is_n[i]
			self.h_hum_r_is_ns[i][n_a] = h_hum_r_is_n[i]
			self.q_hum_is_ns[i][n_a] = q_hum_is_n.AtVec(i)
			self.x_hum_is_ns[i][n_a] = x_hum_is_n.AtVec(i)
			self.v_reak_is_ns[i][n_a] = v_leak_is_n.AtVec(i)
			self.v_ntrl_is_ns[i][n_a] = v_vent_ntr_is_n[i]
		}
	}
}

func (r *Recorder) export_pd() (string, string) {
	var sb_a, sb_i strings.Builder

	// Header
	sb_a.WriteString(strings.Join(r.get_header_a(), ","))
	sb_i.WriteString(strings.Join(r.get_header_i(), ","))
	sb_a.WriteString("\n")
	sb_i.WriteString("\n")

	// 開始日時と終了日時、インターバル
	start, end, freq := r._get_date_index()

	// 現在の日時を開始日時にセット
	current := start

	// 終了日時に到達するまでループ
	n := 0
	for current.Before(end) {
		// 30分を加算
		next := current.Add(time.Duration(freq * float64(time.Minute)))

		// インデックス
		sb_a.WriteString(fmt.Sprintf("%s,%s",
			current.Format("2006-01-02 15:04:05"),
			next.Format("2006-01-02 15:04:05"),
		))
		sb_i.WriteString(current.Format("2006-01-02 15:04:05"))

		// 外気温等
		sb_i.WriteString(fmt.Sprintf(",%g,%g", r.theta_o_ns[n], r.x_o_ns[n]))

		// 室
		for i := range r._id_rm_is {
			// 概要
			sb_a.WriteString(fmt.Sprintf(",OperationMode.%s,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g",
				r.operation_mode_is_ns[i][n],
				r.ac_demand_is_ns[i][n],
				r.h_hum_c_is_ns[i][n],
				r.h_hum_r_is_ns[i][n],
				r.q_gen_is_ns[i][n],
				r.x_gen_is_ns[i][n],
				r.q_hum_is_ns[i][n],
				r.x_hum_is_ns[i][n],
				r.l_cs_is_ns[i][n],
				r.l_rs_is_ns[i][n],
				r.l_cl_is_ns[i][n],
				r.q_frt_is_ns.At(i, n),
				r.q_l_frt_is_ns.At(i, n),
				r.v_reak_is_ns[i][n],
				r.v_ntrl_is_ns[i][n],
				r.v_hum_is_ns[i][n],
				r.clo_is_ns[i][n],
			))

			// 詳細
			sb_i.WriteString(fmt.Sprintf(",%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g",
				r.theta_r_is_ns.At(i, n),
				r.rh_r_is_ns[i][n],
				r.x_r_is_ns.At(i, n),
				r.theta_mrt_hum_is_ns[i][n],
				r.theta_ot[i][n],
				r.q_trs_sol_is_ns[i][n],
				r.theta_frt_is_ns.At(i, n),
				r.q_sol_frt_is_ns[i][n],
				r.x_frt_is_ns.At(i, n),
				r.pmv_is_ns[i][n],
				r.ppd_is_ns[i][n],
			))
		}

		// 境界
		for j := range r._id_bdry_js {
			// 詳細
			sb_i.WriteString(fmt.Sprintf(",%g,%g,%g,%g,%g,%g,%g,%g,%g,%g",
				r.theta_s_js_ns.At(j, n),
				r.theta_ei_js_ns[j][n],
				r.theta_rear_js_ns[j][n],
				r.h_s_r_js_ns[j][n],
				r.q_r_js_ns.At(j, n),
				r.h_s_c_js_ns[j][n],
				r.q_c_js_ns.At(j, n),
				r.q_i_sol_s_ns_js[j][n],
				r.q_s_js_ns[j][n],
				0.0,
			))
		}

		sb_a.WriteString("\n")
		sb_i.WriteString("\n")

		current = next
		n++
	}

	return sb_a.String(), sb_i.String()
}

func (r *Recorder) get_header_i() []string {
	var headers []string

	// インデックス
	headers = append(headers, "start_time")

	// 全体
	headers = append(headers, "out_temp", "out_abs_humid")

	// 室
	for id := range r._id_rm_is {
		for _, item := range r._output_list_room_i {
			headers = append(headers, r._get_room_header_name(id, item.field_name))
		}
	}

	// 境界
	for id := range r._id_bdry_js {
		for _, item := range r._output_list_boudary_i {
			headers = append(headers, r._get_boundary_header_name(id, item.field_name))
		}
	}

	return headers
}

func (r *Recorder) get_header_a() []string {
	var headers []string

	// インデックス
	headers = append(headers, "start_time", "end_time")

	// 室
	for id := range r._id_rm_is {
		for _, item := range r._output_list_room_a {
			headers = append(headers, r._get_room_header_name(id, item.field_name))
		}
	}
	return headers
}

func (r *Recorder) _get_date_index() (time.Time, time.Time, float64) {
	// インターバル(分)
	freq := r._itv.get_pandas_freq()

	// date time index 作成（瞬時値・平均値）
	start := time.Date(r.YEAR, 1, 1, 0, 0, 0, 0, time.Local)
	end := start.Add(time.Duration(float64(r._n_step_i-1) * freq * float64(time.Minute)))

	return start, end, freq
}

func (r *Recorder) _get_room_header_name(id int, name string) string {
	/*room 用のヘッダー名称を取得する。

	Args:
		id: room のID
		name: 出力項目名称

	Returns:
		ヘッダー名称
	*/

	return fmt.Sprintf("rm%d_%s", id, name)
}

func (r *Recorder) _get_room_header_names(name string) []string {
	/*room 用のヘッダー名称を室の数分取得する。

	Args:
		name: 出力項目名称

	Returns:
		ヘッダー名称のリスト
	*/

	var room_headers []string
	for id := range r._id_rm_is {
		room_headers = append(room_headers, r._get_room_header_name(id, name))
	}

	return room_headers
}

func (r *Recorder) _get_boundary_header_name(id int, name string) string {
	/*boundary 用のヘッダ名称を取得する。

	Args:
		pps: PreCalcParameters クラス
		j: boundary の ID
		name: 出力項目名称

	Returns:

	*/

	return fmt.Sprintf("b%d_%s", id, name)
}

func (r *Recorder) _get_boundary_header_names(name string) []string {
	/*boundary 用のヘッダ名称を boundary の数だけ取得する。

	Args:
		name: 出力項目名称

	Returns:

	*/
	var boundary_headers []string
	for id := range r._id_bdry_js {
		boundary_headers = append(boundary_headers, r._get_boundary_header_name(id, name))
	}

	return boundary_headers
}
