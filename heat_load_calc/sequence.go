package main

import (
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

type PreCalcParameters struct {
	// ステップnの室iにおける機械換気量（全般換気量+局所換気量）, m3/s, [i, 8760*4]
	v_vent_mec_is_ns *ScheduleData

	// ステップ n における室 i に設置された備品等による透過日射吸収熱量, W, [i, n+1]
	q_sol_frt_is_ns mat.Matrix

	// ステップnの境界jにおける透過日射熱取得量のうち表面に吸収される日射量, W/m2, [j, 8760*4]
	q_s_sol_js_ns mat.Matrix

	// f_AX のLU分解, [i, i]
	f_ax_js_js *mat.LU

	// 室iの在室者に対する境界j*の形態係数
	f_mrt_hum_is_js mat.Matrix

	// 平均放射温度計算時の境界 j* の表面温度が境界 j に与える重み, [j, j]
	f_mrt_is_js mat.Matrix

	// WSR, WSB の計算 式(24)
	f_wsr_js_is mat.Matrix

	// WSC, W, [j, n]
	f_wsc_js_ns *mat.Dense

	// ステップ n における室 i の在室者表面における放射熱伝達率の総合熱伝達率に対する比, -, [i, 1]
	k_r_is_n mat.Vector

	// ステップnにおける室iの在室者表面における対流熱伝達率の総合熱伝達率に対する比, -, [i, 1]
	k_c_is_n mat.Vector

	// ステップn+1における室iの係数 XOT, [i, i]
	f_xot_is_is_n_pls mat.Matrix
}

type Sequence struct {
	_itv     Interval
	_delta_t float64
	weather  *Weather
	scd      *Schedule
	building *Building
	rms      *Rooms
	bs       *Boundaries // 境界一覧
	mvs      *MechanicalVentilations
	es       *Equipments
	op       *Operation

	/*
		次の係数を求める関数

		ステップ n　からステップ n+1 における係数 f_l_cl_wgt, kg/s(kg/kg(DA)), [i, i]
		ステップ n　からステップ n+1 における係数 f_l_cl_cst, kg/s, [i, 1]
	*/
	get_f_l_cl          func([]float64, []float64, []float64) (*mat.VecDense, *mat.Dense)
	pre_calc_parameters *PreCalcParameters
}

/*
Args:

	itv: 時間間隔
	rd:
	weather:
	scd:
	q_trs_sol_is_ns:
	theta_o_eqv_js_ns:
*/
func NewSequence(
	itv Interval,
	rd *InputJson,
	weather *Weather,
	scd *Schedule,
) *Sequence {
	// 時間間隔, s
	delta_t := itv.get_delta_t()

	// Building Class
	building := CreateBuilding(&rd.Building)

	// Rooms Class
	rms, err := NewRooms(rd.Rooms)
	if err != nil {
		panic(err)
	}

	// Boundaries Class
	bs := NewBoundaries(rms.id_rm_is, rd.Boundaries, weather)

	// MechanicalVentilation Class
	mvs := NewMechanicalVentilations(rd.MechanicalVentilations, rms.n_rm)

	// Equipments Class
	// TODO: Equipments Class を作成するのに Boundaries Class 全部をわたしているのはあまりよくない。
	es := NewEquipments(&rd.Equipments, rms.n_rm, bs.n_b, bs)

	// Operation Class
	op := make_operation(
		&rd.Common,
		scd.ac_setting_is_ns,
		scd.ac_demand_is_ns,
		rms.n_rm,
	)

	// 次の係数を求める関数
	//   ステップ n　からステップ n+1 における係数 f_l_cl_wgt, kg/s(kg/kg(DA)), [i, i]
	//   ステップ n　からステップ n+1 における係数 f_l_cl_cst, kg/s, [i, 1]
	get_f_l_cl := es.make_get_f_l_cl_funcs()

	pre_calc_parameters := _pre_calc(
		scd,
		rms,
		bs,
		mvs,
		es,
		op,
	)

	return &Sequence{

		// 時間間隔クラス
		_itv: itv,

		// 時間間隔, s
		_delta_t: delta_t,

		// Weather Class
		weather: weather,

		// Schedule Class
		scd: scd,

		// Building Class
		building: building,

		// Rooms Class
		rms: rms,

		// Boundaries Class
		bs: bs,

		// MechanicalVentilation Class
		mvs: mvs,

		// Equipments Class
		es: es,

		// Operation Class
		op: op,

		// 次の係数を求める関数
		//   ステップ n　からステップ n+1 における係数 f_l_cl_wgt, kg/s(kg/kg(DA)), [i, i]
		//   ステップ n　からステップ n+1 における係数 f_l_cl_cst, kg/s, [i, 1]
		get_f_l_cl: get_f_l_cl,

		pre_calc_parameters: pre_calc_parameters,
	}
}

func (s *Sequence) run_tick(n int, nn int, nn_plus int, c_n *Conditions, recorder *Recorder) *Conditions {
	ss := s.pre_calc_parameters
	delta_t := s._delta_t

	return _run_tick(s, n, nn, nn_plus, delta_t, ss, c_n, recorder)
}

func (s *Sequence) run_tick_ground(gc_n *GroundConditions, n int, nn int) *GroundConditions {

	pp := s.pre_calc_parameters

	return _run_tick_ground(s, pp, gc_n, n, nn)
}

/*
Args:

	f_ax_js_js: 係数 f_{AX}, -, [j, j]
	f_crx_js_ns: 係数 f_{CRX,n}, degree C, [j, n]

Returns:

	係数 f_{WSC,n}, degree C, [j, n]

Notes:

	式(4.1)
*/
func get_f_wsc_js_ns(f_ax_js_js *mat.LU, f_crx_js_ns mat.Matrix) *mat.Dense {
	var temp1 mat.Dense
	f_ax_js_js.SolveTo(&temp1, false, f_crx_js_ns)
	return &temp1
}

/*
Args:

	f_ax_js_js: 係数 f_AX, -, [j, j]
	f_fia_js_is: 係数 f_FIA, -, [j, i]

Returns:

	係数 f_WSR, -, [j, i]

Notes:

	式(4.2)
*/
func get_f_wsr_js_is(f_ax_js_js *mat.LU, f_fia_js_is mat.Matrix) *mat.Dense {
	var temp1 mat.Dense
	f_ax_js_js.SolveTo(&temp1, false, f_fia_js_is)
	return &temp1
}

/*
Args:

	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	h_s_r_js: 境界 j の室内側放射熱伝達率, W/(m2 K), [j, 1]
	k_ei_js_js: 境界 j の裏面温度に境界　j∗ の等価温度が与える影響, -, [j, j]
	phi_a0_js: 境界 j の吸熱応答係数の初項, m2 K/W, [j, 1]
	phi_t0_js: 境界 j の貫流応答係数の初項, -, [j, 1]
	q_s_sol_js_ns: ステップ n における境界 j の透過日射吸収熱量, W/m2, [j, n]
	k_eo_js: 境界 j の裏面温度に境界 j の相当外気温度が与える影響, -, [j, 1]
	theta_o_eqv_js_ns: ステップ n における境界 j の相当外気温度, degree C, [j, 1]

Returns:

	係数 f_CRX, degree C, [j, n]

Notes:

	式(4.3)
*/
func get_f_crx_js_ns(
	h_s_c_js mat.Vector,
	h_s_r_js mat.Vector,
	k_ei_js_js mat.Matrix,
	phi_a0_js mat.Vector,
	phi_t0_js mat.Vector,
	q_s_sol_js_ns mat.Matrix,
	k_eo_js mat.Vector,
	theta_o_eqv_js_ns mat.Matrix,
) *mat.Dense {
	// q_s_sol_js_ns/(h_s_c_js+h_s_r_js)
	var temp2 mat.Dense
	temp2.Apply(func(j int, i int, v float64) float64 {
		return v / (h_s_c_js.AtVec(j) + h_s_r_js.AtVec(j))
	}, q_s_sol_js_ns)

	// np.dot(k_ei_js_js, temp2)
	var temp3 mat.Dense
	temp3.Mul(k_ei_js_js, &temp2)

	// phi_a0_js*q_s_sol_js_ns + phi_t0js*temp3+phi_t0_js*theta_o_eqv_js_ns*k_eo_js
	var result mat.Dense
	result.Apply(func(j, i int, v float64) float64 {
		return phi_a0_js.AtVec(j)*v +
			phi_t0_js.AtVec(j)*temp3.At(j, i) +
			phi_t0_js.AtVec(j)*theta_o_eqv_js_ns.At(j, i)*k_eo_js.AtVec(j)
	}, q_s_sol_js_ns)

	return &result
}

/*
Args:

	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	h_s_r_js: 境界 j の室内側放射熱伝達率, W/(m2 K), [j, 1]
	k_ei_js_js: 境界 j の裏面温度に境界　j∗ の等価温度が与える影響, -, [j, j]
	p_js_is: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [j, i]
	phi_a0_js: 境界 j の吸熱応答係数の初項, m2 K/W, [j, 1]
	phi_t0_js: 境界 j の貫流応答係数の初項, -, [j, 1]

Returns:

	係数 f_FIA, -, [j, i]

Notes:

	式(4.4)
*/
func get_f_fia_js_is(
	h_s_c_js mat.Vector,
	h_s_r_js mat.Vector,
	k_ei_js_js mat.Matrix,
	p_js_is mat.Matrix,
	phi_a0_js mat.Vector,
	phi_t0_js mat.Vector,
	k_s_r_js_is mat.Matrix,
) *mat.Dense {

	// phi_a0_js * h_s_c_js * p_js_is
	var temp1 mat.Dense
	temp1.Apply(func(j, i int, v float64) float64 {
		return phi_a0_js.AtVec(j) * h_s_c_js.AtVec(j) * v
	}, p_js_is)

	// np.dot(k_ei_js_js, p_js_is) * phi_t0_js * h_s_c_js / (h_s_c_js + h_s_r_js)
	var temp2 mat.Dense
	temp2.Mul(k_ei_js_js, p_js_is)
	temp2.Apply(func(j, i int, v float64) float64 {
		return v * phi_t0_js.AtVec(j) * h_s_c_js.AtVec(j) / (h_s_c_js.AtVec(j) + h_s_r_js.AtVec(j))
	}, &temp2)

	// phi_t0_js * k_s_r_js_is
	var temp3 mat.Dense
	temp3.Apply(func(j, i int, v float64) float64 {
		return phi_t0_js.AtVec(j) * v
	}, k_s_r_js_is)

	temp1.Add(&temp1, &temp2)
	temp1.Add(&temp1, &temp3)

	return &temp1
}

/*
Args:

	f_mrt_is_js: 室 i の微小球に対する境界 j の形態係数, -, [i, j]
	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	h_s_r_js: 境界 j の室内側放射熱伝達率, W/(m2 K), [j, 1]
	k_ei_js_js: 境界 j の裏面温度に境界　j∗ の等価温度が与える影響, -, [j, j]
	p_js_is: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [j, i]
	phi_a0_js: 境界 j の吸熱応答係数の初項, m2 K/W, [j, 1]
	phi_t0_js: 境界 j の貫流応答係数の初項, -, [j, 1]

Returns:

	係数 f_AX, -, [j, j]

Notes:

	式(4.5)
*/
func get_f_ax_js_is(
	f_mrt_is_js mat.Matrix,
	h_s_c_js mat.Vector,
	h_s_r_js mat.Vector,
	k_ei_js_js mat.Matrix,
	p_js_is mat.Matrix,
	phi_a0_js mat.Vector,
	phi_t0_js mat.Vector,
) *mat.LU {
	// 1.0 + phi_a0_js * (h_s_r_js + p_js_is)
	temp1 := make([]float64, phi_a0_js.Len())
	for i := 0; i < phi_a0_js.Len(); i++ {
		temp1[i] = 1.0 + phi_a0_js.AtVec(i)*(h_s_c_js.AtVec(i)+h_s_r_js.AtVec(i))
	}

	// diag of temp1
	temp2 := mat.NewDiagDense(phi_a0_js.Len(), temp1)

	// h_s_r_js * phi_a0_js
	var temp3 mat.VecDense
	temp3.MulElemVec(h_s_r_js, phi_a0_js)

	// np.dot(p_js_is, f_mrt_is_js)
	var temp4 mat.Dense
	temp4.Mul(p_js_is, f_mrt_is_js)

	// temp4 * temp3
	var temp5 mat.Dense
	temp5.Apply(func(i, j int, v float64) float64 {
		return v * temp3.AtVec(i)
	}, &temp4)

	// np.dot(k_ei_js_js, np.dot(p_js_is, f_mrt_is_js))
	var temp6 mat.Dense
	temp6.Mul(k_ei_js_js, &temp4)

	// h_s_r_js * phi_t0_js / (h_s_c_js + h_s_r_js)
	var temp7 mat.VecDense
	temp7.AddVec(h_s_c_js, h_s_r_js)
	temp7.MulElemVec(&temp7, phi_t0_js)
	temp7.MulElemVec(&temp7, h_s_r_js)

	// temp6 * temp7
	temp6.Apply(func(i, j int, v float64) float64 {
		return v * temp7.AtVec(i)
	}, &temp6)

	// temp2 - temp5 - temp6
	var result mat.Dense
	result.Sub(temp2, &temp5)
	result.Sub(&result, &temp6)

	// 高速化のために、ここで LU 分解を行う
	var lu mat.LU
	lu.Factorize(&result)

	return &lu
}

/*
Args:

	v_vent_mec_general_is: ステップ n からステップ n+1 における室 i の機械換気量（全般換気量）, m3/s, [i, 1]
	v_vent_mec_local_is_ns: ステップ n からステップ n+1 における室 i の機械換気量（局所換気量）, m3/s, [i, 1]

Returns:

	ステップ n からステップ n+1 における室 i の機械換気量（全般換気量と局所換気量の合計値）, m3/s, [i, 1]

Notes:

	式(4.7)
*/
func get_v_vent_mec_is_ns(
	v_vent_mec_general_is []float64,
	v_vent_mec_local_is_ns *ScheduleData,
) *ScheduleData {
	v_vent_mec_is_ns := make([]float64, len(v_vent_mec_local_is_ns.Data))
	off := 0
	for i := 0; i < v_vent_mec_local_is_ns.Len(); i++ {
		for j := 0; j < v_vent_mec_local_is_ns.BatchSize; j++ {
			v_vent_mec_is_ns[i] = v_vent_mec_local_is_ns.Data[off] + v_vent_mec_general_is[j]
		}
	}

	return &ScheduleData{
		Data:      v_vent_mec_is_ns,
		BatchSize: v_vent_mec_local_is_ns.BatchSize,
	}
}

/*
助走計算用パラメータの生成

Args:

	scd: Scheduleクラス
	rms: Roomsクラス
	bs: Boundariesクラス
	mvs: MechanicalVentilationsクラス
	es: Equipmenstクラス

Returns:

	PreCalcParametersクラス
*/
func _pre_calc(
	scd *Schedule,
	rms *Rooms,
	bs *Boundaries,
	mvs *MechanicalVentilations,
	es *Equipments,
	op *Operation,
) *PreCalcParameters {
	// 室 i の在室者に対する境界jの形態係数, [i, j]
	f_mrt_hum_is_js := get_f_mrt_hum_js(bs.p_is_js, bs.a_s_js, bs.is_floor_js)

	// 室 i の微小球に対する境界 j の形態係数, -, [i, j]
	f_mrt_is_js := get_f_mrt_is_js(bs.a_s_js, bs.h_s_r_js, bs.p_is_js)

	// ステップ n からステップ n+1 における室 i の機械換気量（全般換気量と局所換気量の合計値）, m3/s, [i, 1]
	v_vent_mec_is_ns := get_v_vent_mec_is_ns(mvs.v_vent_mec_general_is, scd.v_mec_vent_local_is_ns)

	// ステップ n における室 i に設置された備品等による透過日射吸収熱量, W, [i, n+1]
	q_sol_frt_is_ns := get_q_sol_frt_is_ns(bs.q_trs_sol_is_ns)

	// ステップ n における境界 j の透過日射吸収熱量, W/m2, [j, n]
	q_s_sol_js_ns := get_q_s_sol_js_ns(
		bs.p_is_js,
		bs.a_s_js,
		bs.p_s_sol_abs_js,
		bs.p_js_is,
		bs.q_trs_sol_is_ns,
	)

	// 係数 f_AX, -, [j, j]
	f_ax_js_js := get_f_ax_js_is(
		f_mrt_is_js,
		bs.h_s_c_js,
		bs.h_s_r_js,
		bs.k_ei_js_js,
		bs.p_js_is,
		bs.phi_a0_js,
		bs.phi_t0_js,
	)

	// 係数 f_FIA, -, [j, i]
	f_fia_js_is := get_f_fia_js_is(
		bs.h_s_c_js,
		bs.h_s_r_js,
		bs.k_ei_js_js,
		bs.p_js_is,
		bs.phi_a0_js,
		bs.phi_t0_js,
		bs.k_s_r_js,
	)

	// 係数 f_CRX, degree C, [j, n]
	f_crx_js_ns := get_f_crx_js_ns(
		bs.h_s_c_js,
		bs.h_s_r_js,
		bs.k_ei_js_js,
		bs.phi_a0_js,
		bs.phi_t0_js,
		q_s_sol_js_ns,
		bs.k_eo_js,
		bs.theta_o_eqv_js_ns,
	)

	// 係数 f_WSR, -, [j, i]
	f_wsr_js_is := get_f_wsr_js_is(f_ax_js_js, f_fia_js_is)

	// 係数 f_{WSC, n}, degree C, [j, n]
	f_wsc_js_ns := get_f_wsc_js_ns(f_ax_js_js, f_crx_js_ns)

	// ステップnにおける室iの在室者表面における対流熱伝達率の総合熱伝達率に対する比, -, [i, 1]
	// ステップ n における室 i の在室者表面における放射熱伝達率の総合熱伝達率に対する比, -, [i, 1]
	k_c_is_n, k_r_is_n := op.get_k_is()

	// ステップn+1における室iの係数 XOT, [i, i]
	f_xot_is_is_n_pls := get_f_xot_is_is_n_pls(
		f_mrt_hum_is_js,
		f_wsr_js_is,
		k_c_is_n,
		k_r_is_n,
	)

	pre_calc_parameters := PreCalcParameters{
		v_vent_mec_is_ns:  v_vent_mec_is_ns,
		f_mrt_hum_is_js:   f_mrt_hum_is_js,
		f_mrt_is_js:       f_mrt_is_js,
		q_s_sol_js_ns:     q_s_sol_js_ns,
		q_sol_frt_is_ns:   q_sol_frt_is_ns,
		f_wsr_js_is:       f_wsr_js_is,
		f_ax_js_js:        f_ax_js_js,
		f_wsc_js_ns:       f_wsc_js_ns,
		k_r_is_n:          k_r_is_n,
		k_c_is_n:          k_c_is_n,
		f_xot_is_is_n_pls: f_xot_is_is_n_pls,
	}

	return &pre_calc_parameters
}

// /*
// 気象データを読み込む。
// Args:
//     pp (pd.DataFrame): 気象データのDataFrame
// Returns:
//     外気温度, degree C
//     外気絶対湿度, kg/kg(DA)
//     法線面直達日射量, W/m2
//     水平面天空日射量, W/m2
//     夜間放射量, W/m2
//     太陽高度, rad
//     太陽方位角, rad
// */
// func _read_weather_data(pp *WeatherPdFrame) ([]float64, []float64, []float64, []float64, []float64, []float64, []float64) {

// 	// 外気温度, degree C
// 	theta_o_ns := pp.temperature

// 	// 外気絶対湿度, kg/kg(DA)
// 	x_o_ns := pp.absolute_humidity

// 	// 法線面直達日射量, W/m2
// 	i_dn_ns := pp.normal_direct_solar_radiation

// 	// 水平面天空日射量, W/m2
// 	i_sky_ns := pp.horizontal_sky_solar_radiation

// 	// 夜間放射量, W/m2
// 	r_n_ns := pp.outward_radiation

// 	// 太陽高度, rad
// 	h_sun_ns := pp.sun_altitude

// 	// 太陽方位角, rad
// 	a_sun_ns := pp.sun_azimuth

// 	return a_sun_ns, h_sun_ns, i_dn_ns, i_sky_ns, r_n_ns, theta_o_ns, x_o_ns
// }

/*
室の温湿度・熱負荷の計算
Args:

	n: ステップ
	delta_t: 時間間隔, s
	ss: ループ計算前に計算可能なパラメータを含めたクラス
	c_n: 前の時刻からの状態量
	recorder: Recorder クラス

Returns:

	次の時刻にわたす状態量
*/
func _run_tick(self *Sequence, n int, nn int, nn_plus int, delta_t float64, ss *PreCalcParameters, c_n *Conditions, recorder *Recorder) *Conditions {

	// ----------- 人体発熱・人体発湿 -----------

	// ステップnからステップn+1における室iの1人あたりの人体発熱, W, [i, 1]
	q_hum_psn_is_n := get_q_hum_psn_is_n(c_n.theta_r_is_n)

	// ステップ n からステップ n+1 における室 i の人体発熱, W, [i, 1]
	q_hum_is_n := get_q_hum_is_n(
		self.scd.n_hum_is_ns.Get(nn),
		q_hum_psn_is_n,
	)

	// ステップnの室iにおける1人あたりの人体発湿, kg/s, [i, 1]
	x_hum_psn_is_n := get_x_hum_psn_is_n(c_n.theta_r_is_n, q_hum_psn_is_n)

	// ステップnの室iにおける人体発湿, kg/s, [i, 1]
	x_hum_is_n := get_x_hum_is_n(self.scd.n_hum_is_ns.Get(nn), x_hum_psn_is_n)

	// --------------------------------------------

	// ステップnの境界jにおける裏面温度, deg    // ステップ n の境界 j における裏面温度, degree C, [j, 1]
	theta_rear_js_n := get_theta_s_rear_js_n(
		self.bs.k_ei_js_js,
		c_n.theta_ei_js_n,
		self.bs.k_eo_js,
		self.bs.theta_o_eqv_js_ns.ColView(nn_plus),
		self.bs.k_s_r_js,
		c_n.theta_r_is_n,
	)

	// ステップnの室iにおけるすきま風量, m3/s, [i, 1]
	v_leak_is_n := self.building.get_v_leak_is_n(
		c_n.theta_r_is_n,
		self.weather.theta_o_ns_plus[nn_plus],
		self.rms.v_rm_is,
	)

	// ステップ n+1 の境界 j における項別公比法の指数項 m の貫流応答の項別成分, degree C, [j, m] (m=12), eq.(29)
	theta_dsh_s_t_js_ms_n_pls := get_theta_dsh_s_t_js_ms_n_pls(
		self.bs.phi_t1_js_ms,
		self.bs.r_js_ms,
		c_n.theta_dsh_srf_t_js_ms_n,
		theta_rear_js_n,
	)

	// ステップ n+1 の境界 j における項別公比法の指数項 m の吸熱応答の項別成分, degree C, [j, m]
	theta_dsh_s_a_js_ms_n_pls := get_theta_dsh_s_a_js_ms_n_pls(
		self.bs.phi_a1_js_ms,
		c_n.q_s_js_n,
		self.bs.r_js_ms,
		c_n.theta_dsh_srf_a_js_ms_n,
	)

	// ステップ n+1 の境界 j における係数f_CVL, degree C, [j, 1]
	f_cvl_js_n_pls := get_f_cvl_js_n_pls(
		theta_dsh_s_a_js_ms_n_pls,
		theta_dsh_s_t_js_ms_n_pls,
	)

	// ステップ n+1 の境界 j における係数 f_WSV, degree C, [j, 1]
	f_wsv_js_n_pls := get_f_wsv_js_n_pls(
		f_cvl_js_n_pls,
		ss.f_ax_js_js,
	)

	// ステップnからステップn+1における室iの換気・隙間風による外気の流入量, m3/s, [i, 1]
	v_vent_out_non_nv_is_n := get_v_vent_out_non_ntr_is_n(
		v_leak_is_n,
		ss.v_vent_mec_is_ns.Get(nn),
	)

	// ステップ n+1 の室 i における係数 f_BRC, W, [i, 1]
	// TODO: q_sol_frt_is_ns の値は n+1 の値を使用するべき？
	f_brc_non_nv_is_n_pls, f_brc_nv_is_n_pls := get_f_brc_is_n_pls(
		self.bs.a_s_js,
		self.rms.v_rm_is,
		self.rms.c_sh_frt_is,
		delta_t,
		ss.f_wsc_js_ns.ColView(nn_plus+1),
		f_wsv_js_n_pls,
		self.rms.g_sh_frt_is,
		self.bs.h_s_c_js,
		self.bs.p_is_js,
		self.scd.q_gen_is_ns.Get(nn),
		q_hum_is_n,
		ss.q_sol_frt_is_ns.(mat.ColViewer).ColView(nn_plus),
		c_n.theta_frt_is_n,
		self.weather.theta_o_ns_plus[nn_plus+1],
		c_n.theta_r_is_n,
		v_vent_out_non_nv_is_n,
		self.rms.v_vent_ntr_set_is,
	)

	// ステップ n+1 における係数 f_BRM, W/K, [i, i]
	f_brm_non_nv_is_is_n_pls, f_brm_nv_is_is_n_pls := get_f_brm_is_is_n_pls(
		self.bs.a_s_js,
		self.rms.v_rm_is,
		self.rms.c_sh_frt_is,
		delta_t,
		ss.f_wsr_js_is,
		self.rms.g_sh_frt_is,
		self.bs.h_s_c_js,
		self.bs.p_is_js,
		self.bs.p_js_is,
		self.mvs.v_vent_int_is_is,
		v_vent_out_non_nv_is_n,
		self.rms.v_vent_ntr_set_is,
	)

	// ステップn+1における室iの係数 XC, [i, 1]
	f_xc_is_n_pls := get_f_xc_is_n_pls(
		ss.f_mrt_hum_is_js,
		ss.f_wsc_js_ns.ColView(nn_plus+1),
		f_wsv_js_n_pls,
		ss.f_xot_is_is_n_pls,
		ss.k_r_is_n,
	)

	// ステップ n における係数 f_BRM,OT, W/K, [i, i]
	f_brm_ot_non_nv_is_is_n_pls, f_brm_ot_nv_is_is_n_pls := get_f_brm_ot_is_is_n_pls(
		ss.f_xot_is_is_n_pls,
		f_brm_non_nv_is_is_n_pls,
		f_brm_nv_is_is_n_pls,
	)

	// ステップ n における係数 f_BRC,OT, W, [i, 1]
	f_brc_ot_non_nv_is_n_pls, f_brc_ot_nv_is_n_pls := get_f_brc_ot_is_n_pls(
		f_xc_is_n_pls,
		f_brc_non_nv_is_n_pls,
		f_brc_nv_is_n_pls,
		f_brm_non_nv_is_is_n_pls,
		f_brm_nv_is_is_n_pls,
	)

	// ステップnにおける室iの自然風の非利用時の潜熱バランスに関する係数f_h_cst, kg / s, [i, 1]
	// ステップnにおける室iの自然風の利用時の潜熱バランスに関する係数f_h_cst, kg / s, [i, 1]
	f_h_cst_non_nv_is_n, f_h_cst_nv_is_n := get_f_h_cst_is_n(
		self.rms.c_lh_frt_is,
		delta_t,
		self.rms.g_lh_frt_is,
		rho_a,
		self.rms.v_rm_is,
		c_n.x_frt_is_n,
		self.scd.x_gen_is_ns.Get(nn),
		x_hum_is_n,
		self.weather.x_o_ns_plus.AtVec(nn_plus+1),
		c_n.x_r_is_n,
		v_vent_out_non_nv_is_n,
		self.rms.v_vent_ntr_set_is,
	)

	// ステップnにおける自然風非利用時の室i*の絶対湿度が室iの潜熱バランスに与える影響を表す係数,　kg/(s kg/kg(DA)), [i, i]
	// ステップnにおける自然風利用時の室i*の絶対湿度が室iの潜熱バランスに与える影響を表す係数,　kg/(s kg/kg(DA)), [i, i]
	f_h_wgt_non_nv_is_is_n, f_h_wgt_nv_is_is_n := get_f_h_wgt_is_is_n(
		self.rms.c_lh_frt_is,
		delta_t,
		self.rms.g_lh_frt_is,
		self.rms.v_rm_is,
		self.mvs.v_vent_int_is_is,
		v_vent_out_non_nv_is_n,
		self.rms.v_vent_ntr_set_is,
	)

	// ステップn+1における自然風非利用時の自然作用温度, degree C, [i, 1]
	// ステップn+1における自然風利用時の自然作用温度, degree C, [i, 1]
	theta_r_ot_ntr_non_nv_is_n_pls, theta_r_ot_ntr_nv_is_n_pls := get_theta_r_ot_ntr_is_n_pls(
		f_brc_ot_non_nv_is_n_pls,
		f_brc_ot_nv_is_n_pls,
		f_brm_ot_non_nv_is_is_n_pls,
		f_brm_ot_nv_is_is_n_pls,
	)

	_f_xc_is_n_pls := mat.NewVecDense(len(f_xc_is_n_pls), f_xc_is_n_pls)

	var theta_r_ntr_non_nv_is_n_pls, theta_r_ntr_nv_is_n_pls mat.VecDense
	theta_r_ntr_non_nv_is_n_pls.MulVec(ss.f_xot_is_is_n_pls, theta_r_ot_ntr_non_nv_is_n_pls)
	theta_r_ntr_non_nv_is_n_pls.SubVec(&theta_r_ntr_non_nv_is_n_pls, _f_xc_is_n_pls)
	theta_r_ntr_nv_is_n_pls.MulVec(ss.f_xot_is_is_n_pls, theta_r_ot_ntr_nv_is_n_pls)
	theta_r_ntr_nv_is_n_pls.SubVec(&theta_r_ntr_nv_is_n_pls, _f_xc_is_n_pls)

	var theta_s_ntr_non_nv_js_n_pls, theta_s_ntr_nv_js_n_pls mat.VecDense
	theta_s_ntr_non_nv_js_n_pls.MulVec(ss.f_wsr_js_is, &theta_r_ntr_non_nv_is_n_pls)
	theta_s_ntr_non_nv_js_n_pls.AddVec(&theta_s_ntr_non_nv_js_n_pls, ss.f_wsc_js_ns.ColView(nn_plus+1))
	theta_s_ntr_non_nv_js_n_pls.AddVec(&theta_s_ntr_non_nv_js_n_pls, f_wsv_js_n_pls)
	theta_s_ntr_nv_js_n_pls.MulVec(ss.f_wsr_js_is, &theta_r_ntr_nv_is_n_pls)
	theta_s_ntr_nv_js_n_pls.AddVec(&theta_s_ntr_nv_js_n_pls, ss.f_wsc_js_ns.ColView(nn_plus+1))
	theta_s_ntr_nv_js_n_pls.AddVec(&theta_s_ntr_nv_js_n_pls, f_wsv_js_n_pls)

	var theta_mrt_hum_ntr_non_nv_is_n_pls, theta_mrt_hum_ntr_nv_is_n_pls mat.VecDense
	theta_mrt_hum_ntr_non_nv_is_n_pls.MulVec(ss.f_mrt_is_js, &theta_s_ntr_non_nv_js_n_pls)
	theta_mrt_hum_ntr_nv_is_n_pls.MulVec(ss.f_mrt_is_js, &theta_s_ntr_nv_js_n_pls)

	// ステップn+1における室iの自然風非利用時の加湿・除湿を行わない場合の絶対湿度, kg/kg(DA) [i, 1]
	// ステップn+1における室iの自然風利用時の加湿・除湿を行わない場合の絶対湿度, kg/kg(DA) [i, 1]
	x_r_ntr_non_nv_is_n_pls, x_r_ntr_nv_is_n_pls := get_x_r_ntr_is_n_pls(
		f_h_cst_non_nv_is_n,
		f_h_wgt_non_nv_is_is_n,
		f_h_cst_nv_is_n,
		f_h_wgt_nv_is_is_n,
	)

	// ステップ n における室 i の運転モード, [i, 1]
	operation_mode_is_n := self.op.get_operation_mode_is_n(
		n,
		nn,
		self.es.is_radiative_heating_is,
		self.es.is_radiative_cooling_is,
		self.rms.met_is,
		theta_r_ot_ntr_non_nv_is_n_pls,
		theta_r_ot_ntr_nv_is_n_pls,
		theta_r_ntr_non_nv_is_n_pls.RawVector().Data,
		theta_r_ntr_nv_is_n_pls.RawVector().Data,
		theta_mrt_hum_ntr_non_nv_is_n_pls.RawVector().Data,
		theta_mrt_hum_ntr_nv_is_n_pls.RawVector().Data,
		x_r_ntr_non_nv_is_n_pls.RawVector().Data,
		x_r_ntr_nv_is_n_pls.RawVector().Data,
	)

	f_brm_is_is_n_pls := mat.NewDense(self.rms.n_rm, self.rms.n_rm, nil)
	v_vent_ntr_is_n := make([]float64, self.rms.n_rm)
	f_brm_ot_is_is_n_pls := mat.NewDense(self.rms.n_rm, self.rms.n_rm, nil)
	f_brc_ot_is_n_pls := make([]float64, self.rms.n_rm)
	f_h_cst_is_n := make([]float64, self.rms.n_rm)
	f_h_wgt_is_is_n := mat.NewDense(self.rms.n_rm, self.rms.n_rm, nil)
	theta_r_ot_ntr_is_n_pls := make([]float64, self.rms.n_rm)
	theta_r_ntr_is_n_pls := make([]float64, self.rms.n_rm)
	theta_mrt_hum_ntr_is_n_pls := make([]float64, self.rms.n_rm)
	x_r_ntr_is_n_pls := make([]float64, self.rms.n_rm)
	for i := 0; i < self.rms.n_rm; i++ {
		if operation_mode_is_n[i] == STOP_OPEN {
			for j := 0; j < self.rms.n_rm; j++ {
				f_brm_is_is_n_pls.Set(i, j, f_brm_nv_is_is_n_pls.At(i, j))
				f_brm_ot_is_is_n_pls.Set(i, j, f_brm_ot_nv_is_is_n_pls.At(i, j))
				f_h_wgt_is_is_n.Set(i, j, f_h_wgt_nv_is_is_n.At(i, j))
			}
			v_vent_ntr_is_n[i] = self.rms.v_vent_ntr_set_is[i]
			f_brc_ot_is_n_pls[i] = f_brc_ot_nv_is_n_pls.AtVec(i)
			f_h_cst_is_n[i] = f_h_cst_nv_is_n.AtVec(i)
			theta_r_ot_ntr_is_n_pls[i] = theta_r_ot_ntr_nv_is_n_pls.AtVec(i)
			theta_r_ntr_is_n_pls[i] = theta_r_ntr_nv_is_n_pls.AtVec(i)
			theta_mrt_hum_ntr_is_n_pls[i] = theta_mrt_hum_ntr_nv_is_n_pls.AtVec(i)
			x_r_ntr_is_n_pls[i] = x_r_ntr_nv_is_n_pls.AtVec(i)
		} else {
			for j := 0; j < self.rms.n_rm; j++ {
				f_brm_is_is_n_pls.Set(i, j, f_brm_non_nv_is_is_n_pls.At(i, j))
				f_brm_ot_is_is_n_pls.Set(i, j, f_brm_ot_non_nv_is_is_n_pls.At(i, j))
				f_h_wgt_is_is_n.Set(i, j, f_h_wgt_non_nv_is_is_n.At(i, j))
			}
			v_vent_ntr_is_n[i] = 0.0
			f_brc_ot_is_n_pls[i] = f_brc_ot_non_nv_is_n_pls.AtVec(i)
			f_h_cst_is_n[i] = f_h_cst_non_nv_is_n.AtVec(i)
			theta_r_ot_ntr_is_n_pls[i] = theta_r_ot_ntr_non_nv_is_n_pls.AtVec(i)
			theta_r_ntr_is_n_pls[i] = theta_r_ntr_non_nv_is_n_pls.AtVec(i)
			theta_mrt_hum_ntr_is_n_pls[i] = theta_mrt_hum_ntr_non_nv_is_n_pls.AtVec(i)
			x_r_ntr_is_n_pls[i] = x_r_ntr_non_nv_is_n_pls.AtVec(i)
		}
	}

	theta_lower_target_is_n_pls, theta_upper_target_is_n_pls, _, _ :=
		self.op.get_theta_target_is_n(
			operation_mode_is_n,
			theta_r_ntr_is_n_pls,
			theta_mrt_hum_ntr_is_n_pls,
			x_r_ntr_is_n_pls,
			n,
			nn,
			self.es.is_radiative_heating_is,
			self.es.is_radiative_cooling_is,
			self.rms.met_is,
		)

	// ステップ n+1 における係数 f_flr, -, [j, i]
	f_flr_js_is_n := get_f_flr_js_is_n(
		self.es.f_flr_c_js_is,
		self.es.f_flr_h_js_is,
		operation_mode_is_n,
	)

	// ステップ n からステップ n+1 における室 i の放射暖冷房設備の対流成分比率, -, [i, 1]
	beta_is_n := get_beta_is_n(
		self.es.beta_c_is,
		self.es.beta_h_is,
		operation_mode_is_n,
	)

	// ステップ n における係数 f_FLB, K/W, [j, i]
	f_flb_js_is_n_pls := get_f_flb_js_is_n_pls(
		self.bs.a_s_js,
		beta_is_n,
		f_flr_js_is_n,
		self.bs.h_s_c_js,
		self.bs.h_s_r_js,
		self.bs.k_ei_js_js,
		self.bs.phi_a0_js,
		self.bs.phi_t0_js,
	)

	// ステップ n における係数 f_WSB, K/W, [j, i]
	f_wsb_js_is_n_pls := get_f_wsb_js_is_n_pls(
		f_flb_js_is_n_pls,
		ss.f_ax_js_js,
	)

	// ステップ n における係数 f_BRL, -, [i, i]
	f_brl_is_is_n := get_f_brl_is_is_n(
		self.bs.a_s_js,
		beta_is_n,
		f_wsb_js_is_n_pls,
		self.bs.h_s_c_js,
		self.bs.p_is_js,
	)

	// ステップn+1における室iの係数 f_XLR, K/W, [i, i]
	f_xlr_is_is_n_pls := get_f_xlr_is_is_n_pls(
		ss.f_mrt_hum_is_js,
		f_wsb_js_is_n_pls,
		ss.f_xot_is_is_n_pls,
		ss.k_r_is_n,
	)

	// ステップ n における係数 f_BRL_OT, -, [i, i]
	f_brl_ot_is_is_n := get_f_brl_ot_is_is_n(
		f_brl_is_is_n,
		f_brm_is_is_n_pls,
		f_xlr_is_is_n_pls,
	)

	// ステップ n+1 における室 i の作用温度, degree C, [i, 1] (ステップn+1における瞬時値）
	// ステップ n における室 i に設置された対流暖房の放熱量, W, [i, 1] (ステップn～ステップn+1までの平均値）
	// ステップ n における室 i に設置された放射暖房の放熱量, W, [i, 1]　(ステップn～ステップn+1までの平均値）
	theta_ot_is_n_pls, l_cs_is_n, l_rs_is_n := get_next_temp_and_load(
		self.scd.ac_demand_is_ns,
		f_brc_ot_is_n_pls,
		f_brm_ot_is_is_n_pls,
		f_brl_ot_is_is_n,
		theta_lower_target_is_n_pls,
		theta_upper_target_is_n_pls,
		operation_mode_is_n,
		self.es.is_radiative_heating_is,
		self.es.is_radiative_cooling_is,
		self.es.q_rs_h_max_is,
		self.es.q_rs_c_max_is,
		theta_r_ot_ntr_is_n_pls,
		n,
	)

	// ステップ n+1 における室 i の室温, degree C, [i, 1]
	theta_r_is_n_pls := get_theta_r_is_n_pls(
		_f_xc_is_n_pls,
		f_xlr_is_is_n_pls,
		ss.f_xot_is_is_n_pls,
		l_rs_is_n,
		theta_ot_is_n_pls,
	)

	// ステップ n+1 における境界 j の表面温度, degree C, [j, 1]
	theta_s_js_n_pls := get_theta_s_js_n_pls(
		f_wsb_js_is_n_pls,
		ss.f_wsc_js_ns.ColView(nn_plus+1),
		ss.f_wsr_js_is,
		f_wsv_js_n_pls,
		l_rs_is_n,
		theta_r_is_n_pls,
	)

	// ステップ n+1 における室 i　の備品等の温度, degree C, [i, 1]
	// TODO: q_sol_frt_is_ns の値は n+1 の値を使用するべき？
	theta_frt_is_n_pls := get_theta_frt_is_n_pls(
		self.rms.c_sh_frt_is,
		delta_t,
		self.rms.g_sh_frt_is,
		ss.q_sol_frt_is_ns.(mat.ColViewer).ColView(nn_plus),
		c_n.theta_frt_is_n,
		theta_r_is_n_pls,
	)

	// ステップ n+1 における室 i の人体に対する平均放射温度, degree C, [i, 1]
	theta_mrt_hum_is_n_pls := get_theta_mrt_hum_is_n_pls(
		ss.f_mrt_hum_is_js,
		theta_s_js_n_pls,
	)

	// ステップ n+1 における境界 j の等価温度, degree C, [j, 1]
	theta_ei_js_n_pls := get_theta_ei_js_n_pls(
		self.bs.a_s_js,
		beta_is_n,
		ss.f_mrt_is_js,
		f_flr_js_is_n,
		self.bs.h_s_c_js,
		self.bs.h_s_r_js,
		l_rs_is_n,
		self.bs.p_js_is,
		ss.q_s_sol_js_ns.(mat.ColViewer).ColView(nn_plus+1),
		theta_r_is_n_pls,
		theta_s_js_n_pls,
	)

	// ステップ n+1 における境界 j の裏面温度, degree C, [j, 1]
	// TODO: この値は記録にしか使用していないので、ポスト処理にまわせる。
	// theta_rear_js_n_pls := get_theta_s_rear_js_n(
	// 	self.bs.k_ei_js_js,
	// 	theta_ei_js_n_pls,
	// 	self.bs.k_eo_js,
	// 	self.bs.theta_o_eqv_js_ns.ColView(nn+1),
	// 	self.bs.k_s_r_js,
	// 	theta_r_is_n_pls,
	// )

	// ステップ n+1 における境界 j の表面熱流（壁体吸熱を正とする）, W/m2, [j, 1]
	q_s_js_n_pls := get_q_s_js_n_pls(
		self.bs.h_s_c_js,
		self.bs.h_s_r_js,
		theta_ei_js_n_pls,
		theta_s_js_n_pls,
	)

	// ステップ n+1 における室 i∗ の絶対湿度がステップ n から n+1 における室 i の潜熱負荷に与える影響を表す係数, kg/(s (kg/kg(DA))), [i, i*]
	// ステップ n から n+1 における室 i の潜熱負荷に与える影響を表す係数, kg/s, [i, 1]
	f_l_cl_cst_is_n, f_l_cl_wgt_is_is_n := self.get_f_l_cl(
		l_cs_is_n,
		theta_r_is_n_pls,
		x_r_ntr_is_n_pls,
	)

	// ステップ n+1 における室 i の 絶対湿度, kg/kg(DA), [i, 1]
	x_r_is_n_pls := get_x_r_is_n_pls(
		f_h_cst_is_n,
		f_h_wgt_is_is_n,
		f_l_cl_cst_is_n,
		f_l_cl_wgt_is_is_n,
	)

	// ステップ n から ステップ n+1 における室 i の潜熱負荷（加湿を正・除湿を負とする）, W, [i, 1]
	// l_cl_is_n := get_l_cl_is_n(
	// 	f_l_cl_wgt_is_is_n,
	// 	f_l_cl_cst_is_n,
	// 	get_l_wtr(),
	// 	x_r_is_n_pls,
	// )

	// ステップ n+1 における室 i の備品等等の絶対湿度, kg/kg(DA), [i, 1]
	x_frt_is_n_pls := get_x_frt_is_n_pls(
		self.rms.c_lh_frt_is,
		delta_t,
		self.rms.g_lh_frt_is,
		c_n.x_frt_is_n,
		x_r_is_n_pls,
	)

	// if recorder is not None:
	//     recorder.recording(
	//         n=n,
	//         theta_r_is_n_pls=theta_r_is_n_pls,
	//         theta_mrt_hum_is_n_pls=theta_mrt_hum_is_n_pls,
	//         x_r_is_n_pls=x_r_is_n_pls,
	//         theta_frt_is_n_pls=theta_frt_is_n_pls,
	//         x_frt_is_n_pls=x_frt_is_n_pls,
	//         theta_ei_js_n_pls=theta_ei_js_n_pls,
	//         q_s_js_n_pls=q_s_js_n_pls,
	//         theta_ot_is_n_pls=theta_ot_is_n_pls,
	//         theta_s_js_n_pls=theta_s_js_n_pls,
	//         theta_rear_js_n=theta_rear_js_n_pls,
	//         f_cvl_js_n_pls=f_cvl_js_n_pls,
	//         operation_mode_is_n=operation_mode_is_n,
	//         l_cs_is_n=l_cs_is_n,
	//         l_rs_is_n=l_rs_is_n,
	//         l_cl_is_n=l_cl_is_n,
	//         h_hum_c_is_n=h_hum_c_is_n,
	//         h_hum_r_is_n=h_hum_r_is_n,
	//         q_hum_is_n=q_hum_is_n,
	//         x_hum_is_n=x_hum_is_n,
	//         v_leak_is_n=v_leak_is_n,
	//         v_vent_ntr_is_n=v_vent_ntr_is_n
	//     )

	c_n.operation_mode_is_n = operation_mode_is_n
	c_n.theta_r_is_n = theta_r_is_n_pls
	c_n.theta_mrt_hum_is_n = theta_mrt_hum_is_n_pls
	c_n.x_r_is_n = x_r_is_n_pls
	c_n.theta_dsh_srf_a_js_ms_n = theta_dsh_s_a_js_ms_n_pls
	c_n.theta_dsh_srf_t_js_ms_n = theta_dsh_s_t_js_ms_n_pls
	c_n.q_s_js_n = q_s_js_n_pls
	c_n.theta_frt_is_n = theta_frt_is_n_pls
	c_n.x_frt_is_n = x_frt_is_n_pls
	c_n.theta_ei_js_n = theta_ei_js_n_pls

	return c_n
	// return NewConditions(
	// 	operation_mode_is_n,
	// 	theta_r_is_n_pls,
	// 	theta_mrt_hum_is_n_pls,
	// 	x_r_is_n_pls,
	// 	theta_dsh_s_a_js_ms_n_pls,
	// 	theta_dsh_s_t_js_ms_n_pls,
	// 	q_s_js_n_pls,
	// 	theta_frt_is_n_pls,
	// 	x_frt_is_n_pls,
	// 	theta_ei_js_n_pls,
	// )

}

// 地盤の計算（n+1ステップを計算する）
func _run_tick_ground(self *Sequence, pp *PreCalcParameters, gc_n *GroundConditions, n int, nn int) *GroundConditions {
	if gc_n == nil {
		return nil
	}

	is_ground := self.bs.is_ground_js
	ground_idx := make([]int, 0, self.bs.n_b)
	for i, v := range is_ground {
		if v {
			ground_idx = append(ground_idx, i)
		}
	}

	_, l := self.bs.theta_o_eqv_js_ns.Dims()

	theta_o_eqv_js_ns := mat.NewDense(len(ground_idx), l, nil)
	for i := 0; i < len(ground_idx); i++ {
		theta_o_eqv_js_ns.SetRow(i, self.bs.theta_o_eqv_js_ns.RawRowView(ground_idx[i]))
	}

	h_i_js := mat.NewVecDense(len(ground_idx), nil)
	for i := 0; i < len(ground_idx); i++ {
		gidx := ground_idx[i]
		h_i_js.SetVec(i, self.bs.h_s_r_js.AtVec(gidx)+self.bs.h_s_c_js.AtVec(gidx))
	}

	theta_dsh_srf_a_js_ms_npls := make([][]float64, len(ground_idx))
	for j := 0; j < len(ground_idx); j++ {
		gidx := ground_idx[j]
		r := self.bs.r_js_ms[gidx]
		phi_a1 := self.bs.phi_a1_js_ms[gidx]
		theta_dsh_srf_a_js_ms_npls[j] = make([]float64, 12)
		for i := 0; i < 12; i++ {
			theta_dsh_srf_a_js_ms_npls[j][i] = phi_a1[i]*gc_n.q_srf_js_n[j] +
				r[i]*gc_n.theta_dsh_srf_a_js_ms_n[j][i]
		}
	}

	theta_dsh_srf_t_js_ms_npls := make([][]float64, len(ground_idx))
	for j := 0; j < len(ground_idx); j++ {
		gidx := ground_idx[j]
		phi_t1 := self.bs.phi_t1_js_ms[gidx]
		r := self.bs.r_js_ms[gidx]
		theta_dsh_srf_t_js_ms_npls[j] = make([]float64, 12)
		for i := 0; i < 12; i++ {
			theta_dsh_srf_t_js_ms_npls[j][i] =
				phi_t1[i]*self.bs.k_eo_js.AtVec(gidx)*self.bs.theta_o_eqv_js_ns.At(gidx, nn) +
					r[i]*gc_n.theta_dsh_srf_t_js_ms_n[j][i]
		}
	}

	theta_s_js_npls := make([]float64, len(ground_idx))
	for j := 0; j < len(ground_idx); j++ {
		sum_a := floats.Sum(theta_dsh_srf_a_js_ms_npls[j])
		sum_t := floats.Sum(theta_dsh_srf_t_js_ms_npls[j])

		gidx := ground_idx[j]
		theta_s_js_npls[j] =
			(self.bs.phi_a0_js.AtVec(gidx)*h_i_js.AtVec(j)*self.weather.theta_o_ns_plus[nn+1] +
				self.bs.phi_t0_js.AtVec(gidx)*self.bs.k_eo_js.AtVec(gidx)*self.bs.theta_o_eqv_js_ns.At(gidx, nn+1) +
				sum_a + sum_t) / (1.0 + self.bs.phi_a0_js.AtVec(gidx)*h_i_js.AtVec(j))
	}

	q_srf_js_n := make([]float64, len(ground_idx))
	for j := 0; j < len(ground_idx); j++ {
		q_srf_js_n[j] = h_i_js.AtVec(j) * (self.weather.theta_o_ns_plus[nn+1] - theta_s_js_npls[j])
	}

	return &GroundConditions{
		theta_dsh_srf_a_js_ms_n: theta_dsh_srf_a_js_ms_npls,
		theta_dsh_srf_t_js_ms_n: theta_dsh_srf_t_js_ms_npls,
		q_srf_js_n:              q_srf_js_n,
	}
}

/*
Args:

	c_lh_frt_is: 室 i の備品等の湿気容量, kg/(kg/kg(DA)), [i, 1]
	delta_t: 1ステップの時間間隔, s
	g_lh_frt_is: 室 i の備品等と空気間の湿気コンダクタンス, kg/(s kg/kg(DA)), [i, 1]
	x_frt_is_n: ステップ n における室 i の備品等の絶対湿度, kg/kg(DA), [i, 1]
	x_r_is_n_pls: ステップ n+1 における室 i の絶対湿度, kg/kg(DA), [i, 1]

Returns:

	ステップ n+1 における室 i の備品等等の絶対湿度, kg/kg(DA), [i, 1]

Notes:

	式(1.1)
*/
func get_x_frt_is_n_pls(
	c_lh_frt_is mat.Vector,
	delta_t float64,
	g_lh_frt_is mat.Vector,
	x_frt_is_n []float64,
	x_r_is_n_pls []float64,
) []float64 {
	var temp1, temp2, temp3 mat.VecDense

	// c_lh_frt_is * x_frt_is_n
	temp1.MulElemVec(c_lh_frt_is, mat.NewVecDense(len(x_frt_is_n), x_frt_is_n))

	// g_lh_frt_is * x_r_is_n_pls
	temp2.MulElemVec(g_lh_frt_is, mat.NewVecDense(len(x_r_is_n_pls), x_r_is_n_pls))

	// c_lh_frt_is * x_frt_is_n + delta_t * (g_lh_frt_is * x_r_is_n_pls)
	temp1.AddScaledVec(&temp1, delta_t, &temp2)

	// c_lh_frt_is + delta_t * g_lh_frt_is
	temp3.AddScaledVec(c_lh_frt_is, delta_t, g_lh_frt_is)

	// temp1 / temp3
	temp1.DivElemVec(&temp1, &temp3)

	slice := make([]float64, temp1.Len())
	copy(slice, temp1.RawVector().Data)

	return slice
}

/*
Args:

	f_l_cl_wgt_is_is_n: 係数, kg/(s (kg/kg(DA)))
	f_l_cl_cst_is_n: 係数, kg/s
	l_wtr: 水の蒸発潜熱, J/kg
	x_r_is_n_pls: ステップ n+1 における室 i の絶対湿度, kg/kg(DA)

Returns:

	ステップ n から ステップ n+1 における室 i の潜熱負荷（加湿を正・除湿を負とする）, W

Notes:

	式(1.2)
*/
func get_l_cl_is_n(
	f_l_cl_wgt_is_is_n mat.Matrix,
	f_l_cl_cst_is_n mat.Vector,
	l_wtr float64,
	x_r_is_n_pls mat.Vector,
) *mat.Dense {
	var temp1 mat.Dense
	temp1.Mul(f_l_cl_wgt_is_is_n, x_r_is_n_pls)
	temp1.Add(&temp1, f_l_cl_cst_is_n)
	temp1.Scale(l_wtr, &temp1)

	return &temp1
}

/*
Args:

	f_h_cst_is_n: 係数, kg/s
	f_h_wgt_is_is_n: 係数, kg/(s (kg/kg(DA)))
	f_l_cl_cst_is_n: 係数, kg/s
	f_l_cl_wgt_is_is_n: 係数, kg/(s (kg/kg(DA)))

Returns:

	ステップ n+1 における室 i の 絶対湿度, kg/kg(DA), [i, 1]

Notes:

	式(1.3)
*/
func get_x_r_is_n_pls(
	f_h_cst_is_n []float64,
	f_h_wgt_is_is_n mat.Matrix,
	f_l_cl_cst_is_n mat.Vector,
	f_l_cl_wgt_is_is_n mat.Matrix,
) []float64 {
	//f_h_cst_is_n + f_l_cl_cst_is_n
	var temp1 mat.VecDense
	temp1.AddVec(mat.NewVecDense(len(f_h_cst_is_n), f_h_cst_is_n), f_l_cl_cst_is_n)

	// f_h_wgt_is_is_n - f_l_cl_wgt_is_is_n
	var temp2 mat.Dense
	temp2.Sub(f_h_wgt_is_is_n, f_l_cl_wgt_is_is_n)

	var result mat.VecDense
	result.SolveVec(&temp2, &temp1)

	slice := make([]float64, result.Len())
	copy(slice, result.RawVector().Data)

	return slice
}

/*
Args:

	f_h_cst_non_nv_is_n: ステップnにおける自然風非利用時の室iの係数f_h_cst, kg/s
	f_h_wgt_non_nv_is_is_n: ステップnにおける自然風利用時の室iの係数f_h_wgt, kg/(s (kg/kg(DA)))
	f_h_cst_nv_is_n: ステップnにおける自然風利用時の室iの係数f_h_cst, kg/s
	f_h_wgt_nv_is_is_n: ステップnにおける自然風利用時の室iの係数f_wgt, kg/(s (kg/kg(DA)))

Returns:

	ステップ n+1 における室 i の加湿・除湿を行わない場合の絶対湿度, kg/kg(DA) [i, 1]

Notes:

	式(1.4)
*/
func get_x_r_ntr_is_n_pls(
	f_h_cst_non_nv_is_n mat.Vector,
	f_h_wgt_non_nv_is_is_n mat.Matrix,
	f_h_cst_nv_is_n mat.Vector,
	f_h_wgt_nv_is_is_n mat.Matrix,
) (*mat.VecDense, *mat.VecDense) {
	var result1, result2 mat.VecDense
	result1.SolveVec(f_h_wgt_non_nv_is_is_n, f_h_cst_non_nv_is_n)
	result2.SolveVec(f_h_wgt_nv_is_is_n, f_h_cst_nv_is_n)
	return &result1, &result2
}

/*
Args:

	c_lh_frt_is: 室 i の備品等の湿気容量, kg/(kg/kg(DA)), [i, 1]
	delta_t: 1ステップの時間間隔, s
	g_lh_frt_is: 室 i の備品等と空気間の湿気コンダクタンス, kg/(s kg/kg(DA)), [i, 1]
	v_rm_is: 室 i の容量, m3, [i, 1]
	v_vent_int_is_is_n:　ステップ n から ステップ n+1 における室 i* から室 i への室間の空気移動量（流出換気量を含む）, m3/s
	v_vent_out_is_n: ステップ n から ステップ n+1 における室 i の換気・すきま風・自然風の利用による外気の流入量, m3/s

Returns:

	ステップ n における室 i* の絶対湿度が室 i の潜熱バランスに与える影響を表す係数,　kg/(s kg/kg(DA)), [i, i]

Notes:

	式(1.5)
*/
func get_f_h_wgt_is_is_n(
	c_lh_frt_is mat.Vector,
	delta_t float64,
	g_lh_frt_is mat.Vector,
	v_rm_is []float64,
	v_vent_int_is_is_n mat.Matrix,
	v_vent_out_is_n []float64,
	v_vent_ntr_is []float64,
) (*mat.Dense, *mat.Dense) {

	// rho_a*(v_rm_is/delta_t+v_vent_out_is_n)
	_temp1 := make([]float64, len(v_rm_is))
	floats.AddScaledTo(_temp1, v_vent_out_is_n, 1/delta_t, v_rm_is)
	floats.Scale(rho_a, _temp1)
	temp1 := mat.NewVecDense(len(_temp1), _temp1)

	// c_lh_frt_is + delta_t * g_lh_frt_is)
	var temp2 mat.VecDense
	temp2.AddScaledVec(c_lh_frt_is, delta_t, g_lh_frt_is)

	// c_lh_frt_is * g_lh_frt_is / (c_lh_frt_is + delta_t * g_lh_frt_is)
	var temp3 mat.VecDense
	temp3.MulElemVec(c_lh_frt_is, g_lh_frt_is)
	temp3.DivElemVec(&temp3, &temp2)

	// diag of temp1 + temp3
	temp1.AddVec(temp1, &temp3)
	temp4 := mat.NewDiagDense(temp1.Len(), nil)
	for i := 0; i < temp1.Len(); i++ {
		temp4.SetDiag(i, temp1.AtVec(i))
	}

	var result1 mat.Dense
	result1.Apply(func(i, j int, v float64) float64 {
		return v - rho_a*v_vent_int_is_is_n.At(i, j)
	}, temp4)

	var result2, temp5 mat.Dense
	temp5.Scale(rho_a, mat.NewDiagDense(len(v_vent_ntr_is), v_vent_ntr_is))
	result2.Add(&result1, &temp5)

	return &result1, &result2
}

/*
Args:

	c_lh_frt_is: 室 i の備品等の湿気容量, kg/(kg/kg(DA)), [i, 1]
	delta_t: 1ステップの時間間隔, s
	g_lh_frt_is: 室 i の備品等と空気間の湿気コンダクタンス, kg/(s kg/kg(DA)), [i, 1]
	rho_a: 空気の密度, kg/m3
	v_rm_is: 室 i の容量, m3, [i, 1]
	x_frt_is_n: ステップ n における室 i の備品等の絶対湿度, kg/kg(DA), [i, 1]
	x_gen_is_n: ステップ n からステップ n+1 における室 i の人体発湿を除く内部発湿, kg/s
	x_hum_is_n: ステップ n からステップ n+1 における室 i の人体発湿, kg/s
	x_o_n_pls: ステップ n における外気絶対湿度, kg/kg(DA)
	x_r_is_n: ステップ n における室 i の絶対湿度, kg/kg(DA)
	v_vent_out_non_nv_is_n: ステップnからステップn+1における室iの換気・隙間風による外気の流入量, m3/s, [i, 1]
	v_vent_ntr_is: 室iの自然風利用時の換気量, m3/s, [i, 1]

Returns:
ステップnにおける室iの自然風の非利用時の潜熱バランスに関する係数f_h_cst, kg/s, [i, 1]
ステップnにおける室iの自然風の利用時の潜熱バランスに関する係数f_h_cst, kg/s, [i, 1]

Notes:

	式(1.6)
*/
func get_f_h_cst_is_n(
	c_lh_frt_is mat.Vector,
	delta_t float64,
	g_lh_frt_is mat.Vector,
	rho_a float64,
	v_rm_is []float64,
	x_frt_is_n []float64,
	x_gen_is_n []float64,
	x_hum_is_n []float64,
	x_o_n_pls float64,
	x_r_is_n []float64,
	v_vent_out_non_nv_is_n []float64,
	v_vent_ntr_is []float64,
) (*mat.VecDense, *mat.VecDense) {
	// Python: rho_a * v_rm_is / delta_t * x_r_is_n
	temp1 := make([]float64, len(v_rm_is))
	floats.ScaleTo(temp1, rho_a/delta_t, v_rm_is)
	floats.Mul(temp1, x_r_is_n)

	// Python: rho_a * v_vent_out_non_nv_is_n * x_o_n_pls
	temp2 := make([]float64, len(v_rm_is))
	floats.ScaleTo(temp2, rho_a*x_o_n_pls, v_vent_out_non_nv_is_n)

	// Python: c_lh_frt_is * g_lh_frt_is / (c_lh_frt_is + delta_t * g_lh_frt_is) * x_frt_is_n
	temp3 := &mat.VecDense{}
	temp4 := &mat.VecDense{}
	temp3.MulElemVec(c_lh_frt_is, g_lh_frt_is)            //(c_lh_frt_is * g_lh_frt_is)
	temp4.AddScaledVec(c_lh_frt_is, delta_t, g_lh_frt_is) //(c_lh_frt_is + delta_t * g_lh_frt_is)
	temp3.DivElemVec(temp3, temp4)                        //(c_lh_frt_is * g_lh_frt_is / (c_lh_frt_is + delta_t * g_lh_frt_is)
	temp3.MulElemVec(temp3, mat.NewVecDense(len(x_frt_is_n), x_frt_is_n))

	// Python: temp1 + temp2 + temp3 + x_gen_is_n + x_hum_is_n
	var result1 mat.VecDense
	result1.AddVec(mat.NewVecDense(len(temp1), temp1), mat.NewVecDense(len(temp2), temp2))
	result1.AddVec(&result1, temp3)
	result1.AddVec(&result1, mat.NewVecDense(len(x_gen_is_n), x_gen_is_n))
	result1.AddVec(&result1, mat.NewVecDense(len(x_hum_is_n), x_hum_is_n))

	var result2 mat.VecDense
	result2.AddScaledVec(&result1, rho_a*x_o_n_pls, mat.NewVecDense(len(v_vent_ntr_is), v_vent_ntr_is))

	return &result1, &result2
}

/*
Args:

	n_hum_is_n: ステップ n からステップ n+1 における室 i の在室人数, -
	x_hum_psn_is_n: ステップ n からステップ n+1 における室 i の1人あたりの人体発湿, kg/s

Returns:

	ステップnの室iにおける人体発湿, kg/s, [i, 1]

Notes:

	式(1.7)
*/
func get_x_hum_is_n(n_hum_is_n []float64, x_hum_psn_is_n []float64) []float64 {
	x_hum_is_n := make([]float64, len(x_hum_psn_is_n))
	floats.MulTo(x_hum_is_n, n_hum_is_n, x_hum_psn_is_n)
	return x_hum_is_n
}

/*
Args:

	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	h_s_r_js: 境界 j の室内側放射熱伝達率, W/(m2 K), [j, 1]
	theta_ei_js_n_pls: ステップ n+1 における境界 j の等価温度, degree C, [j, 1]
	theta_s_js_n_pls: ステップ n+1 における境界 j の表面温度, degree C, [j, 1]

Returns:

	ステップ n+1 における境界 j の表面熱流（壁体吸熱を正とする）, W/m2, [j, 1]

Notes:

	式(2.1)
*/
func get_q_s_js_n_pls(
	h_s_c_js mat.Vector,
	h_s_r_js mat.Vector,
	theta_ei_js_n_pls []float64,
	theta_s_js_n_pls mat.Vector,
) []float64 {
	var temp1, temp2 mat.VecDense
	temp1.SubVec(mat.NewVecDense(len(theta_ei_js_n_pls), theta_ei_js_n_pls), theta_s_js_n_pls)
	temp2.AddVec(h_s_c_js, h_s_r_js)
	temp1.MulElemVec(&temp1, &temp2)

	slice := make([]float64, temp1.Len())
	copy(slice, temp1.RawVector().Data)
	return slice
}

/*
Args:

	a_s_js: 境界 j の面積, m2, [j, 1]
	beta_is_n: ステップ n からステップ n+1 における室 i の放射暖冷房設備の対流成分比率, -, [i, 1]
	f_mrt_is_js: 室 i の微小球に対する境界 j の形態係数, -, [i, j]
	f_flr_js_is_n: ステップ n からステップ n+1 における室 i の放射暖冷房設備の放熱量の放射成分に対する境界 j の室内側表面の吸収比率, -, [j, i]
	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	h_s_r_js: 境界 j の室内側放射熱伝達率, W/(m2 K), [j, 1]
	l_rs_is_n: ステップ n からステップ n+1 における室 i に放射暖冷房設備の顕熱処理量（暖房を正・冷房を負とする）, W, [i, 1]
	p_js_is: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [j, i]
	q_s_sol_js_n_pls: ステップ n+1 における境界 j の透過日射吸収熱量, W/m2, [j, 1]
	theta_r_is_n_pls: ステップ n+1 における室 i の温度, degree C, [i, 1]
	theta_s_js_n_pls: ステップ n+1 における境界 j の表面温度, degree C, [j, 1]

Returns:

	ステップ n+1 における境界 j の等価温度, degree C, [j, 1]

Notes:

	式(2.2)
*/
func get_theta_ei_js_n_pls(
	a_s_js mat.Vector,
	beta_is_n mat.Vector,
	f_mrt_is_js mat.Matrix,
	f_flr_js_is_n mat.Matrix,
	h_s_c_js mat.Vector,
	h_s_r_js mat.Vector,
	l_rs_is_n []float64,
	p_js_is mat.Matrix,
	q_s_sol_js_n_pls mat.Vector,
	theta_r_is_n_pls []float64,
	theta_s_js_n_pls mat.Vector,
) []float64 {
	// Python: h_s_c_js*np.dot(p_js_is, theta_r_is_n_pls)
	var temp1 mat.VecDense
	temp1.MulVec(p_js_is, mat.NewVecDense(len(theta_r_is_n_pls), theta_r_is_n_pls))
	temp1.MulElemVec(h_s_c_js, &temp1)

	// Python: h_s_r_js*np.dot(np.dot(p_js_is, f_mrt_is_js), theta_s_js_n_pls)
	var temp2 mat.Dense
	var temp3 mat.VecDense
	temp2.Mul(p_js_is, f_mrt_is_js)
	temp3.MulVec(&temp2, theta_s_js_n_pls)
	temp3.MulElemVec(h_s_r_js, &temp3)

	// Python: np.dot(f_flr_js_is_n, (1.0-beta_is_n)*l_rs_is_n)/a_s_js
	var temp4, temp5 mat.VecDense
	temp4.MulElemVec(beta_is_n, mat.NewVecDense(len(l_rs_is_n), l_rs_is_n))
	temp4.SubVec(mat.NewVecDense(len(l_rs_is_n), l_rs_is_n), &temp4)
	temp5.MulVec(f_flr_js_is_n, &temp4)
	temp5.DivElemVec(&temp5, a_s_js)

	// Python: (temp1 + temp3 + q_s_sol_js_n_pls + temp5) / (h_s_c_js + h_s_r_js)
	var result mat.VecDense
	result.AddVec(&temp1, &temp3)
	result.AddVec(&result, q_s_sol_js_n_pls)
	result.AddVec(&result, &temp5)
	var temp6 mat.VecDense
	temp6.AddVec(h_s_c_js, h_s_r_js)
	result.DivElemVec(&result, &temp6)

	slice := make([]float64, result.Len())
	copy(slice, result.RawVector().Data)

	return slice
}

/*
Args:

	f_mrt_hum_is_js: 室 i の人体に対する境界 j の形態係数, -, [i, j]
	theta_s_js_n_pls: ステップ n+1 における境界 j の表面温度, degree C, [j, 1]

Returns:

	ステップ n+1 における室 i の人体に対する平均放射温度, degree C, [i, 1]

Notes:

	式(2.3)
*/
func get_theta_mrt_hum_is_n_pls(
	f_mrt_hum_is_js mat.Matrix,
	theta_s_js_n_pls mat.Vector,
) []float64 {
	var result mat.VecDense
	result.MulVec(f_mrt_hum_is_js, theta_s_js_n_pls)
	slice := make([]float64, result.Len())
	copy(slice, result.RawVector().Data)
	return slice
}

/*
Args:

	c_sh_frt_is: 室 i の備品等の熱容量, J/K, [i, 1]
	delta_t: 1ステップの時間間隔, s
	g_sh_frt_is: 室 i の備品等と空気間の熱コンダクタンス, W/K, [i, 1]
	q_sol_frt_is_n: ステップ n からステップ n+1 における室 i に設置された備品等による透過日射吸収熱量時間平均値, W, [i, 1]
	theta_frt_is_n: ステップ n における室 i の備品等の温度, degree C, [i, 1]
	theta_r_is_n_pls: ステップ n+1 における室 i の温度, degree C, [i, 1]

Returns:

	ステップ n+1 における室 i　の備品等の温度, degree C, [i, 1]

Notes:

	式(2.4)
*/
func get_theta_frt_is_n_pls(
	c_sh_frt_is mat.Vector,
	delta_t float64,
	g_sh_frt_is mat.Vector,
	q_sol_frt_is_n mat.Vector,
	theta_frt_is_n []float64,
	theta_r_is_n_pls []float64,
) []float64 {
	var temp1, temp2, temp3 mat.VecDense

	// c_sh_frt_is*theta_frt_is_n
	temp1.MulElemVec(c_sh_frt_is, mat.NewVecDense(len(theta_frt_is_n), theta_frt_is_n))

	// delta_t*g_sh_frt_is*theta_r_is_n_pls
	temp2.MulElemVec(g_sh_frt_is, mat.NewVecDense(len(theta_r_is_n_pls), theta_r_is_n_pls))
	temp2.ScaleVec(delta_t, &temp2)

	// q_sol_frt_is_n*delta_t
	temp3.ScaleVec(delta_t, q_sol_frt_is_n)

	// (c_sh_frt_is*theta_frt_is_n + delta_t*g_sh_frt_is*theta_r_is_n_pls + q_sol_frt_is_n*delta_t)
	temp1.AddVec(&temp1, &temp2)
	temp1.AddVec(&temp1, &temp3)

	// c_sh_frt_is + delta_t*g_sh_frt_is
	temp2.ScaleVec(delta_t, g_sh_frt_is)
	temp2.AddVec(c_sh_frt_is, &temp2)

	temp1.DivElemVec(&temp1, &temp2)

	slice := make([]float64, temp1.Len())
	copy(slice, temp1.RawVector().Data)

	return slice
}

/*
Args:

	f_wsb_js_is_n_pls: ステップ n+1 における係数 f_WSB, K/W, [j, 1]
	f_wsc_js_n_pls: ステップ n+1 における係数 f_WSC, degree C, [j, 1]
	f_wsr_js_is: 係数 f_WSR, - [j, i]
	f_wsv_js_n_pls: ステップ n+1 における係数 f_WSV, degree C, [j, 1]
	l_rs_is_n: ステップ n からステップ n+1 における室 i に放射暖冷房設備の顕熱処理量（暖房を正・冷房を負とする）, W, [i, 1]
	theta_r_is_n_pls: ステップ n+1 における室 i の温度, degree C, [i, 1]

Returns:

	ステップ n+1 における境界 j の表面温度, degree C, [j, 1]

Notes:

	式(2.5)
*/
func get_theta_s_js_n_pls(
	f_wsb_js_is_n_pls mat.Matrix,
	f_wsc_js_n_pls mat.Vector,
	f_wsr_js_is mat.Matrix,
	f_wsv_js_n_pls mat.Vector,
	l_rs_is_n []float64,
	theta_r_is_n_pls []float64,
) *mat.VecDense {
	var temp1, temp2 mat.VecDense
	temp1.MulVec(f_wsr_js_is, mat.NewVecDense(len(theta_r_is_n_pls), theta_r_is_n_pls)) //np.dot(f_wsr_js_is, theta_r_is_n_pls)
	temp2.MulVec(f_wsb_js_is_n_pls, mat.NewVecDense(len(l_rs_is_n), l_rs_is_n))         //np.dot(f_wsb_js_is_n_pls, l_rs_is_n)

	temp1.AddVec(&temp1, f_wsc_js_n_pls)
	temp1.AddVec(&temp1, f_wsv_js_n_pls)
	temp1.AddVec(&temp1, &temp2)

	return &temp1
}

/*
Args:

	f_xc_is_n_pls: ステップ n+1 における係数 f_XC, degree C, [i, 1]
	f_xlr_is_is_n_pls: ステップ n+1 における係数 f_XLR, K/W, [i, i]
	f_xot_is_is_n_pls: ステップ n+1 における係数 f_XOT, -, [i, i]
	l_rs_is_n: ステップ n からステップ n+1 における室 i に放射暖冷房設備の顕熱処理量（暖房を正・冷房を負とする）, W, [i, 1]
	theta_ot_is_n_pls: ステップ n+1 における室 i の作用温度, ℃

Returns:

	ステップ n+1 における室 i の室温, degree C, [i, 1]

Notes:

	式(2.6)
*/
func get_theta_r_is_n_pls(
	f_xc_is_n_pls mat.Vector,
	f_xlr_is_is_n_pls mat.Matrix,
	f_xot_is_is_n_pls mat.Matrix,
	l_rs_is_n []float64,
	theta_ot_is_n_pls []float64,
) []float64 {
	// Python: np.dot(f_xot_is_is_n_pls, theta_ot_is_n_pls)
	var temp1 mat.VecDense
	temp1.MulVec(f_xot_is_is_n_pls, mat.NewVecDense(len(theta_ot_is_n_pls), theta_ot_is_n_pls))

	// Python: np.dot(f_xlr_is_is_n_pls, l_rs_is_n)
	var temp2 mat.VecDense
	temp2.MulVec(f_xlr_is_is_n_pls, mat.NewVecDense(len(l_rs_is_n), l_rs_is_n))

	// Python: temp1 - temp2 - f_xc_is_n_pls
	var temp3 mat.VecDense
	temp3.SubVec(&temp1, &temp2)
	temp3.SubVec(&temp3, f_xc_is_n_pls)

	slice := make([]float64, temp3.Len())
	copy(slice, temp3.RawVector().Data)

	return slice
}

/*
Args:

	f_brl_is_is_n: ステップ n における係数 f_BRL, -, [i, i]
	f_brm_is_is_n_pls: ステップ n+1 における係数 f_BRM, W/K, [i, i]
	f_xlr_is_is_n_pls: ステップ n+1 における係数 f_XLR, K/W, [i, i]

Returns:

	ステップ n における係数 f_BRL,OT, -, [i, i]

Notes:

	式(2.8)
*/
func get_f_brl_ot_is_is_n(
	f_brl_is_is_n mat.Matrix,
	f_brm_is_is_n_pls mat.Matrix,
	f_xlr_is_is_n_pls mat.Matrix,
) *mat.Dense {
	var result mat.Dense
	result.Mul(f_brm_is_is_n_pls, f_xlr_is_is_n_pls)
	result.Add(f_brl_is_is_n, &result)
	return &result
}

/*
Args:

	f_mrt_hum_is_js: 室 i の人体に対する境界 j の形態係数, -, [i, j]
	f_wsb_js_is_n_pls: ステップ n+1 における係数 f_WSB, K/W, [j, 1]
	f_xot_is_is_n_pls: ステップ n+1 における係数 f_XOT, -, [i, i]
	k_r_is_n: ステップ n における室 i の人体表面の放射熱伝達率が総合熱伝達率に占める割合, -, [i, 1]

Returns:

	ステップ n+1 における係数 f_XLR, K/W, [i, i]

Notes:

	式(2.9)
*/
func get_f_xlr_is_is_n_pls(
	f_mrt_hum_is_js mat.Matrix,
	f_wsb_js_is_n_pls mat.Matrix,
	f_xot_is_is_n_pls mat.Matrix,
	k_r_is_n mat.Vector,
) *mat.Dense {
	// Python: np.dot(f_mrt_hum_is_js, f_wsb_js_is_n_pls)
	var temp1 mat.Dense
	temp1.Mul(f_mrt_hum_is_js, f_wsb_js_is_n_pls)

	// Python: k_r_is_n * temp1
	temp1.Apply(func(i, j int, v float64) float64 {
		return k_r_is_n.AtVec(i) * v
	}, &temp1)

	// Python: np.dot(f_xot_is_is_n_pls, temp2)
	var result mat.Dense
	result.Mul(f_xot_is_is_n_pls, &temp1)

	return &result
}

/*
Args:

	a_s_js: 境界 j の面積, m2, [j, 1]
	beta_is_n: ステップ n からステップ n+1 における室 i の放射暖冷房設備の対流成分比率, -, [i, 1]
	f_wsb_js_is_n_pls: ステップ n+1 における係数 f_WSB, K/W, [j, 1]
	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	p_is_js: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [i, j]

Returns:

	ステップ n における係数 f_BRL, -, [i, i]

Notes:

	式(2.10)
*/
func get_f_brl_is_is_n(
	a_s_js mat.Vector,
	beta_is_n mat.Vector,
	f_wsb_js_is_n_pls mat.Matrix,
	h_s_c_js mat.Vector,
	p_is_js mat.Matrix,
) *mat.Dense {
	// f_wsb_js_is_n_pls*h_s_c_js*a_s_js
	var f_wsb_js_is_n_pls_h_s_c_js_a_s_js mat.Dense
	f_wsb_js_is_n_pls_h_s_c_js_a_s_js.Apply(func(i, j int, v float64) float64 {
		return v * h_s_c_js.AtVec(j) * a_s_js.AtVec(j)
	}, f_wsb_js_is_n_pls)

	// np.dot(p_is_js, f_wsb_js_is_n_pls*h_s_c_js*a_s_js)
	var p_is_js_f_wsb_js_is_n_pls_h_s_c_js_a_s_js mat.Dense
	p_is_js_f_wsb_js_is_n_pls_h_s_c_js_a_s_js.Mul(p_is_js, &f_wsb_js_is_n_pls_h_s_c_js_a_s_js)

	// v_diag(beta_is_n)
	beta_is_n_diag := mat.NewDiagDense(beta_is_n.Len(), nil)
	for i := 0; i < beta_is_n.Len(); i++ {
		beta_is_n_diag.SetDiag(i, beta_is_n.AtVec(i))
	}

	// np.dot(p_is_js, f_wsb_js_is_n_pls*h_s_c_js*a_s_js) + v_diag(beta_is_n)
	var result mat.Dense
	result.Add(&p_is_js_f_wsb_js_is_n_pls_h_s_c_js_a_s_js, beta_is_n_diag)

	return &result
}

/*
Args:

	f_flb_js_is_n_pls: ステップ n+1 における係数 f_FLB, K/W, [j, i]
	f_ax_js_js: 係数 f_AX, -, [j, j]

Returns:

	ステップ n+1 における係数 f_WSB, K/W, [j, i]

Notes:

	式(2.11)
*/
func get_f_wsb_js_is_n_pls(
	f_flb_js_is_n_pls mat.Matrix,
	f_ax_js_js *mat.LU,
) *mat.Dense {

	var result mat.Dense
	f_ax_js_js.SolveTo(&result, false, f_flb_js_is_n_pls)

	return &result
}

/*
Args:

	a_s_js: 境界 j の面積, m2, [j, 1]
	beta_is_n: ステップ n からステップ n+1 における室 i の放射暖冷房設備の対流成分比率, -, [i, 1]
	f_flr_js_is_n: ステップ n からステップ n+1 における室 i の放射暖冷房設備の放熱量の放射成分に対する境界 j の室内側表面の吸収比率, -, [j, i]
	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	h_s_r_js: 境界 j の室内側放射熱伝達率, W/(m2 K), [j, 1]
	k_ei_js_js: 境界 j の裏面温度に境界　j* の等価温度が与える影響, -, [j*, j]
	phi_a0_js: 境界 j の吸熱応答係数の初項, m2 K/W, [j]
	phi_t0_js: 境界 |j| の貫流応答係数の初項, -, [j]

Returns:

	ステップ n+1 における係数 f_FLB, K/W, [j, i]

Notes:

	式(2.12)
*/
func get_f_flb_js_is_n_pls(
	a_s_js mat.Vector,
	beta_is_n mat.Vector,
	f_flr_js_is_n mat.Matrix,
	h_s_c_js mat.Vector,
	h_s_r_js mat.Vector,
	k_ei_js_js mat.Matrix,
	phi_a0_js mat.Vector,
	phi_t0_js mat.Vector,
) *mat.Dense {
	var term1, term0 mat.Dense

	// f_flr_js_is_n*(1.0-beta_is_n.T)
	term1.Apply(func(i, j int, v float64) float64 {
		return 1.0 - v
	}, beta_is_n.T())
	term0.Apply(func(i, j int, v float64) float64 {
		return v * term1.At(0, j)
	}, f_flr_js_is_n)

	var term2 mat.VecDense
	// phi_a0_js/a_s_js
	term2.DivElemVec(phi_a0_js, a_s_js)

	// f_flr_js_is_n*(1.0-beta_is_n.T) * phi_a0_js/a_s_js
	term0.Apply(func(i, j int, v float64) float64 {
		return v * term2.AtVec(i)
	}, &term0)

	var term3 mat.Dense
	// np.dot(k_ei_js_js, f_flr_js_is_n*(1.0-beta_is_n.T))
	term3.Mul(k_ei_js_js, &term0)

	var term4 mat.VecDense
	// phi_t0_js/(h_s_c_js+h_s_r_js)/a_s_js
	term4.AddVec(h_s_c_js, h_s_r_js)
	term4.DivElemVec(phi_t0_js, &term4)
	term4.DivElemVec(&term4, a_s_js)

	// np.dot(k_ei_js_js, f_flr_js_is_n*(1.0-beta_is_n.T))*phi_t0_js/(h_s_c_js+h_s_r_js)/a_s_js
	term3.Apply(func(i, j int, v float64) float64 {
		return v * term4.AtVec(i)
	}, &term3)

	term0.Add(&term0, &term3)

	return &term0
}

/*
Args:

	beta_c_is: 室 i の放射冷房設備の対流成分比率, -, [i, 1]
	beta_h_is: 室 i の放射暖房設備の対流成分比率, -, [i, 1]
	operation_mode_is_n: ステップnにおける室iの運転モード, [i, 1]

Returns:

	ステップ n からステップ n+1 における室 i の放射暖冷房設備の対流成分比率, -, [i, 1]

Notes:

	式(2.13)
*/
func get_beta_is_n(
	beta_c_is []float64,
	beta_h_is []float64,
	operation_mode_is_n []OperationMode,
) *mat.VecDense {
	n := len(beta_c_is)
	result := mat.NewVecDense(n, nil)

	for i := 0; i < n; i++ {
		var beta float64
		if operation_mode_is_n[i] == HEATING {
			beta += beta_c_is[i]
		} else if operation_mode_is_n[i] == COOLING {
			beta += beta_h_is[i]
		} else {
			//PASS
		}
		result.SetVec(i, beta)
	}

	return result
}

/*
Args:

	f_flr_c_js_is: 室 i の放射冷房設備の放熱量の放射成分に対する境界 j の室内側表面の吸収比率, -, [j, i]
	f_flr_h_js_is: 室 i の放射暖房設備の放熱量の放射成分に対する境界 j の室内側表面の吸収比率, -, [j, i]
	is_cooling_is_n: 「ステップ n から n+1 における室 i の運転が冷房運転時の場合」かの有無, -, [i, 1]
	is_heating_is_n: 「ステップ n から n+1 における室 i の運転が暖房運転時の場合」かの有無, -, [i, 1]

Returns:

	ステップ n からステップ n+1 における室 i の放射暖冷房設備の放熱量の放射成分に対する境界 j の室内側表面の吸収比率, -, [j, i]

Notes:

	式(2.14)
*/
func get_f_flr_js_is_n(
	f_flr_c_js_is mat.Matrix,
	f_flr_h_js_is mat.Matrix,
	operation_mode_is_n []OperationMode,
) *mat.Dense {
	r, c := f_flr_c_js_is.Dims()
	result := make([]float64, r*c)
	off := 0
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			var f_flr float64
			if operation_mode_is_n[j] == HEATING {
				f_flr = f_flr_c_js_is.At(i, j)
			} else if operation_mode_is_n[j] == COOLING {
				f_flr = f_flr_h_js_is.At(i, j)
			} else {
				//PASS
			}
			result[off] = f_flr
			off++
		}
	}
	return mat.NewDense(r, c, result)
}

/*
Args:

	f_xc_is_n_pls: ステップ n+1 における係数 f_XC, degree C, [i, 1]
	f_brc_non_nv_is_n_pls: ステップn+1における自然風非利用時の係数f_BRC,OT, W, [i,1]
	f_brc_nv_is_n_pls: ステップn+1における自然風利用時の係数f_BRC,OT, W, [i,1]
	f_brm_non_nv_is_is_n_pls: ステップn+1における自然風非利用時の係数f_BRM,OT, W, [i,i]
	f_brm_nv_is_is_n_pls: ステップn+1における自然風利用時の係数f_BRM,OT, W, [i,i]

Returns:

	ステップn+1における自然風非利用時の係数f_BRC,OT, W, [i, 1]
	ステップn+1における自然風利用時の係数f_BRC,OT, W, [i, 1]

Notes:

	式(2.17)
*/
func get_f_brc_ot_is_n_pls(
	f_xc_is_n_pls []float64,
	f_brc_non_nv_is_n_pls mat.Vector,
	f_brc_nv_is_n_pls mat.Vector,
	f_brm_non_nv_is_is_n_pls mat.Matrix,
	f_brm_nv_is_is_n_pls mat.Matrix,
) (mat.Vector, mat.Vector) {
	var temp1, temp2 mat.VecDense

	_f_xc_is_n_pls := mat.NewVecDense(len(f_xc_is_n_pls), f_xc_is_n_pls)

	// f_brc_non_nv_is_n_pls + np.dot(f_brm_non_nv_is_is_n_pls, f_xc_is_n_pls)
	temp1.MulVec(f_brm_non_nv_is_is_n_pls, _f_xc_is_n_pls)
	temp1.AddVec(f_brc_non_nv_is_n_pls, &temp1)

	// f_brc_nv_is_n_pls + np.dot(f_brm_nv_is_is_n_pls, f_xc_is_n_pls)
	temp2.MulVec(f_brm_nv_is_is_n_pls, _f_xc_is_n_pls)
	temp2.AddVec(f_brc_nv_is_n_pls, &temp2)

	return &temp1, &temp2
}

/*
Args:

	f_xot_is_is_n_pls: ステップ n+1 における係数 f_XOT, -, [i, i]
	f_brm_non_nv_is_is_n_pls: ステップ n+1 における自然風非利用時の係数 f_BRM, W/K, [i, i]
	f_brm_nv_is_is_n_pls: ステップ n+1 における自然風利用時の係数 f_BRM, W/K, [i, i]

Returns:

	ステップn+1における自然風非利用時の係数f_BRM,OT, W/K, [i, 1]
	ステップn+1における自然風利用時の係数f_BRM,OT, W/K, [i, 1]

Notes:

	式(2.18)
*/
func get_f_brm_ot_is_is_n_pls(f_xot_is_is_n_pls mat.Matrix, f_brm_non_nv_is_is_n_pls mat.Matrix, f_brm_nv_is_is_n_pls mat.Matrix) (mat.Matrix, mat.Matrix) {
	var temp1, temp2 mat.Dense
	temp1.Mul(f_brm_non_nv_is_is_n_pls, f_xot_is_is_n_pls)
	temp2.Mul(f_brm_nv_is_is_n_pls, f_xot_is_is_n_pls)
	return &temp1, &temp2
}

/*
	    Args:
	        f_brc_ot_non_nv_is_n_pls: ステップ n+1 における自然風の利用なし時の係数 f_BRC,OT, W, [i, 1]
	        f_brc_ot_nv_is_n_pls: ステップ n+1 における自然風の利用時の係数 f_BRC,OT, W, [i, 1]
	        f_brm_ot_non_nv_is_is_n_pls: ステップ n+1 における自然風の利用なし時の係数 f_BRM,OT, W/K, [i, 1]
	        f_brm_ot_nv_is_is_n_pls: ステップ n+1 における自然風の利用時の係数 f_BRM,OT, W/K, [i, 1]


	    Returns:
	        ステップn+1における自然風非利用時の室iの自然作用温度, degree C, [i, 1]
	        ステップn+1における自然風利用時の室iの自然作用温度, degree C, [i, 1]

		Notes:
			式(2.16)
*/
func get_theta_r_ot_ntr_is_n_pls(
	f_brc_ot_non_nv_is_n_pls mat.Vector,
	f_brc_ot_nv_is_n_pls mat.Vector,
	f_brm_ot_non_nv_is_is_n_pls mat.Matrix,
	f_brm_ot_nv_is_is_n_pls mat.Matrix,
) (*mat.VecDense, *mat.VecDense) {

	var result1, result2 mat.VecDense

	result1.SolveVec(f_brm_ot_non_nv_is_is_n_pls, f_brc_ot_non_nv_is_n_pls)
	result2.SolveVec(f_brm_ot_nv_is_is_n_pls, f_brc_ot_nv_is_n_pls)

	return &result1, &result2
}

/*
Args:

	f_mrt_hum_is_js: 室 i の人体に対する境界 j の形態係数, -, [i, j]
	f_wsc_js_n_pls: ステップ n+1 における係数 f_WSC, degree C, [j, 1]
	f_wsv_js_n_pls: ステップ n+1 における係数 f_WSV, degree C, [j, 1]
	f_xot_is_is_n_pls: ステップ n+1 における係数 f_XOT, -, [i, i]
	k_r_is_n: ステップ n における室 i の人体表面の放射熱伝達率が総合熱伝達率に占める割合, -, [i, 1]

Returns:

	ステップ n+1 における係数 f_XC, degree C, [i, 1]

Notes:

	式(2.19)
*/
func get_f_xc_is_n_pls(
	f_mrt_hum_is_js mat.Matrix,
	f_wsc_js_n_pls mat.Vector,
	f_wsv_js_n_pls mat.Vector,
	f_xot_is_is_n_pls mat.Matrix,
	k_r_is_n mat.Vector,
) []float64 {
	// f_wsc_js_n_pls + f_wsv_js_n_pls
	var f_wscPlusWsv mat.VecDense
	f_wscPlusWsv.AddVec(f_wsc_js_n_pls, f_wsv_js_n_pls)

	// np.dot(f_mrt_hum_is_js, (f_wsc_js_n_pls + f_wsv_js_n_pls))
	var dotProduct mat.VecDense
	dotProduct.MulVec(f_mrt_hum_is_js, &f_wscPlusWsv)

	// k_r_is_n * np.dot(f_mrt_hum_is_js, (f_wsc_js_n_pls + f_wsv_js_n_pls))
	var elementwiseProduct mat.VecDense
	elementwiseProduct.MulElemVec(k_r_is_n, &dotProduct)

	// np.dot(f_xot_is_is_n_pls, k_r_is_n * np.dot(f_mrt_hum_is_js, (f_wsc_js_n_pls + f_wsv_js_n_pls)))
	var f_xc_is_n_pls mat.VecDense
	f_xc_is_n_pls.MulVec(f_xot_is_is_n_pls, &elementwiseProduct)

	slice := make([]float64, f_xc_is_n_pls.Len())
	copy(slice, f_xc_is_n_pls.RawVector().Data)

	return slice
}

/*
Args:

	f_mrt_hum_is_js: 室 i の人体に対する境界 j の形態係数, -, [i, j]
	f_wsr_js_is: 係数 f_WSR, - [j, i]
	k_c_is_n: ステップ n における室 i の人体表面の対流熱伝達率が総合熱伝達率に占める割合, -, [i, 1]
	k_r_is_n: ステップ n における室 i の人体表面の放射熱伝達率が総合熱伝達率に占める割合, -, [i, 1]

Returns:

	ステップ n+1 における係数 f_XOT, -, [i, i]

Notes:

	式(2.20)
*/
func get_f_xot_is_is_n_pls(
	f_mrt_hum_is_js mat.Matrix,
	f_wsr_js_is mat.Matrix,
	k_c_is_n mat.Vector,
	k_r_is_n mat.Vector,
) *mat.Dense {
	// ベクトルを対角行列に変換
	vDiag := func(v mat.Vector) *mat.DiagDense {
		n := v.Len()
		diag := mat.NewDiagDense(n, nil)
		for i := 0; i < n; i++ {
			diag.SetDiag(i, v.AtVec(i))
		}
		return diag
	}

	// Python: v_diag(k_c_is_n)
	temp1 := vDiag(k_c_is_n)

	// Python: np.dot(f_mrt_hum_is_js, f_wsr_js_is)
	temp2 := &mat.Dense{}
	temp2.Mul(f_mrt_hum_is_js, f_wsr_js_is)

	// Python: k_r_is_n * temp2
	temp3 := &mat.Dense{}
	temp3.Apply(func(_, j int, v float64) float64 {
		return v * k_r_is_n.AtVec(j)
	}, temp2)

	// Python: temp1 + temp3
	temp4 := &mat.Dense{}
	temp4.Add(temp1, temp3)

	// Python: np.linalg.inv(temp4)
	result := &mat.Dense{}
	err := result.Inverse(temp4)
	if err != nil {
		panic(err)
	}

	return result
}

/*
Args:

	h_hum_c_is_n: ステップ n における室 i の人体表面の対流熱伝達率, W/(m2 K)
	h_hum_r_is_n: ステップ n における室 i の人体表面の放射熱伝達率, W/(m2 K)

Returns:

	ステップ n における室 i の人体表面の対流熱伝達率が総合熱伝達率に占める割合, -, [i, 1]

Notes:

	式(2.21)
*/
func get_k_c_is_n(h_hum_c_is_n mat.Vector, h_hum_r_is_n mat.Vector) *mat.VecDense {
	// h_hum_c_is_n + h_hum_r_is_n
	var sum mat.VecDense
	sum.AddVec(h_hum_c_is_n, h_hum_r_is_n)

	// h_hum_c_is_n / (h_hum_c_is_n + h_hum_r_is_n)
	var result mat.VecDense
	result.DivElemVec(h_hum_c_is_n, &sum)

	return &result
}

/*
Args:

	h_hum_c_is_n: ステップ n における室 i の人体表面の対流熱伝達率, W/(m2 K)
	h_hum_r_is_n: ステップ n における室 i の人体表面の放射熱伝達率, W/(m2 K)

Returns:

	ステップ n における室 i の人体表面の放射熱伝達率が総合熱伝達率に占める割合, -, [i, 1]

Notes:

	式(2.22)
*/
func get_k_r_is_n(h_hum_c_is_n mat.Vector, h_hum_r_is_n mat.Vector) *mat.VecDense {
	// h_hum_c_is_n + h_hum_r_is_n
	var sum mat.VecDense
	sum.AddVec(h_hum_c_is_n, h_hum_r_is_n)

	// h_hum_r_is_n / (h_hum_c_is_n + h_hum_r_is_n)
	var result mat.VecDense
	result.DivElemVec(h_hum_r_is_n, &sum)

	return &result
}

/*
Args:

	a_s_js: 境界 j の面積, m2, [j, 1]
	v_rm_is: 室 i の容積, m3, [i, 1]
	c_sh_frt_is: 室 i の備品等の熱容量, J/K, [i, 1]
	delta_t: 1ステップの時間間隔, s
	f_wsr_js_is: 係数 f_WSR, - [j, i]
	g_sh_frt_is: 室 i の備品等と空気間の熱コンダクタンス, W/K, [i, 1]
	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	p_is_js: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [i, j]
	p_js_is: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [j, i]
	v_vent_int_is_is_n: ステップ n から ステップ n+1 における室 i* から室 i への室間の空気移動量（流出換気量を含む）, m3/s
	v_vent_out_is_n: ステップ n からステップ n+1 における室 i の換気・すきま風・自然風の利用による外気の流入量, m3/s

Returns:

	ステップ n+1 における係数 f_BRM, W/K, [i, i]

Notes:

	式(2.23)
*/
func get_f_brm_is_is_n_pls(
	a_s_js mat.Vector,
	v_rm_is []float64,
	c_sh_frt_is mat.Vector,
	delta_t float64,
	f_wsr_js_is mat.Matrix,
	g_sh_frt_is mat.Vector,
	h_s_c_js mat.Vector,
	p_is_js mat.Matrix,
	p_js_is mat.Matrix,
	v_vent_int_is_is_n mat.Matrix,
	v_vent_out_is_n []float64,
	v_vent_ntr_set_is []float64,
) (*mat.Dense, *mat.Dense) {
	vDiag := func(v *mat.VecDense) *mat.DiagDense {
		n := v.Len()
		diag := mat.NewDiagDense(n, v.RawVector().Data)
		return diag
	}

	// Python: v_rm_is * rho_a * c_a / delta_t
	temp1 := make([]float64, len(v_rm_is))
	floats.ScaleTo(temp1, rho_a*c_a/delta_t, v_rm_is)

	// Python: (p_js_is - f_wsr_js_is) * a_s_js * h_s_c_js
	var temp2 mat.Dense
	temp2.Apply(func(j, i int, v float64) float64 {
		return (p_js_is.At(j, i) - v) * a_s_js.AtVec(j) * h_s_c_js.AtVec(j)
	}, f_wsr_js_is)

	// Python: np.dot(p_is_js, temp2)
	var temp3 mat.Dense
	temp3.Mul(p_is_js, &temp2)

	// Python: c_sh_frt_is * g_sh_frt_is / (c_sh_frt_is + g_sh_frt_is * delta_t)
	var temp4 mat.VecDense
	temp4.MulElemVec(c_sh_frt_is, g_sh_frt_is)
	var temp5 mat.VecDense
	temp5.ScaleVec(delta_t, g_sh_frt_is)
	temp5.AddVec(c_sh_frt_is, &temp5)
	temp4.DivElemVec(&temp4, &temp5)

	// Python: v_diag(temp4)
	temp6 := vDiag(&temp4)

	// Python: c_a * rho_a * (v_diag(v_vent_out_is_n) - v_vent_int_is_is_n)
	temp7 := mat.NewDiagDense(len(v_vent_out_is_n), v_vent_out_is_n)
	var temp8 mat.Dense
	temp8.Sub(temp7, v_vent_int_is_is_n)
	temp8.Scale(c_a*rho_a, &temp8)

	// Python: temp1 + temp3 + temp6 + temp8
	var result1 mat.Dense
	result1.Add(mat.NewDiagDense(len(temp1), temp1), &temp3)
	result1.Add(&result1, temp6)
	result1.Add(&result1, &temp8)

	var result2, temp9 mat.Dense
	temp9.Scale(c_a*rho_a, mat.NewDiagDense(len(v_vent_ntr_set_is), v_vent_ntr_set_is))
	result2.Add(&result1, &temp9)

	return &result1, &result2
}

/*
Args:

	a_s_js: 境界 j の面積, m2, [j, 1]
	v_rm_is: 室容量, m3, [i, 1]
	c_sh_frt_is: 室 i の備品等の熱容量, J/K, [i, 1]
	delta_t: 1ステップの時間間隔, s
	f_wsc_js_n_pls: ステップ n+1 における係数 f_WSC, degree C, [j, 1]
	f_wsv_js_n_pls: ステップ n+1 における係数 f_WSV, degree C, [j, 1]
	g_sh_frt_is: 室 i の備品等と空気間の熱コンダクタンス, W/K, [i, 1]
	h_s_c_js: 境界 j の室内側対流熱伝達率, W/(m2 K), [j, 1]
	p_is_js: 室 i と境界 j の接続に関する係数（境界 j が室 i に接している場合は 1 とし、それ以外の場合は 0 とする。）, -, [i, j]
	q_gen_is_n: ステップ n からステップ n+1 における室 i の人体発熱を除く内部発熱, W, [i, 1]
	q_hum_is_n: ステップ n からステップ n+1 における室 i の人体発熱, W, [i, 1]
	q_sol_frt_is_n: ステップ n からステップ n+1 における室 i に設置された備品等による透過日射吸収熱量時間平均値, W, [i, 1]
	theta_frt_is_n: ステップ n における室 i の備品等の温度, degree C, [i, 1]
	theta_o_n_pls: ステップ n+1 における外気温度, ℃
	theta_r_is_n: ステップ n における室 i の温度, ℃
	v_vent_out_non_nv_is_n: ステップ n からステップ n+1 における室 i の換気・すきま風・自然風の利用による外気の流入量, m3/s
	v_vent_ntr_is_n: ステップnからステップn+1における室iの自然風の利用による外気の流入量, m3/s

Returns:

	ステップ n+1 における係数 f_BRC,OT, W, [i, 1]

Notes:

	式(2.24)
*/
func get_f_brc_is_n_pls(
	a_s_js mat.Vector,
	v_rm_is []float64,
	c_sh_frt_is mat.Vector,
	delta_t float64,
	f_wsc_js_n_pls mat.Vector,
	f_wsv_js_n_pls mat.Vector,
	g_sh_frt_is mat.Vector,
	h_s_c_js mat.Vector,
	p_is_js mat.Matrix,
	q_gen_is_n []float64,
	q_hum_is_n []float64,
	q_sol_frt_is_n mat.Vector,
	theta_frt_is_n []float64,
	theta_o_n_pls float64,
	theta_r_is_n []float64,
	v_vent_out_non_nv_is_n []float64,
	v_vent_ntr_is_n []float64,
) (*mat.VecDense, *mat.VecDense) {
	// v_rm_is * c_a * rho_a / delta_t * theta_r_is_n
	_result1_slice := make([]float64, len(v_rm_is))
	floats.MulTo(_result1_slice, v_rm_is, theta_r_is_n)
	floats.Scale(c_a*rho_a/delta_t, _result1_slice)
	result1 := mat.NewVecDense(len(_result1_slice), _result1_slice)

	// (f_wsc_js_n_pls + f_wsv_js_n_pls)
	var tmpVec1 mat.VecDense
	tmpVec1.AddVec(f_wsc_js_n_pls, f_wsv_js_n_pls)

	// h_s_c_js * a_s_js * (f_wsc_js_n_pls + f_wsv_js_n_pls)
	var tmpVec2 mat.VecDense
	tmpVec2.MulElemVec(h_s_c_js, a_s_js)
	tmpVec2.MulElemVec(&tmpVec2, &tmpVec1)

	// + np.dot(p_is_js, h_s_c_js * a_s_js * (f_wsc_js_n_pls + f_wsv_js_n_pls))
	var tmpMat1 mat.VecDense
	tmpMat1.MulVec(p_is_js, &tmpVec2)
	result1.AddVec(result1, &tmpMat1)

	// + c_a * rho_a * v_vent_out_non_nv_is_n * theta_o_n_pls
	tmpVec3 := make([]float64, len(v_vent_out_non_nv_is_n))
	floats.ScaleTo(tmpVec3, c_a*rho_a*theta_o_n_pls, v_vent_out_non_nv_is_n)
	result1.AddVec(result1, mat.NewVecDense(len(tmpVec3), tmpVec3))

	// + q_gen_is_n + q_hum_is_n
	result1.AddVec(result1, mat.NewVecDense(len(q_gen_is_n), q_gen_is_n))
	result1.AddVec(result1, mat.NewVecDense(len(q_hum_is_n), q_hum_is_n))

	// c_sh_frt_is * theta_frt_is_n
	var tmpVec4 mat.VecDense
	tmpVec4.MulElemVec(c_sh_frt_is, mat.NewVecDense(len(theta_frt_is_n), theta_frt_is_n))

	// q_sol_frt_is_n * delta_t
	var tmpVec5 mat.VecDense
	tmpVec5.ScaleVec(delta_t, q_sol_frt_is_n)

	// c_sh_frt_is * theta_frt_is_n + q_sol_frt_is_n * delta_t
	tmpVec4.AddVec(&tmpVec4, &tmpVec5)

	// c_sh_frt_is + delta_t * g_sh_frt_is
	var tmpVec6 mat.VecDense
	tmpVec6.AddScaledVec(c_sh_frt_is, delta_t, g_sh_frt_is)

	// g_sh_frt_is * (c_sh_frt_is * theta_frt_is_n + q_sol_frt_is_n * delta_t) / (c_sh_frt_is + delta_t * g_sh_frt_is)
	tmpVec4.DivElemVec(&tmpVec4, &tmpVec6)
	tmpVec4.MulElemVec(g_sh_frt_is, &tmpVec4)

	// + g_sh_frt_is * (c_sh_frt_is * theta_frt_is_n + q_sol_frt_is_n * delta_t) / (c_sh_frt_is + delta_t * g_sh_frt_is)
	result1.AddVec(result1, &tmpVec4)

	//  result1 + c_a * rho_a * v_vent_ntr_is_n * theta_o_n_pls
	var result2 mat.VecDense
	result2.AddScaledVec(result1, c_a*rho_a*theta_o_n_pls, mat.NewVecDense(len(v_vent_ntr_is_n), v_vent_ntr_is_n))

	return result1, &result2
}

/*
Args:

	v_leak_is_n: ステップ n からステップ n+1 における室 i のすきま風量, m3/s, [i, 1]
	v_vent_mec_is_n: ステップ n からステップ n+1 における室 i の機械換気量（全般換気量と局所換気量の合計値）, m3/s, [i, 1]

Returns:

	ステップ n からステップ n+1 における室 i の換気・すきま風・自然風の利用による外気の流入量, m3/s

Notes:

	式(2.25)
*/
func get_v_vent_out_non_ntr_is_n(
	v_leak_is_n []float64,
	v_vent_mec_is_n []float64,
) []float64 {
	result := make([]float64, len(v_leak_is_n))
	floats.AddTo(result, v_leak_is_n, v_vent_mec_is_n)
	return result
}

/*
Args:

	operation_mode_is_n: ステップ n からステップ n+1 における室 i の運転モード, [i, 1]
	v_vent_ntr_set_is: 室 i の自然風利用時の換気量, m3/s

Returns:

	ステップ n からステップ n+1 における室 i の自然風利用による換気量, m3/s, [i, 1]

Notes:

	式(2.26)
*/
func get_v_vent_ntr_is_n(operation_mode_is_n []OperationMode, v_vent_ntr_set_is mat.Vector) *mat.VecDense {
	result := mat.NewVecDense(len(operation_mode_is_n), nil)
	for i := 0; i < len(operation_mode_is_n); i++ {
		if operation_mode_is_n[i] == STOP_OPEN {
			result.SetVec(i, v_vent_ntr_set_is.AtVec(i))
		} else {
			result.SetVec(i, 0.0)
		}
	}
	return result
}

/*
Args:

	f_cvl_js_n_pls: ステップ n+1 における係数 f_CVL, degree C, [j, 1]
	f_ax_js_js: 係数 f_AX, -, [j, j]

Returns:

	ステップ n+1 の係数 f_WSV, degree C, [j, 1]

Notes:

	式(2.27)
*/
func get_f_wsv_js_n_pls(
	f_cvl_js_n_pls []float64,
	f_ax_js_js *mat.LU,
) *mat.VecDense {
	var result mat.VecDense
	f_ax_js_js.SolveVecTo(&result, false, mat.NewVecDense(len(f_cvl_js_n_pls), f_cvl_js_n_pls))

	return &result
}

/*
Args:

	theta_dsh_s_a_js_ms_n_pls: ステップ n+1 における境界 j の項別公比法の指数項 m の吸熱応答の項別成分, degree C, [j, m]
	theta_dsh_s_t_js_ms_n_pls: ステップ n+1 における境界 j の項別公比法の指数項 m の貫流応答の項別成分, degree C, [j, m]

Returns:

	ステップ n+1 における係数 f_CVL, degree C, [j, 1]

Notes:

	式(2.28)
*/
func get_f_cvl_js_n_pls(
	theta_dsh_s_a_js_ms_n_pls [][]float64,
	theta_dsh_s_t_js_ms_n_pls [][]float64,
) []float64 {
	// Pythonコード: np.sum(theta_dsh_s_t_js_ms_n_pls + theta_dsh_s_a_js_ms_n_pls, axis=1, keepdims=True)
	l := len(theta_dsh_s_t_js_ms_n_pls)
	result := make([]float64, l)
	for i := 0; i < l; i++ {
		rowSum := 0.0
		c := len(theta_dsh_s_t_js_ms_n_pls[i])
		for j := 0; j < c; j++ {
			rowSum += theta_dsh_s_t_js_ms_n_pls[i][j] + theta_dsh_s_a_js_ms_n_pls[i][j]
		}
		result[i] = rowSum
	}

	return result
}

/*
Args:

	phi_a1_js_ms: 境界 j の項別公比法の指数項 m の吸熱応答係数, m2 K/W, [j, m]
	q_s_js_n: ステップ n における境界 j の表面熱流（壁体吸熱を正とする）, W/m2, [j, 1]
	r_js_ms: 境界 j の項別公比法の指数項 m の公比, -, [j, m]
	theta_dsh_srf_a_js_ms_n: ステップ n における境界 j の項別公比法の指数項 m の吸熱応答の項別成分, degree C, [j, m]

Returns:

	ステップ n+1 における境界 j の項別公比法の指数項 m の吸熱応答の項別成分, degree C, [j, m]

Notes:

	式(2.29)
*/
func get_theta_dsh_s_a_js_ms_n_pls(
	phi_a1_js_ms [][]float64,
	q_s_js_n []float64,
	r_js_ms [][]float64,
	theta_dsh_srf_a_js_ms_n [][]float64,
) [][]float64 {
	result := make([][]float64, len(phi_a1_js_ms))
	for i := 0; i < len(phi_a1_js_ms); i++ {
		result[i] = make([]float64, len(phi_a1_js_ms[i]))
		for j := 0; j < len(phi_a1_js_ms[i]); j++ {
			result[i][j] = phi_a1_js_ms[i][j]*q_s_js_n[i] + r_js_ms[i][j]*theta_dsh_srf_a_js_ms_n[i][j]
		}
	}
	return result
}

/*
Args:

	phi_t1_js_ms: 境界 j の項別公比法の指数項 m の貫流応答係数, -, [j, m]
	r_js_ms: 境界 j の項別公比法の指数項 m の公比, -, [j, m]
	theta_dsh_srf_t_js_ms_n: ステップ n における境界 j の項別公比法の指数項 m の貫流応答の項別成分, degree C, [j, m]
	theta_rear_js_n: ステップ n における境界 j の裏面温度, degree C, [j, 1]

Returns:

	ステップ n+1 における境界 j の項別公比法の指数項 m の貫流応答の項別成分, degree C, [j, m]

Notes:

	式(2.30)
*/
func get_theta_dsh_s_t_js_ms_n_pls(
	phi_t1_js_ms [][]float64,
	r_js_ms [][]float64,
	theta_dsh_srf_t_js_ms_n [][]float64,
	theta_rear_js_n []float64,
) [][]float64 {
	result := make([][]float64, len(phi_t1_js_ms))
	for j := 0; j < len(phi_t1_js_ms); j++ {
		result[j] = make([]float64, len(phi_t1_js_ms[j]))
		for m := 0; m < len(phi_t1_js_ms[j]); m++ {
			result[j][m] = phi_t1_js_ms[j][m]*theta_rear_js_n[j] + r_js_ms[j][m]*theta_dsh_srf_t_js_ms_n[j][m]
		}
	}

	return result
}

/*
Args:

	n_hum_is_n: ステップ n からステップ n+1 における室 i の在室人数, -, [i, 1]
	q_hum_psn_is_n: ステップ n からステップ n+1 における室 i の1人あたりの人体発熱, W, [i, 1]

Returns:

	ステップ n からステップ n+1 における室 i の人体発熱, W, [i, 1]

Notes:

	式(2.31)
*/
func get_q_hum_is_n(n_hum_is_n []float64, q_hum_psn_is_n []float64) []float64 {
	result := make([]float64, len(q_hum_psn_is_n))
	floats.MulTo(result, n_hum_is_n, q_hum_psn_is_n)
	return result
}

/*
Args:

	k_ei_js_js: 境界 j の裏面温度に境界　j* の等価温度が与える影響, -, [j*, j]
	theta_ei_js_n: ステップ n における境界 j の等価温度, degree C, [j, 1]
	k_eo_js: 温度差係数, -, [j, 1]
	theta_o_eqv_js_n: ステップ n の境界 j における相当外気温度, ℃, [j, n]

Returns:

	ステップ n における境界 j の裏面温度, degree C, [j, 1]

Notes:

	式(2.32)
*/
func get_theta_s_rear_js_n(
	k_s_er_js_js mat.Matrix,
	theta_er_js_n []float64,
	k_s_eo_js mat.Vector,
	theta_eo_js_n mat.Vector,
	k_s_r_js_is mat.Matrix,
	theta_r_is_n []float64,
) []float64 {

	var result1, result2, result3 mat.VecDense

	// np.dot(k_s_er_js_js, theta_er_js_n)
	result1.MulVec(k_s_er_js_js, mat.NewVecDense(len(theta_er_js_n), theta_er_js_n))

	// k_s_eo_js*theta_eo_js_n
	result2.MulElemVec(k_s_eo_js, theta_eo_js_n)

	// np.dot(k_s_r_js_is, theta_r_is_n)
	result3.MulVec(k_s_r_js_is, mat.NewVecDense(len(theta_r_is_n), theta_r_is_n))

	//np.dot(k_s_er_js_js, theta_er_js_n) + k_s_eo_js * theta_eo_js_n + np.dot(k_s_r_js_is, theta_r_is_n)
	var finalResult mat.VecDense
	finalResult.AddVec(&result1, &result2)
	finalResult.AddVec(&finalResult, &result3)

	result := finalResult.RawVector()
	if result.Inc != 1 {
		panic("result.Inc != 1")
	}

	return result.Data
}
