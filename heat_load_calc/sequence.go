package heat_load_calc

import (
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

// 事前計算結果
// See: https://hc-energy.readthedocs.io/ja/latest/contents/03_09_eval_sequence.html#id12
type PreCalcParameters struct {
	// ステップnの室iにおける機械換気量（全般換気量+局所換気量）, m3/s, [i, 8760*4]
	v_vent_mec_is_ns *ScheduleData

	// ステップ n における室 i に設置された備品等による透過日射吸収熱量, W, [i, n+1]
	// NOTE: 計算仕様では扱わない
	q_sol_frt_is_ns mat.Matrix

	// ステップnの境界jにおける透過日射熱取得量のうち表面に吸収される日射量, W/m2, [j, 8760*4]
	q_s_sol_js_ns mat.Matrix

	// 係数 f_AX (LU分解済み) , [i, i] 式(4.5)
	f_ax_js_js *mat.LU

	// 係数 f_AX^1, [i, i] 式(4.5)
	f_ax_js_js_inv *mat.Dense

	// 室iの在室者に対する境界j*の形態係数
	f_mrt_hum_is_js mat.Matrix

	// 平均放射温度計算時の境界 j* の表面温度が境界 j に与える重み, [j, j]
	f_mrt_is_js *mat.Dense

	// WSR, WSB の計算 式(4.1)
	f_wsr_js_is mat.Matrix

	// WSC, W, [j, n] 式(4.2)
	f_wsc_js_ns *mat.Dense

	// ステップ n における室 i の在室者表面における放射熱伝達率の総合熱伝達率に対する比, -, [i, 1] 式(2.21)
	// NOTE: 計算仕様における実装と異なり、決め打ちの値が入る
	k_r_is_n *mat.VecDense

	// ステップnにおける室iの在室者表面における対流熱伝達率の総合熱伝達率に対する比, -, [i, 1] 式(2.22)
	// NOTE: 計算仕様における実装と異なり、決め打ちの値が入る
	k_c_is_n *mat.VecDense

	// ステップn+1における室iの係数 XOT, [i, i] (逆行列,密行列) 式(2.20)
	f_xot_is_is_n_pls *mat.Dense
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
	get_f_l_cl          func(*mat.VecDense, *mat.VecDense, []float64) (*mat.VecDense, *mat.Dense)
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

// ----------------------------------------------------------------------------------
// 3 繰り返し計算（建物全般）
// ----------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
// 3.1 湿度と潜熱処理量
// ----------------------------------------------------------------------------------

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
	// ----------- ここから人体発熱・人体発湿 -----------

	// ステップnからステップn+1における室iの1人あたりの人体発熱, W, [i, 1]
	q_hum_psn_is_n := get_q_hum_psn_is_n(c_n.theta_r_is_n)

	// 式(2.31)
	// ステップ n からステップ n+1 における室 i の人体発熱, W, [i, 1]
	q_hum_is_n := get_q_hum_is_n(
		self.scd.n_hum_is_ns.Get(nn),
		q_hum_psn_is_n,
	)

	// ステップnの室iにおける1人あたりの人体発湿, kg/s, [i, 1]
	x_hum_psn_is_n := get_x_hum_psn_is_n(q_hum_psn_is_n)

	// 式(1.7) 人体発熱計算
	// ステップnの室iにおける人体発湿, kg/s, [i, 1]
	x_hum_is_n := get_x_hum_is_n(self.scd.n_hum_is_ns.Get(nn), x_hum_psn_is_n)

	// ----------- ここまで人体発熱・人体発湿 -----------

	// 式(2.32) 裏面温度計算
	// ステップ n の境界 j における裏面温度, degree C, [j, 1]
	theta_rear_js_n := get_theta_s_rear_js_n(
		self.bs.k_ei_js_js,
		c_n.theta_ei_js_n,
		self.bs.k_eo_js,
		self.bs.theta_o_eqv_js_ns.ColView(nn_plus),
		self.bs.k_s_r_js,
		c_n.theta_r_is_n,
	)

	// ステップnの室iにおけるすきま風量, m3/s, [i, 1]
	// 式(建物全般のパラメータ:2)
	v_leak_is_n := self.building.get_v_leak_is_n(
		c_n.theta_r_is_n,
		self.weather.theta_o_ns_plus[nn_plus],
		self.rms.v_rm_is,
	)

	// 式(2.30)
	// ステップ n+1 の境界 j における項別公比法の指数項 m の貫流応答の項別成分, degree C, [j, m] (m=12), eq.(29)
	theta_dsh_s_t_js_ms_n_pls := get_theta_dsh_s_t_js_ms_n_pls(
		self.bs.phi_t1_js_ms,
		self.bs.r_js_ms,
		c_n.theta_dsh_srf_t_js_ms_n,
		theta_rear_js_n,
	)

	// 式(2.29)
	// ステップ n+1 の境界 j における項別公比法の指数項 m の吸熱応答の項別成分, degree C, [j, m]
	theta_dsh_s_a_js_ms_n_pls := get_theta_dsh_s_a_js_ms_n_pls(
		self.bs.phi_a1_js_ms,
		c_n.q_s_js_n,
		self.bs.r_js_ms,
		c_n.theta_dsh_srf_a_js_ms_n,
	)

	// 式(2.28)
	// ステップ n+1 の境界 j における係数f_CVL, degree C, [j, 1]
	f_cvl_js_n_pls := get_f_cvl_js_n_pls(
		theta_dsh_s_a_js_ms_n_pls,
		theta_dsh_s_t_js_ms_n_pls,
	)

	// 式(2.27)
	// ステップ n+1 の境界 j における係数 f_WSV, degree C, [j, 1]
	f_wsv_js_n_pls := get_f_wsv_js_n_pls(
		f_cvl_js_n_pls,
		ss.f_ax_js_js_inv,
		ss.f_ax_js_js,
	)

	// 式(2.25)
	// ステップnからステップn+1における室iの換気・隙間風による外気の流入量, m3/s, [i, 1]
	v_vent_out_non_nv_is_n := get_v_vent_out_non_ntr_is_n(
		v_leak_is_n,
		ss.v_vent_mec_is_ns.Get(nn),
	)

	// 式(2.24)
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

	// 式(2.23)
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

	// 式(2.19)
	// ステップn+1における室iの係数 XC, [i, 1]
	f_xc_is_n_pls := get_f_xc_is_n_pls(
		ss.f_mrt_hum_is_js,
		ss.f_wsc_js_ns.ColView(nn_plus+1),
		f_wsv_js_n_pls,
		ss.f_xot_is_is_n_pls, //NOTE: 計算仕様では毎時計算になっている
		ss.k_r_is_n,
	)

	// ステップn+1における自然風の非利用時・利用時の係数f_BRM,OT, W/K, [i, i] 式(2.18)
	f_brm_ot_non_nv_is_is_n_pls, f_brm_ot_nv_is_is_n_pls := get_f_brm_ot_is_is_n_pls(
		ss.f_xot_is_is_n_pls,
		f_brm_non_nv_is_is_n_pls,
		f_brm_nv_is_is_n_pls,
	)

	// ステップ n における自然風の非利用時・利用時の係数 f_BRC,OT, W, [i, 1] 式(2.17)
	f_brc_ot_non_nv_is_n_pls, f_brc_ot_nv_is_n_pls := get_f_brc_ot_is_n_pls(
		f_xc_is_n_pls,
		f_brc_non_nv_is_n_pls,
		f_brc_nv_is_n_pls,
		f_brm_non_nv_is_is_n_pls,
		f_brm_nv_is_is_n_pls,
	)

	// 式(1.6)
	// ステップnにおける室iの自然風の非利用時・利用時の潜熱バランスに関する係数f_h_cst, kg / s, [i, 1]
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

	// 式(1.5)
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

	// 式(2.16)
	// ステップn+1における自然風非利用時の自然作用温度, degree C, [i, 1]
	// ステップn+1における自然風利用時の自然作用温度, degree C, [i, 1]
	theta_r_ot_ntr_non_nv_is_n_pls, theta_r_ot_ntr_nv_is_n_pls := get_theta_r_ot_ntr_is_n_pls(
		f_brc_ot_non_nv_is_n_pls,
		f_brc_ot_nv_is_n_pls,
		f_brm_ot_non_nv_is_is_n_pls,
		f_brm_ot_nv_is_is_n_pls,
	)

	var theta_r_ntr_non_nv_is_n_pls, theta_r_ntr_nv_is_n_pls mat.VecDense
	theta_r_ntr_non_nv_is_n_pls.MulVec(ss.f_xot_is_is_n_pls, theta_r_ot_ntr_non_nv_is_n_pls)
	theta_r_ntr_non_nv_is_n_pls.SubVec(&theta_r_ntr_non_nv_is_n_pls, f_xc_is_n_pls)
	theta_r_ntr_nv_is_n_pls.MulVec(ss.f_xot_is_is_n_pls, theta_r_ot_ntr_nv_is_n_pls)
	theta_r_ntr_nv_is_n_pls.SubVec(&theta_r_ntr_nv_is_n_pls, f_xc_is_n_pls)

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

	// 式(1.4)
	// ステップn+1における室iの自然風非利用時の加湿・除湿を行わない場合の絶対湿度, kg/kg(DA) [i, 1]
	// ステップn+1における室iの自然風利用時の加湿・除湿を行わない場合の絶対湿度, kg/kg(DA) [i, 1]
	x_r_ntr_non_nv_is_n_pls, x_r_ntr_nv_is_n_pls := get_x_r_ntr_is_n_pls(
		f_h_cst_non_nv_is_n,
		f_h_wgt_non_nv_is_is_n,
		f_h_cst_nv_is_n,
		f_h_wgt_nv_is_is_n,
	)

	// ステップ n における室 i の運転モード, [i, 1]
	operation_mode_is_n, all_stop := self.op.get_operation_mode_is_n(
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
	f_h_cst_is_n := mat.NewVecDense(self.rms.n_rm, nil)
	f_h_wgt_is_is_n := mat.NewDense(self.rms.n_rm, self.rms.n_rm, nil)
	theta_r_ot_ntr_is_n_pls := make([]float64, self.rms.n_rm)
	theta_r_ntr_is_n_pls := make([]float64, self.rms.n_rm)
	theta_mrt_hum_ntr_is_n_pls := make([]float64, self.rms.n_rm)
	x_r_ntr_is_n_pls := make([]float64, self.rms.n_rm)
	for i := 0; i < self.rms.n_rm; i++ {
		if operation_mode_is_n[i] == STOP_OPEN {
			f_brm_is_is_n_pls.SetRow(i, f_brm_nv_is_is_n_pls.RawRowView(i))
			f_brm_ot_is_is_n_pls.SetRow(i, f_brm_ot_nv_is_is_n_pls.RawRowView(i))
			f_h_wgt_is_is_n.SetRow(i, f_h_wgt_nv_is_is_n.RawRowView(i))
			v_vent_ntr_is_n[i] = self.rms.v_vent_ntr_set_is[i]
			f_brc_ot_is_n_pls[i] = f_brc_ot_nv_is_n_pls.AtVec(i)
			f_h_cst_is_n.SetVec(i, f_h_cst_nv_is_n.AtVec(i))
			theta_r_ot_ntr_is_n_pls[i] = theta_r_ot_ntr_nv_is_n_pls.AtVec(i)
			theta_r_ntr_is_n_pls[i] = theta_r_ntr_nv_is_n_pls.AtVec(i)
			theta_mrt_hum_ntr_is_n_pls[i] = theta_mrt_hum_ntr_nv_is_n_pls.AtVec(i)
			x_r_ntr_is_n_pls[i] = x_r_ntr_nv_is_n_pls.AtVec(i)
		} else {
			f_brm_is_is_n_pls.SetRow(i, f_brm_non_nv_is_is_n_pls.RawRowView(i))
			f_brm_ot_is_is_n_pls.SetRow(i, f_brm_ot_non_nv_is_is_n_pls.RawRowView(i))
			f_h_wgt_is_is_n.SetRow(i, f_h_wgt_non_nv_is_is_n.RawRowView(i))
			v_vent_ntr_is_n[i] = 0.0
			f_brc_ot_is_n_pls[i] = f_brc_ot_non_nv_is_n_pls.AtVec(i)
			f_h_cst_is_n.SetVec(i, f_h_cst_non_nv_is_n.AtVec(i))
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

	// ---- ここから放射暖冷房の計算 ----

	var beta_is_n *mat.VecDense
	var f_flr_js_is_n *mat.VecDense
	var f_wsb_js_is_n_pls *mat.Dense
	var f_xlr_is_is_n_pls *mat.Dense
	var f_brl_ot_is_is_n *mat.Dense
	if !all_stop && self.es.has_flr {
		// 放射暖冷房設備が設置の場合 and 運転

		// 式(2.14)
		// ステップ n からステップ n+1 における室 i の放射暖冷房設備の放熱量の放射成分に対する境界 j の室内側表面の吸収比率 f_flr, -, [j, i]
		// NOTE: 放射暖冷房設備が未設置の場合 or 運転時間外の場合はゼロ行列になる
		f_flr_js_is_n := get_f_flr_js_is_n(
			self.es.f_flr_c_js_is,
			self.es.f_flr_h_js_is,
			operation_mode_is_n,
		)

		// 式(2.13)
		// ステップ n からステップ n+1 における室 i の放射暖冷房設備の対流成分比率, -, [i, 1]
		// NOTE: 放射暖冷房設備が未設置の場合 or 運転時間外の場合はゼロ行列になる
		beta_is_n = get_beta_is_n(
			self.es.beta_c_is,
			self.es.beta_h_is,
			operation_mode_is_n,
		)

		// 式(2.12)
		// ステップ n における係数 f_WSB, K/W, [j, i]
		// ステップ n における係数 f_FLB, K/W, [j, i]
		// NOTE: f_flr_js_is_n の全ての要素が0である場合、常に結果はゼロである。
		//       すなわち、放射暖冷房設備が未設置の場合 or 運転時間外の場合はゼロ行列になる
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

		// 式(2.11)
		// ステップ n+1 における係数 f_WSB, K/W, [j, i]
		// NOTE: f_flb_js_is_n_plsがゼロ行列の場合、 f_WSBもまたゼロ行列となる。
		//       すなわち、放射暖冷房設備が未設置の場合 or 運転時間外の場合はゼロ行列になる
		f_wsb_js_is_n_pls = get_f_wsb_js_is_n_pls(
			f_flb_js_is_n_pls,
			ss.f_ax_js_js_inv,
			ss.f_ax_js_js,
		)

		// 式(2.10)
		// ステップ n における係数 f_BRL, -, [i, i]
		// NOTE: f_WSBがゼロ行列の場合、f_BRLは beta_is_nの対角化行列と等しい。
		//       f_WSBがゼロ行列であるのは、放射暖冷房設備が未設置の場合 or 運転時間外の場合である。
		//       この場合は、beta_is_nもまたゼロ行列である。
		//       結果的に、放射暖冷房設備が未設置の場合 or 運転時間外の場合にはf_BRLもゼロ行列になる。
		f_brl_is_is_n := get_f_brl_is_is_n(
			self.bs.a_s_js,
			beta_is_n,
			f_wsb_js_is_n_pls,
			self.bs.h_s_c_js,
			self.bs.p_is_js,
		)

		// 式(2.9)
		// ステップn+1における室iの係数 f_XLR, K/W, [i, i]
		// NOTE: 放射暖冷房設備が未設置の場合 or 運転時間外の場合にはf_XLR もゼロ行列になる。
		f_xlr_is_is_n_pls = get_f_xlr_is_is_n_pls(
			ss.f_mrt_hum_is_js,
			f_wsb_js_is_n_pls,
			ss.f_xot_is_is_n_pls,
			ss.k_r_is_n,
		)

		// 式(2.8)
		// ステップ n における係数 f_BRL_OT, -, [i, i]
		// NOTE: 放射暖冷房設備が未設置の場合 or 運転時間外の場合にはf_BRL_ot もゼロ行列になる。
		//       なぜなら、f_BRLとf_XLR がゼロであるため、 f_BRM * f_XLR + f_BRL はゼロである。
		f_brl_ot_is_is_n = get_f_brl_ot_is_is_n(
			f_brl_is_is_n,
			f_brm_is_is_n_pls,
			f_xlr_is_is_n_pls,
		)
	} else {
		f_brl_ot_is_is_n = nil
	}

	// ---- ここまで放射暖冷房の計算 ----

	// ---- ここから温度計算 ----

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

	// 式(2.6) 室温計算
	// ステップ n+1 における室 i の室温, degree C, [i, 1]
	theta_r_is_n_pls := get_theta_r_is_n_pls(
		f_xc_is_n_pls,
		f_xlr_is_is_n_pls,
		ss.f_xot_is_is_n_pls,
		l_rs_is_n,
		theta_ot_is_n_pls,
	)

	// 式(2.5) 表面温度計算
	// ステップ n+1 における境界 j の表面温度, degree C, [j, 1]
	theta_s_js_n_pls := get_theta_s_js_n_pls(
		f_wsb_js_is_n_pls,
		ss.f_wsc_js_ns.ColView(nn_plus+1),
		ss.f_wsr_js_is,
		f_wsv_js_n_pls,
		l_rs_is_n,
		theta_r_is_n_pls,
	)

	// 式(2.4) 備品等の温度計算
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

	// 式(2.3) 人体の平均放射温度計算
	// ステップ n+1 における室 i の人体に対する平均放射温度, degree C, [i, 1]
	theta_mrt_hum_is_n_pls := get_theta_mrt_hum_is_n_pls(
		ss.f_mrt_hum_is_js,
		theta_s_js_n_pls,
	)

	// 式(2.2) 等価温度計算
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

	// 式(2.1) 表面熱流計算
	// ステップ n+1 における境界 j の表面熱流（壁体吸熱を正とする）, W/m2, [j, 1]
	q_s_js_n_pls := get_q_s_js_n_pls(
		self.bs.h_s_c_js,
		self.bs.h_s_r_js,
		theta_ei_js_n_pls,
		theta_s_js_n_pls,
	)

	// ---- ここまで温度計算 ----

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
	c_n.theta_r_is_n.CopyVec(theta_r_is_n_pls)
	c_n.theta_mrt_hum_is_n = theta_mrt_hum_is_n_pls
	c_n.x_r_is_n = x_r_is_n_pls
	c_n.theta_dsh_srf_a_js_ms_n = theta_dsh_s_a_js_ms_n_pls
	c_n.theta_dsh_srf_t_js_ms_n = theta_dsh_s_t_js_ms_n_pls
	c_n.q_s_js_n = q_s_js_n_pls
	c_n.theta_frt_is_n = theta_frt_is_n_pls
	c_n.x_frt_is_n = x_frt_is_n_pls
	c_n.theta_ei_js_n = theta_ei_js_n_pls

	return c_n
}

var __theta_s_js_npls []float64
var __theta_dsh_srf_a_js_ms_npls *mat.Dense
var __theta_dsh_srf_t_js_ms_npls *mat.Dense
var __theta_o_eqv_js_ns *mat.Dense
var __h_i_js *mat.VecDense
var __ground_idx []int

// 地盤の計算（n+1ステップを計算する）
func _run_tick_ground(self *Sequence, pp *PreCalcParameters, gc_n *GroundConditions, n int, nn int) *GroundConditions {
	if gc_n == nil {
		return nil
	}

	if __ground_idx == nil {
		__ground_idx = make([]int, 0, self.bs.n_b)

		is_ground := self.bs.is_ground_js
		for i, v := range is_ground {
			if v {
				__ground_idx = append(__ground_idx, i)
			}
		}
	}
	ground_idx := __ground_idx
	_, l := self.bs.theta_o_eqv_js_ns.Dims()

	if __theta_o_eqv_js_ns == nil {
		__theta_o_eqv_js_ns = mat.NewDense(len(ground_idx), l, nil)
	}
	if __theta_dsh_srf_a_js_ms_npls == nil {
		__theta_dsh_srf_a_js_ms_npls = mat.NewDense(len(ground_idx), 12, nil)
	}
	if __theta_dsh_srf_t_js_ms_npls == nil {
		__theta_dsh_srf_t_js_ms_npls = mat.NewDense(len(ground_idx), 12, nil)
	}
	if __theta_s_js_npls == nil {
		__theta_s_js_npls = make([]float64, len(ground_idx))
	}
	if __h_i_js == nil {
		__h_i_js = mat.NewVecDense(len(ground_idx), nil)
	}

	theta_o_eqv_js_ns := __theta_o_eqv_js_ns
	for i := 0; i < len(ground_idx); i++ {
		theta_o_eqv_js_ns.SetRow(i, self.bs.theta_o_eqv_js_ns.RawRowView(ground_idx[i]))
	}

	h_i_js := __h_i_js
	for i := 0; i < len(ground_idx); i++ {
		gidx := ground_idx[i]
		h_i_js.SetVec(i, self.bs.h_s_r_js.AtVec(gidx)+self.bs.h_s_c_js.AtVec(gidx))
	}

	theta_dsh_srf_a_js_ms_npls := __theta_dsh_srf_a_js_ms_npls
	for j := 0; j < len(ground_idx); j++ {
		gidx := ground_idx[j]
		for i := 0; i < 12; i++ {
			theta_dsh_srf_a_js_ms_npls.Set(j, i, self.bs.phi_a1_js_ms.At(gidx, i)*gc_n.q_srf_js_n[j]+
				self.bs.r_js_ms.At(gidx, i)*gc_n.theta_dsh_srf_a_js_ms_n.At(j, i))
		}
	}

	theta_dsh_srf_t_js_ms_npls := __theta_dsh_srf_t_js_ms_npls
	for j := 0; j < len(ground_idx); j++ {
		gidx := ground_idx[j]
		for i := 0; i < 12; i++ {
			theta_dsh_srf_t_js_ms_npls.Set(j, i,
				self.bs.phi_t1_js_ms.At(gidx, i)*self.bs.k_eo_js.AtVec(gidx)*self.bs.theta_o_eqv_js_ns.At(gidx, nn)+
					self.bs.r_js_ms.At(gidx, i)*gc_n.theta_dsh_srf_t_js_ms_n.At(j, i))
		}
	}

	theta_s_js_npls := __theta_s_js_npls
	for j := 0; j < len(ground_idx); j++ {
		sum_a := floats.Sum(theta_dsh_srf_a_js_ms_npls.RawRowView(j))
		sum_t := floats.Sum(theta_dsh_srf_t_js_ms_npls.RawRowView(j))

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

// ----------------------------------------------------------------------------------
// 3.1 湿度と潜熱処理量
// ----------------------------------------------------------------------------------

var __x_frt_is_n_pls__temp1 mat.VecDense
var __x_frt_is_n_pls__temp2 mat.VecDense
var __x_frt_is_n_pls__temp3 mat.VecDense

/*
備品等の絶対湿度を求める。

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
	c_lh_frt_is *mat.VecDense,
	delta_t float64,
	g_lh_frt_is *mat.VecDense,
	x_frt_is_n *mat.VecDense,
	x_r_is_n_pls *mat.VecDense,
) *mat.VecDense {
	temp1 := &__x_frt_is_n_pls__temp1
	temp2 := &__x_frt_is_n_pls__temp2
	temp3 := &__x_frt_is_n_pls__temp3

	// c_lh_frt_is * x_frt_is_n
	temp1.MulElemVec(c_lh_frt_is, x_frt_is_n)

	// g_lh_frt_is * x_r_is_n_pls
	temp2.MulElemVec(g_lh_frt_is, x_r_is_n_pls)

	// c_lh_frt_is * x_frt_is_n + delta_t * (g_lh_frt_is * x_r_is_n_pls)
	temp1.AddScaledVec(temp1, delta_t, temp2)

	// c_lh_frt_is + delta_t * g_lh_frt_is
	temp3.AddScaledVec(c_lh_frt_is, delta_t, g_lh_frt_is)

	// temp1 / temp3
	temp1.DivElemVec(temp1, temp3)

	return temp1
}

var __l_cl_is_n mat.Dense

/*
対流暖冷房設備の潜熱処理量を求める。

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
	f_l_cl_cst_is_n *mat.VecDense,
	l_wtr float64,
	x_r_is_n_pls *mat.VecDense,
) *mat.Dense {
	temp1 := &__l_cl_is_n
	temp1.Mul(f_l_cl_wgt_is_is_n, x_r_is_n_pls)
	temp1.Add(temp1, f_l_cl_cst_is_n)
	temp1.Scale(l_wtr, temp1)

	return temp1
}

var __x_r_is_n_pls__temp1 mat.VecDense
var __x_r_is_n_pls__temp2 mat.Dense
var __x_r_is_n_pls__result mat.VecDense

/*
絶対湿度を求める。

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
	f_h_cst_is_n *mat.VecDense,
	f_h_wgt_is_is_n *mat.Dense,
	f_l_cl_cst_is_n *mat.VecDense,
	f_l_cl_wgt_is_is_n *mat.Dense,
) *mat.VecDense {
	temp1 := &__x_r_is_n_pls__temp1
	temp2 := &__x_r_is_n_pls__temp2
	result := &__x_r_is_n_pls__result

	//f_h_cst_is_n + f_l_cl_cst_is_n
	temp1.AddVec(f_h_cst_is_n, f_l_cl_cst_is_n)

	// f_h_wgt_is_is_n - f_l_cl_wgt_is_is_n
	temp2.Sub(f_h_wgt_is_is_n, f_l_cl_wgt_is_is_n)

	result.SolveVec(temp2, temp1)

	return result
}

var __x_r_ntr_is_n_pls__result1 mat.VecDense
var __x_r_ntr_is_n_pls__result2 mat.VecDense

/*
加湿・除湿を行わない場合の絶対湿度を求める。

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
	f_h_cst_non_nv_is_n *mat.VecDense,
	f_h_wgt_non_nv_is_is_n *mat.Dense,
	f_h_cst_nv_is_n *mat.VecDense,
	f_h_wgt_nv_is_is_n *mat.Dense,
) (*mat.VecDense, *mat.VecDense) {
	result1 := &__x_r_ntr_is_n_pls__result1
	result2 := &__x_r_ntr_is_n_pls__result2
	result1.SolveVec(f_h_wgt_non_nv_is_is_n, f_h_cst_non_nv_is_n)
	result2.SolveVec(f_h_wgt_nv_is_is_n, f_h_cst_nv_is_n)
	return result1, result2
}

var __f_h_wgt_is_is_n__temp1 mat.VecDense
var __f_h_wgt_is_is_n__temp2 mat.VecDense
var __f_h_wgt_is_is_n__temp3 mat.VecDense
var __f_h_wgt_is_is_n__temp5 mat.Dense
var __f_h_wgt_is_is_n__result1 *mat.Dense
var __f_h_wgt_is_is_n__result2 mat.Dense

/*
係数 f_h,wgt を求める

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
	c_lh_frt_is *mat.VecDense,
	delta_t float64,
	g_lh_frt_is *mat.VecDense,
	v_rm_is *mat.VecDense,
	v_vent_int_is_is_n mat.Matrix,
	v_vent_out_is_n *mat.VecDense,
	v_vent_ntr_is []float64,
) (*mat.Dense, *mat.Dense) {
	temp1 := &__f_h_wgt_is_is_n__temp1
	temp2 := &__f_h_wgt_is_is_n__temp2
	temp3 := &__f_h_wgt_is_is_n__temp3
	temp5 := &__f_h_wgt_is_is_n__temp5
	result2 := &__f_h_wgt_is_is_n__result2

	// rho_a*(v_rm_is/delta_t+v_vent_out_is_n)
	temp1.AddScaledVec(v_vent_out_is_n, 1/delta_t, v_rm_is)
	temp1.ScaleVec(rho_a, temp1)

	// c_lh_frt_is + delta_t * g_lh_frt_is)
	temp2.AddScaledVec(c_lh_frt_is, delta_t, g_lh_frt_is)

	// c_lh_frt_is * g_lh_frt_is / (c_lh_frt_is + delta_t * g_lh_frt_is)
	temp3.MulElemVec(c_lh_frt_is, g_lh_frt_is)
	temp3.DivElemVec(temp3, temp2)

	// result1 = diag of temp1 + temp3
	temp1.AddVec(temp1, temp3)
	if __f_h_wgt_is_is_n__result1 == nil {
		l := temp1.Len()
		__f_h_wgt_is_is_n__result1 = mat.NewDense(l, l, nil)
	}
	result1 := __f_h_wgt_is_is_n__result1
	for i := 0; i < temp1.Len(); i++ {
		result1.Set(i, i, temp1.AtVec(i))
	}

	result1.Apply(func(i, j int, v float64) float64 {
		return v - rho_a*v_vent_int_is_is_n.At(i, j)
	}, result1)

	temp5.Scale(rho_a, NewDiagAsDenseFromFloat64(v_vent_ntr_is))
	result2.Add(result1, temp5)

	return result1, result2
}

var __f_h_cst_is_n__temp1 mat.VecDense
var __f_h_cst_is_n__temp2 mat.VecDense
var __f_h_cst_is_n__temp3 mat.VecDense
var __f_h_cst_is_n__temp4 mat.VecDense
var __f_h_cst_is_n__result1 mat.VecDense
var __f_h_cst_is_n__result2 mat.VecDense

/*
係数 f_h,cst を求める

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
	c_lh_frt_is *mat.VecDense,
	delta_t float64,
	g_lh_frt_is *mat.VecDense,
	rho_a float64,
	v_rm_is *mat.VecDense,
	x_frt_is_n *mat.VecDense,
	x_gen_is_n mat.Vector,
	x_hum_is_n *mat.VecDense,
	x_o_n_pls float64,
	x_r_is_n *mat.VecDense,
	v_vent_out_non_nv_is_n *mat.VecDense,
	v_vent_ntr_is []float64,
) (*mat.VecDense, *mat.VecDense) {
	temp1 := &__f_h_cst_is_n__temp1
	temp2 := &__f_h_cst_is_n__temp2
	temp3 := &__f_h_cst_is_n__temp3
	temp4 := &__f_h_cst_is_n__temp4
	result1 := &__f_h_cst_is_n__result1
	result2 := &__f_h_cst_is_n__result2

	// Python: rho_a * v_rm_is / delta_t * x_r_is_n
	temp1.ScaleVec(rho_a/delta_t, v_rm_is)
	temp1.MulElemVec(temp1, x_r_is_n)

	// Python: rho_a * v_vent_out_non_nv_is_n * x_o_n_pls
	temp2.ScaleVec(rho_a*x_o_n_pls, v_vent_out_non_nv_is_n)

	// Python: c_lh_frt_is * g_lh_frt_is / (c_lh_frt_is + delta_t * g_lh_frt_is) * x_frt_is_n
	temp3.MulElemVec(c_lh_frt_is, g_lh_frt_is)            //(c_lh_frt_is * g_lh_frt_is)
	temp4.AddScaledVec(c_lh_frt_is, delta_t, g_lh_frt_is) //(c_lh_frt_is + delta_t * g_lh_frt_is)
	temp3.DivElemVec(temp3, temp4)                        //(c_lh_frt_is * g_lh_frt_is / (c_lh_frt_is + delta_t * g_lh_frt_is)
	temp3.MulElemVec(temp3, x_frt_is_n)

	// Python: temp1 + temp2 + temp3 + x_gen_is_n + x_hum_is_n
	result1.AddVec(temp1, temp2)
	result1.AddVec(result1, temp3)
	result1.AddVec(result1, x_gen_is_n)
	result1.AddVec(result1, x_hum_is_n)

	result2.AddScaledVec(result1, rho_a*x_o_n_pls, mat.NewVecDense(len(v_vent_ntr_is), v_vent_ntr_is))

	return result1, result2
}

var __x_hum_is_n mat.VecDense

/*
人体発湿を求める。
Args:

	n_hum_is_n: ステップ n からステップ n+1 における室 i の在室人数, -
	x_hum_psn_is_n: ステップ n からステップ n+1 における室 i の1人あたりの人体発湿, kg/s

Returns:

	ステップnの室iにおける人体発湿, kg/s, [i, 1]

Notes:

	式(1.7)
*/
func get_x_hum_is_n(n_hum_is_n *mat.VecDense, x_hum_psn_is_n *mat.VecDense) *mat.VecDense {
	x_hum_is_n := &__x_hum_is_n
	x_hum_is_n.MulElemVec(n_hum_is_n, x_hum_psn_is_n)
	return x_hum_is_n
}

// ----------------------------------------------------------------------------------
// 3.2 温度と顕熱処理量
// ----------------------------------------------------------------------------------

var __q_s_js_n_pls__temp1 mat.VecDense
var __q_s_js_n_pls__temp2 mat.VecDense

/*
表面熱流を求める。

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
	h_s_c_js *mat.VecDense,
	h_s_r_js *mat.VecDense,
	theta_ei_js_n_pls *mat.VecDense,
	theta_s_js_n_pls *mat.VecDense,
) *mat.VecDense {
	temp1 := &__q_s_js_n_pls__temp1
	temp2 := &__q_s_js_n_pls__temp2

	temp1.SubVec(theta_ei_js_n_pls, theta_s_js_n_pls)
	temp2.AddVec(h_s_c_js, h_s_r_js)
	temp1.MulElemVec(temp1, temp2)

	return temp1
}

var __theta_ei_js_n_pls__temp1 mat.VecDense
var __theta_ei_js_n_pls__temp2 mat.Dense
var __theta_ei_js_n_pls__temp3 mat.VecDense
var __theta_ei_js_n_pls__temp4 mat.VecDense
var __theta_ei_js_n_pls__temp5 mat.VecDense
var __theta_ei_js_n_pls__temp6 mat.VecDense

/*
等価温度を求める。

Args:

	a_s_js: 境界 j の面積, m2, [j, 1]
	beta_is_n: ステップ n からステップ n+1 における室 i の放射暖冷房設備の対流成分比率, -, [i, 1]
	f_mrt_is_js_T: 室 i の微小球に対する境界 j の形態係数(転地済み), -, [i, j]
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
	a_s_js *mat.VecDense,
	beta_is_n *mat.VecDense,
	f_mrt_is_js mat.Matrix,
	f_flr_js_is_n *mat.VecDense,
	h_s_c_js *mat.VecDense,
	h_s_r_js *mat.VecDense,
	l_rs_is_n *mat.VecDense,
	p_js_is mat.Matrix,
	q_s_sol_js_n_pls mat.Vector,
	theta_r_is_n_pls *mat.VecDense,
	theta_s_js_n_pls *mat.VecDense,
) *mat.VecDense {
	// Python: h_s_c_js*np.dot(p_js_is, theta_r_is_n_pls)
	temp1 := &__theta_ei_js_n_pls__temp1
	temp1.MulVec(p_js_is, theta_r_is_n_pls)
	temp1.MulElemVec(h_s_c_js, temp1)

	// Python: h_s_r_js*np.dot(np.dot(p_js_is, f_mrt_is_js), theta_s_js_n_pls)
	temp2 := &__theta_ei_js_n_pls__temp2
	temp3 := &__theta_ei_js_n_pls__temp3
	temp2.Mul(p_js_is, f_mrt_is_js)       // temp2 <= np.dot(p_js_is, f_mrt_is_js)
	temp3.MulVec(temp2, theta_s_js_n_pls) // temp3 <= np.dot(temp2, theta_s_js_n_pls)
	temp3.MulElemVec(h_s_r_js, temp3)     // temp3 <= h_s_r_js * temp3

	// Python: np.dot(f_flr_js_is_n, (1.0-beta_is_n)*l_rs_is_n)/a_s_js
	temp4 := &__theta_ei_js_n_pls__temp4
	temp5 := &__theta_ei_js_n_pls__temp5
	if beta_is_n != nil && f_flr_js_is_n != nil {
		temp4.MulElemVec(beta_is_n, l_rs_is_n)
		temp4.SubVec(l_rs_is_n, temp4)
		temp5.MulVec(f_flr_js_is_n, temp4)
		temp4.SubVec(l_rs_is_n, temp4)
		temp5.MulVec(f_flr_js_is_n, temp4)
		temp5.DivElemVec(temp5, a_s_js)
	} else {
		temp4 = nil
		temp5 = nil
	}

	// temp1 <= (temp1 + temp3 + q_s_sol_js_n_pls + temp5) / (h_s_c_js + h_s_r_js)
	temp1.AddVec(temp1, temp3)
	temp1.AddVec(temp1, q_s_sol_js_n_pls)
	if temp5 != nil {
		temp1.AddVec(temp1, temp5)
	}

	temp6 := &__theta_ei_js_n_pls__temp6
	temp6.AddVec(h_s_c_js, h_s_r_js)
	temp1.DivElemVec(temp1, temp6)

	return temp1
}

var __theta_mrt_hum_is_n_pls__resut mat.VecDense

/*
人体の平均放射温度を求める。

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
	theta_s_js_n_pls *mat.VecDense,
) *mat.VecDense {
	result := &__theta_mrt_hum_is_n_pls__resut
	result.MulVec(f_mrt_hum_is_js, theta_s_js_n_pls)
	return result
}

var __theta_frt_is_n_pls__temp1 mat.VecDense
var __theta_frt_is_n_pls__temp2 mat.VecDense

/*
備品等の温度を求める。

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
	c_sh_frt_is *mat.VecDense,
	delta_t float64,
	g_sh_frt_is *mat.VecDense,
	q_sol_frt_is_n mat.Vector,
	theta_frt_is_n *mat.VecDense,
	theta_r_is_n_pls *mat.VecDense,
) *mat.VecDense {
	temp1 := &__theta_frt_is_n_pls__temp1
	temp2 := &__theta_frt_is_n_pls__temp2

	// temp1 <= c_sh_frt_is*theta_frt_is_n
	temp1.MulElemVec(c_sh_frt_is, theta_frt_is_n)

	// temp2 <= g_sh_frt_is*theta_r_is_n_pls
	temp2.MulElemVec(g_sh_frt_is, theta_r_is_n_pls)

	// temp1 <= (temp1 + delta_t*temp2 + q_sol_frt_is_n*delta_t)
	temp1.AddScaledVec(temp1, delta_t, temp2)
	temp1.AddScaledVec(temp1, delta_t, q_sol_frt_is_n)

	// temp2 <= c_sh_frt_is + delta_t*g_sh_frt_is
	temp2.AddScaledVec(c_sh_frt_is, delta_t, g_sh_frt_is)

	temp1.DivElemVec(temp1, temp2)

	return temp1
}

var __theta_s_js_n_pls__temp1 mat.VecDense
var __theta_s_js_n_pls__temp2 mat.VecDense

/*
表面温度を求める。

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
	f_wsb_js_is_n_pls *mat.Dense,
	f_wsc_js_n_pls mat.Vector,
	f_wsr_js_is mat.Matrix,
	f_wsv_js_n_pls *mat.VecDense,
	l_rs_is_n *mat.VecDense,
	theta_r_is_n_pls *mat.VecDense,
) *mat.VecDense {
	temp1 := &__theta_s_js_n_pls__temp1
	temp2 := &__theta_s_js_n_pls__temp2

	temp1.MulVec(f_wsr_js_is, theta_r_is_n_pls) //np.dot(f_wsr_js_is, theta_r_is_n_pls)
	temp1.AddVec(temp1, f_wsc_js_n_pls)
	temp1.AddVec(temp1, f_wsv_js_n_pls)

	if f_wsb_js_is_n_pls != nil {
		temp2.MulVec(f_wsb_js_is_n_pls, l_rs_is_n) //np.dot(f_wsb_js_is_n_pls, l_rs_is_n)
		temp1.AddVec(temp1, temp2)
	}

	return temp1
}

var __theta_r_is_n_pls mat.VecDense

/*
室温を求める。

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
	f_xc_is_n_pls *mat.VecDense,
	f_xlr_is_is_n_pls *mat.Dense,
	f_xot_is_is_n_pls mat.Matrix,
	l_rs_is_n *mat.VecDense,
	theta_ot_is_n_pls *mat.VecDense,
) *mat.VecDense {
	// Python: np.dot(f_xot_is_is_n_pls, theta_ot_is_n_pls)
	temp1 := &__theta_r_is_n_pls
	temp1.MulVec(f_xot_is_is_n_pls, theta_ot_is_n_pls)

	if f_xlr_is_is_n_pls != nil {
		// Python: np.dot(f_xlr_is_is_n_pls, l_rs_is_n)
		var temp2 mat.VecDense
		temp2.MulVec(f_xlr_is_is_n_pls, l_rs_is_n)

		temp1.SubVec(temp1, &temp2)
	}

	// Python: temp1 - temp2 - f_xc_is_n_pls
	temp1.SubVec(temp1, f_xc_is_n_pls)

	return temp1
}

var __f_brl_ot_is_is_n mat.Dense

/*
係数 f^_BRL,OT を求める。

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
	result := &__f_brl_ot_is_is_n
	result.Mul(f_brm_is_is_n_pls, f_xlr_is_is_n_pls)
	result.Add(f_brl_is_is_n, result)
	return result
}

var __f_xlr_is_is_n_pls__temp1 mat.Dense
var __f_xlr_is_is_n_pls__result mat.Dense

/*
係数 f_XLR を求める。

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
	k_r_is_n *mat.VecDense,
) *mat.Dense {
	// Python: np.dot(f_mrt_hum_is_js, f_wsb_js_is_n_pls)
	temp1 := &__f_xlr_is_is_n_pls__temp1
	temp1.Mul(f_mrt_hum_is_js, f_wsb_js_is_n_pls)

	// Python: k_r_is_n * temp1
	__ScaleRows(temp1, k_r_is_n)

	// Python: np.dot(f_xot_is_is_n_pls, temp2)
	result := &__f_xlr_is_is_n_pls__result
	result.Mul(f_xot_is_is_n_pls, temp1)

	return result
}

var __f_wsb_js_is_n_pls_h_s_c_js_a_s_js *mat.Dense
var __p_is_js_f_wsb_js_is_n_pls_h_s_c_js_a_s_js mat.Dense
var __f_brl_is_is_n__temp mat.VecDense
var __f_brl_is_is_n mat.Dense
var __beta_is_n_diag *mat.Dense

/*
係数 f_BRL を求める。

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
	a_s_js *mat.VecDense,
	beta_is_n *mat.VecDense,
	f_wsb_js_is_n_pls *mat.Dense,
	h_s_c_js *mat.VecDense,
	p_is_js *mat.Dense,
) mat.Matrix {
	// v_diag(beta_is_n)
	// ^Beta,n は ^Beta,i,n を要素に持つ対角化行列
	if __beta_is_n_diag == nil {
		__beta_is_n_diag = mat.NewDense(beta_is_n.Len(), beta_is_n.Len(), nil)
	}
	beta_is_n_diag := __beta_is_n_diag
	Diag(beta_is_n_diag, beta_is_n)

	// f_WSB,n+1 がゼロ行列である場合、 f_BRL,n = ^Beta,n
	if f_wsb_js_is_n_pls == nil {
		return beta_is_n_diag
	}

	// f_wsb_js_is_n_pls*h_s_c_js*a_s_js
	temp := &__f_brl_is_is_n__temp
	temp.AddVec(h_s_c_js, a_s_js)

	if __f_wsb_js_is_n_pls_h_s_c_js_a_s_js == nil {
		r, c := f_wsb_js_is_n_pls.Dims()
		__f_wsb_js_is_n_pls_h_s_c_js_a_s_js = mat.NewDense(r, c, nil)
	}
	f_wsb_js_is_n_pls_h_s_c_js_a_s_js := __f_wsb_js_is_n_pls_h_s_c_js_a_s_js

	__ScaleRowsTo(f_wsb_js_is_n_pls_h_s_c_js_a_s_js, temp, f_wsb_js_is_n_pls)

	// np.dot(p_is_js, f_wsb_js_is_n_pls*h_s_c_js*a_s_js)
	p_is_js_f_wsb_js_is_n_pls_h_s_c_js_a_s_js := &__p_is_js_f_wsb_js_is_n_pls_h_s_c_js_a_s_js
	p_is_js_f_wsb_js_is_n_pls_h_s_c_js_a_s_js.Mul(p_is_js, f_wsb_js_is_n_pls_h_s_c_js_a_s_js)

	// np.dot(p_is_js, f_wsb_js_is_n_pls*h_s_c_js*a_s_js) + v_diag(beta_is_n)
	result := &__f_brl_is_is_n
	result.Add(p_is_js_f_wsb_js_is_n_pls_h_s_c_js_a_s_js, beta_is_n_diag)

	return result
}

func __ScaleRows(x *mat.Dense, a *mat.VecDense) {
	mat := x.RawMatrix()

	r := mat.Rows
	c := mat.Cols

	for i := 0; i < r; i++ {
		floats.Scale(a.AtVec(i), mat.Data[i*mat.Stride:i*mat.Stride+c])
	}
}

func __ScaleRowsTo(dst *mat.Dense, a *mat.VecDense, src *mat.Dense) {
	dst_mat := dst.RawMatrix()
	src_mat := src.RawMatrix()

	dr := dst_mat.Rows
	dc := dst_mat.Cols

	sr := src_mat.Rows
	sc := src_mat.Cols

	if dr != sr || dc != sc {
		panic("Invalid Matrix dimensions")
	}

	_a := a.RawVector().Data
	for i := 0; i < dr; i++ {
		floats.ScaleTo(dst_mat.Data[i*dst_mat.Stride:i*dst_mat.Stride+dc], _a[i], src_mat.Data[i*src_mat.Stride:i*src_mat.Stride+sc])
	}
}

func __ScaleColumn(a *mat.VecDense, dst *mat.Dense) {
	mat := dst.RawMatrix()

	r := mat.Rows
	c := mat.Cols

	for i := 0; i < r; i++ {
		floats.MulTo(mat.Data[i*mat.Stride:i*mat.Stride+c], a.RawVector().Data, mat.Data[i*mat.Stride:i*mat.Stride+c])
	}
}

func __ScaleColumnTo(dst *mat.Dense, a *mat.VecDense, src *mat.Dense) {
	dst_mat := dst.RawMatrix()
	src_mat := src.RawMatrix()

	dr := dst_mat.Rows
	dc := dst_mat.Cols

	sr := src_mat.Rows
	sc := src_mat.Cols

	if dr != sr || dc != sc {
		panic("Invalid Matrix dimensions")
	}

	for i := 0; i < dr; i++ {
		floats.MulTo(dst_mat.Data[i*dst_mat.Stride:i*dst_mat.Stride+dc], a.RawVector().Data, src_mat.Data[i*src_mat.Stride:i*src_mat.Stride+sc])
	}
}

// dst = 1 - src
func __InvertVecTo(dst *mat.VecDense, src *mat.VecDense) {
	// dst = -src
	floats.ScaleTo(dst.RawVector().Data, -1.0, src.RawVector().Data)
	// dst = 1 + dst = 1 - src
	floats.AddConst(1.0, dst.RawVector().Data)
}

var __f_wsb_js_is_n_pls mat.Dense

/*
係数 f_WSB を求める。

Args:

	f_flb_js_is_n_pls: ステップ n+1 における係数 f_FLB, K/W, [j, i]
	f_ax_js_js_revert: 係数 f_AX^-1, -, [j, j]

Returns:

	ステップ n+1 における係数 f_WSB, K/W, [j, i]

Notes:

	式(2.11)
*/
func get_f_wsb_js_is_n_pls(
	f_flb_js_is_n_pls mat.Matrix,
	f_ax_js_js_inv *mat.Dense,
	f_ax_js_js *mat.LU,
) *mat.Dense {
	// 事前にLU分解した結果を用いて解を得る場合
	// ※総じて、LU分解のほうが安定して速度が出る
	f_ax_js_js.SolveTo(&__f_wsb_js_is_n_pls, false, f_flb_js_is_n_pls)

	// 事前にもとめた逆行列を使って解を得る場合
	//__f_wsb_js_is_n_pls.Mul(f_ax_js_js_inv, f_flb_js_is_n_pls)

	return &__f_wsb_js_is_n_pls
}

var __f_flb_js_is_n_pls_temp0 *mat.Dense
var __f_flb_js_is_n_pls_temp1 *mat.VecDense
var __f_flb_js_is_n_pls_temp2 mat.VecDense
var __f_flb_js_is_n_pls_temp3 mat.Dense
var __f_flb_js_is_n_pls_temp4 mat.VecDense

/*
係数 f_FLB を求める。

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
	beta_is_n の全ての要素が1である場合、常に結果はゼロである。
*/
func get_f_flb_js_is_n_pls(
	a_s_js *mat.VecDense,
	beta_is_n *mat.VecDense,
	f_flr_js_is_n *mat.Dense,
	h_s_c_js *mat.VecDense,
	h_s_r_js *mat.VecDense,
	k_ei_js_js mat.Matrix,
	phi_a0_js *mat.VecDense,
	phi_t0_js *mat.VecDense,
) mat.Matrix {
	if __f_flb_js_is_n_pls_temp0 == nil {
		r, c := f_flr_js_is_n.Dims()
		__f_flb_js_is_n_pls_temp0 = mat.NewDense(r, c, nil)
	}
	if __f_flb_js_is_n_pls_temp1 == nil {
		__f_flb_js_is_n_pls_temp1 = mat.NewVecDense(beta_is_n.Len(), nil)
	}

	term0 := __f_flb_js_is_n_pls_temp0
	term1 := __f_flb_js_is_n_pls_temp1

	// f_flr_js_is_n*(1.0-beta_is_n.T)
	__InvertVecTo(term1, beta_is_n) // (1.0-beta_is_n)  ※ .T が無いことに注意
	__ScaleColumnTo(term0, term1, f_flr_js_is_n)

	term2 := &__f_flb_js_is_n_pls_temp2
	// phi_a0_js/a_s_js
	term2.DivElemVec(phi_a0_js, a_s_js)

	// f_flr_js_is_n*(1.0-beta_is_n.T) * phi_a0_js/a_s_js
	__ScaleRows(term0, term2)

	term3 := &__f_flb_js_is_n_pls_temp3
	// np.dot(k_ei_js_js, f_flr_js_is_n*(1.0-beta_is_n.T))
	term3.Mul(k_ei_js_js, term0)

	term4 := &__f_flb_js_is_n_pls_temp4
	// phi_t0_js/(h_s_c_js+h_s_r_js)/a_s_js
	term4.AddVec(h_s_c_js, h_s_r_js)
	term4.DivElemVec(phi_t0_js, term4)
	term4.DivElemVec(term4, a_s_js)

	// np.dot(k_ei_js_js, f_flr_js_is_n*(1.0-beta_is_n.T))*phi_t0_js/(h_s_c_js+h_s_r_js)/a_s_js
	__ScaleRows(term3, term4)
	term0.Add(term0, term3)

	return term0
}

/*
放射暖冷房設備の対流成分比率を求める。

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
		if operation_mode_is_n[i] == HEATING {
			result.SetVec(i, beta_h_is[i])
		} else if operation_mode_is_n[i] == COOLING {
			result.SetVec(i, beta_c_is[i])
		} else {
			//PASS
		}
	}

	return result
}

var __f_flr_js_is_n *mat.Dense

/*
室内側表面の吸収比率を求める。

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
	if __f_flr_js_is_n == nil {
		__f_flr_js_is_n = mat.NewDense(r, c, nil)
	}

	result := __f_flr_js_is_n
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			if operation_mode_is_n[j] == HEATING {
				result.Set(i, j, f_flr_h_js_is.At(i, j))
			} else if operation_mode_is_n[j] == COOLING {
				result.Set(i, j, f_flr_c_js_is.At(i, j))
			} else {
				//PASS
			}
		}
	}
	return result
}

var __f_brc_ot_is_n_pls_temp1 mat.VecDense
var __f_brc_ot_is_n_pls_temp2 mat.VecDense

/*
係数 f^_BRC_OT を求める。

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
	f_xc_is_n_pls *mat.VecDense,
	f_brc_non_nv_is_n_pls *mat.VecDense,
	f_brc_nv_is_n_pls *mat.VecDense,
	f_brm_non_nv_is_is_n_pls mat.Matrix,
	f_brm_nv_is_is_n_pls mat.Matrix,
) (*mat.VecDense, *mat.VecDense) {
	temp1 := &__f_brc_ot_is_n_pls_temp1
	temp2 := &__f_brc_ot_is_n_pls_temp2

	// f_brc_non_nv_is_n_pls + np.dot(f_brm_non_nv_is_is_n_pls, f_xc_is_n_pls)
	temp1.MulVec(f_brm_non_nv_is_is_n_pls, f_xc_is_n_pls)
	temp1.AddVec(f_brc_non_nv_is_n_pls, temp1)

	// f_brc_nv_is_n_pls + np.dot(f_brm_nv_is_is_n_pls, f_xc_is_n_pls)
	temp2.MulVec(f_brm_nv_is_is_n_pls, f_xc_is_n_pls)
	temp2.AddVec(f_brc_nv_is_n_pls, temp2)

	return temp1, temp2
}

var __f_brm_ot_is_is_n_pls__result1 mat.Dense
var __f_brm_ot_is_is_n_pls__result2 mat.Dense

/*
係数 f^_BRM_OT を求める。

Args:

	f_xot_is_is_n_pls: ステップ n+1 における係数 f_XOT, -, [i, i]
	f_brm_non_nv_is_is_n_pls: ステップ n+1 における自然風非利用時の係数 f_BRM, W/K, [i, i]
	f_brm_nv_is_is_n_pls: ステップ n+1 における自然風利用時の係数 f_BRM, W/K, [i, i]

Returns:

	ステップn+1における自然風非利用時の係数f_BRM,OT, W/K, [i, i]
	ステップn+1における自然風利用時の係数f_BRM,OT, W/K, [i, i]

Notes:

	式(2.18)
*/
func get_f_brm_ot_is_is_n_pls(f_xot_is_is_n_pls mat.Matrix, f_brm_non_nv_is_is_n_pls mat.Matrix, f_brm_nv_is_is_n_pls mat.Matrix) (*mat.Dense, *mat.Dense) {
	temp1 := &__f_brm_ot_is_is_n_pls__result1
	temp2 := &__f_brm_ot_is_is_n_pls__result2
	temp1.Mul(f_brm_non_nv_is_is_n_pls, f_xot_is_is_n_pls)
	temp2.Mul(f_brm_nv_is_is_n_pls, f_xot_is_is_n_pls)
	return temp1, temp2
}

var __theta_r_ot_ntr_is_n_pls__result1 mat.VecDense
var __theta_r_ot_ntr_is_n_pls__result2 mat.VecDense

/*
自然作用温度を求める。

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
	f_brc_ot_non_nv_is_n_pls *mat.VecDense,
	f_brc_ot_nv_is_n_pls *mat.VecDense,
	f_brm_ot_non_nv_is_is_n_pls mat.Matrix,
	f_brm_ot_nv_is_is_n_pls mat.Matrix,
) (*mat.VecDense, *mat.VecDense) {

	result1 := &__theta_r_ot_ntr_is_n_pls__result1
	result2 := &__theta_r_ot_ntr_is_n_pls__result2

	result1.SolveVec(f_brm_ot_non_nv_is_is_n_pls, f_brc_ot_non_nv_is_n_pls)
	result2.SolveVec(f_brm_ot_nv_is_is_n_pls, f_brc_ot_nv_is_n_pls)

	return result1, result2
}

var __f_xc_is_n_pls__f_wscPlusWsv mat.VecDense
var __f_xc_is_n_pls__dotProduct mat.VecDense
var __f_xc_is_n_pls__elementwiseProduct mat.VecDense
var __f_xc_is_n_pls__f_xc_is_n_pls mat.VecDense

/*
係数 f_XC を求める。

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
	f_wsv_js_n_pls *mat.VecDense,
	f_xot_is_is_n_pls mat.Matrix,
	k_r_is_n *mat.VecDense,
) *mat.VecDense {
	// f_wsc_js_n_pls + f_wsv_js_n_pls
	f_wscPlusWsv := &__f_xc_is_n_pls__f_wscPlusWsv
	f_wscPlusWsv.AddVec(f_wsc_js_n_pls, f_wsv_js_n_pls)

	// np.dot(f_mrt_hum_is_js, (f_wsc_js_n_pls + f_wsv_js_n_pls))
	dotProduct := &__f_xc_is_n_pls__dotProduct
	dotProduct.MulVec(f_mrt_hum_is_js, f_wscPlusWsv)

	// k_r_is_n * np.dot(f_mrt_hum_is_js, (f_wsc_js_n_pls + f_wsv_js_n_pls))
	elementwiseProduct := &__f_xc_is_n_pls__elementwiseProduct
	elementwiseProduct.MulElemVec(k_r_is_n, dotProduct)

	// np.dot(f_xot_is_is_n_pls, k_r_is_n * np.dot(f_mrt_hum_is_js, (f_wsc_js_n_pls + f_wsv_js_n_pls)))
	f_xc_is_n_pls := &__f_xc_is_n_pls__f_xc_is_n_pls
	f_xc_is_n_pls.MulVec(f_xot_is_is_n_pls, elementwiseProduct)

	return f_xc_is_n_pls
}

/*
係数 f_XOT を求める。

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
	k_c_is_n *mat.VecDense,
	k_r_is_n *mat.VecDense,
) *mat.Dense {
	// Python: np.dot(f_mrt_hum_is_js, f_wsr_js_is)
	temp2 := &mat.Dense{}
	temp2.Mul(f_mrt_hum_is_js, f_wsr_js_is)

	// temp2 := k_r_is_n * temp2
	__ScaleColumn(k_r_is_n, temp2)

	// temp2 := diag(k_c_is_n) + temp2
	AddAsDiag(temp2, k_c_is_n)

	// result = temp2^-1
	result := &mat.Dense{}
	err := result.Inverse(temp2)
	if err != nil {
		panic(err)
	}

	return result
}

var __k_c_is_n__sum mat.VecDense
var __k_c_is_n__result mat.VecDense

/*
人体表面の対流熱伝達率が総合熱伝達率に占める割合を求める。

Args:

	h_hum_c_is_n: ステップ n における室 i の人体表面の対流熱伝達率, W/(m2 K)
	h_hum_r_is_n: ステップ n における室 i の人体表面の放射熱伝達率, W/(m2 K)

Returns:

	ステップ n における室 i の人体表面の対流熱伝達率が総合熱伝達率に占める割合, -, [i, 1]

Notes:

	式(2.21)
*/
func get_k_c_is_n(h_hum_c_is_n *mat.VecDense, h_hum_r_is_n *mat.VecDense) *mat.VecDense {
	// h_hum_c_is_n + h_hum_r_is_n
	sum := &__k_c_is_n__sum
	sum.AddVec(h_hum_c_is_n, h_hum_r_is_n)

	// h_hum_c_is_n / (h_hum_c_is_n + h_hum_r_is_n)
	result := &__k_c_is_n__result
	result.DivElemVec(h_hum_c_is_n, sum)

	return result
}

var __k_r_is_n__sum mat.VecDense
var __k_r_is_n__result mat.VecDense

/*
人体表面の放射熱伝達率が総合熱伝達率に占める割合を求める。

Args:

	h_hum_c_is_n: ステップ n における室 i の人体表面の対流熱伝達率, W/(m2 K)
	h_hum_r_is_n: ステップ n における室 i の人体表面の放射熱伝達率, W/(m2 K)

Returns:

	ステップ n における室 i の人体表面の放射熱伝達率が総合熱伝達率に占める割合, -, [i, 1]

Notes:

	式(2.22)
*/
func get_k_r_is_n(h_hum_c_is_n *mat.VecDense, h_hum_r_is_n *mat.VecDense) *mat.VecDense {
	// h_hum_c_is_n + h_hum_r_is_n
	sum := &__k_r_is_n__sum
	sum.AddVec(h_hum_c_is_n, h_hum_r_is_n)

	// h_hum_r_is_n / (h_hum_c_is_n + h_hum_r_is_n)
	result := &__k_r_is_n__result
	result.DivElemVec(h_hum_r_is_n, sum)

	return result
}

var __f_brm_is_is_n_pls__temp1 mat.VecDense
var __f_brm_is_is_n_pls__temp2 mat.Dense
var __f_brm_is_is_n_pls__temp3 mat.Dense
var __f_brm_is_is_n_pls__temp4 mat.VecDense
var __f_brm_is_is_n_pls__temp5 mat.VecDense
var __f_brm_is_is_n_pls__temp6 *mat.Dense
var __f_brm_is_is_n_pls__temp7 mat.Dense
var __f_brm_is_is_n_pls__temp8 mat.Dense
var __f_brm_is_is_n_pls__temp9 mat.Dense
var __f_brm_is_is_n_pls__temp10 *mat.Dense

// var __f_brm_is_is_n_pls__result1 mat.Dense
// var __f_brm_is_is_n_pls__result2 mat.Dense

/*
係数 f^_BRMを求める。

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
	a_s_js *mat.VecDense,
	v_rm_is *mat.VecDense,
	c_sh_frt_is *mat.VecDense,
	delta_t float64,
	f_wsr_js_is mat.Matrix,
	g_sh_frt_is *mat.VecDense,
	h_s_c_js *mat.VecDense,
	p_is_js mat.Matrix,
	p_js_is mat.Matrix,
	v_vent_int_is_is_n mat.Matrix,
	v_vent_out_is_n *mat.VecDense,
	v_vent_ntr_set_is []float64,
) (*mat.Dense, *mat.Dense) {
	// Python: v_rm_is * rho_a * c_a / delta_t
	temp1 := &__f_brm_is_is_n_pls__temp1
	temp1.ScaleVec(rho_a*c_a/delta_t, v_rm_is)

	// Python: (p_js_is - f_wsr_js_is) * a_s_js * h_s_c_js
	temp2 := &__f_brm_is_is_n_pls__temp2
	temp2.Sub(p_js_is, f_wsr_js_is) // (p_js_is - f_wsr_js_is)
	var temp2_1 mat.VecDense
	temp2_1.MulElemVec(a_s_js, h_s_c_js) // a_s_js * h_s_c_js
	__ScaleRows(temp2, &temp2_1)

	// Python: np.dot(p_is_js, temp2)
	temp3 := &__f_brm_is_is_n_pls__temp3
	temp3.Mul(p_is_js, temp2)

	// Python: c_sh_frt_is * g_sh_frt_is / (c_sh_frt_is + g_sh_frt_is * delta_t)
	temp4 := &__f_brm_is_is_n_pls__temp4
	temp4.MulElemVec(c_sh_frt_is, g_sh_frt_is)

	temp5 := &__f_brm_is_is_n_pls__temp5
	temp5.AddScaledVec(c_sh_frt_is, delta_t, g_sh_frt_is)
	temp4.DivElemVec(temp4, temp5)

	// Python: c_a * rho_a * (v_diag(v_vent_out_is_n) - v_vent_int_is_is_n)
	if __f_brm_is_is_n_pls__temp6 == nil {
		n := v_vent_out_is_n.Len()
		__f_brm_is_is_n_pls__temp6 = mat.NewDense(n, n, nil)
	}
	temp6 := __f_brm_is_is_n_pls__temp6
	temp7 := &__f_brm_is_is_n_pls__temp7

	Diag(temp6, v_vent_out_is_n)         // temp6 <= diag(v_vent_out_is_n)
	temp7.Sub(temp6, v_vent_int_is_is_n) // temp7 <= temp6 - v_vent_int_is_is_n
	temp7.Scale(c_a*rho_a, temp7)        // temp7 <= c_a * rho_a * temp7

	// Python: temp3 := diag(temp1) + temp3 + diag(temp4) + temp8
	AddAsDiag(temp3, temp1)
	AddAsDiag(temp3, temp4)
	temp3.Add(temp3, temp7)

	temp9 := &__f_brm_is_is_n_pls__temp9
	//temp8 := &__f_brm_is_is_n_pls__temp8
	if __f_brm_is_is_n_pls__temp10 == nil {
		n := len(v_vent_ntr_set_is)
		__f_brm_is_is_n_pls__temp10 = mat.NewDense(n, n, nil)
	}
	temp10 := __f_brm_is_is_n_pls__temp10
	DiagFromFloat64(temp10, v_vent_ntr_set_is) // temp10 <= diag
	temp10.Scale(c_a*rho_a, temp10)
	temp9.Add(temp3, temp10)

	return temp3, temp9
}

func Diag(dst *mat.Dense, src *mat.VecDense) {
	n := src.Len()
	for i := 0; i < n; i++ {
		dst.Set(i, i, src.AtVec(i))
	}
}

func DiagFromFloat64(dst *mat.Dense, src []float64) {
	n := len(src)
	for i := 0; i < n; i++ {
		dst.Set(i, i, src[i])
	}
}

func NewDiagAsDenseFromFloat64(src []float64) *mat.Dense {
	n := len(src)
	mat := mat.NewDense(n, n, nil)
	for i := 0; i < n; i++ {
		mat.Set(i, i, src[i])
	}
	return mat
}

// 密行列に対角化行列を加算する
// NOTE: AddはDense+Dense以外はSet/Atを使うのがパフォーマンスのペナルティが大きいため、専用関数を作った
func AddAsDiag(dst *mat.Dense, src *mat.VecDense) {
	n := src.Len()
	for i := 0; i < n; i++ {
		dst.Set(i, i, dst.At(i, i)+src.AtVec(i))
	}
}

// var __f_brc_is_n_pls__result1 mat.VecDense
var __f_brc_is_n_pls__tmpVec1 mat.VecDense
var __f_brc_is_n_pls__tmpMat1 mat.VecDense
var __f_brc_is_n_pls__tmpVec4 mat.VecDense
var __f_brc_is_n_pls__tmpVec6 mat.VecDense

/*
係数 f^_BRC を求める。

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
	a_s_js *mat.VecDense,
	v_rm_is *mat.VecDense,
	c_sh_frt_is *mat.VecDense,
	delta_t float64,
	f_wsc_js_n_pls mat.Vector,
	f_wsv_js_n_pls *mat.VecDense,
	g_sh_frt_is *mat.VecDense,
	h_s_c_js *mat.VecDense,
	p_is_js mat.Matrix,
	q_gen_is_n mat.Vector,
	q_hum_is_n mat.Vector,
	q_sol_frt_is_n mat.Vector,
	theta_frt_is_n *mat.VecDense,
	theta_o_n_pls float64,
	theta_r_is_n *mat.VecDense,
	v_vent_out_non_nv_is_n *mat.VecDense,
	v_vent_ntr_is_n []float64,
) (*mat.VecDense, *mat.VecDense) {
	// v_rm_is * c_a * rho_a / delta_t * theta_r_is_n
	var result1 mat.VecDense
	result1.ScaleVec(c_a*rho_a/delta_t, v_rm_is) // v_rm_is * c_a * rho_a / delta_t
	result1.MulElemVec(&result1, theta_r_is_n)   // ~~~ * theta_r_is_n

	// (h_s_c_js * a_s_js * (f_wsc_js_n_pls + f_wsv_js_n_pls)
	tmpVec1 := &__f_brc_is_n_pls__tmpVec1
	tmpVec1.AddVec(f_wsc_js_n_pls, f_wsv_js_n_pls)
	tmpVec1.MulElemVec(tmpVec1, h_s_c_js)
	tmpVec1.MulElemVec(tmpVec1, a_s_js)

	// + np.dot(p_is_js, h_s_c_js * a_s_js * (f_wsc_js_n_pls + f_wsv_js_n_pls))
	tmpMat1 := &__f_brc_is_n_pls__tmpMat1
	tmpMat1.MulVec(p_is_js, tmpVec1)
	result1.AddVec(&result1, tmpMat1)

	// + c_a * rho_a * v_vent_out_non_nv_is_n * theta_o_n_pls
	result1.AddScaledVec(&result1, c_a*rho_a*theta_o_n_pls, v_vent_out_non_nv_is_n)

	// + q_gen_is_n + q_hum_is_n
	result1.AddVec(&result1, q_gen_is_n)
	result1.AddVec(&result1, q_hum_is_n)

	// c_sh_frt_is * theta_frt_is_n
	tmpVec4 := &__f_brc_is_n_pls__tmpVec4
	tmpVec4.MulElemVec(c_sh_frt_is, theta_frt_is_n)

	// NOTE: 計算仕様では q_sol_frt_is_n を加算しない
	// c_sh_frt_is * theta_frt_is_n + q_sol_frt_is_n * delta_t
	tmpVec4.AddScaledVec(tmpVec4, delta_t, q_sol_frt_is_n)

	// c_sh_frt_is + delta_t * g_sh_frt_is
	tmpVec6 := &__f_brc_is_n_pls__tmpVec6
	tmpVec6.AddScaledVec(c_sh_frt_is, delta_t, g_sh_frt_is)

	// g_sh_frt_is * (c_sh_frt_is * theta_frt_is_n + q_sol_frt_is_n * delta_t) / (c_sh_frt_is + delta_t * g_sh_frt_is)
	tmpVec4.DivElemVec(tmpVec4, tmpVec6)
	tmpVec4.MulElemVec(g_sh_frt_is, tmpVec4)

	// + g_sh_frt_is * (c_sh_frt_is * theta_frt_is_n + q_sol_frt_is_n * delta_t) / (c_sh_frt_is + delta_t * g_sh_frt_is)
	result1.AddVec(&result1, tmpVec4)

	//  result1 + c_a * rho_a * v_vent_ntr_is_n * theta_o_n_pls
	var result2 mat.VecDense
	result2.AddScaledVec(&result1, c_a*rho_a*theta_o_n_pls, mat.NewVecDense(len(v_vent_ntr_is_n), v_vent_ntr_is_n))

	return &result1, &result2
}

var __v_vent_out_non_ntr_is_n mat.VecDense

/*
換気・すきま風・自然風の利用による外気の流入量を求める。

Args:

	v_leak_is_n: ステップ n からステップ n+1 における室 i のすきま風量, m3/s, [i, 1]
	v_vent_mec_is_n: ステップ n からステップ n+1 における室 i の機械換気量（全般換気量と局所換気量の合計値）, m3/s, [i, 1]

Returns:

	ステップ n からステップ n+1 における室 i の換気・すきま風・自然風の利用による外気の流入量, m3/s

Notes:

	式(2.25)
*/
func get_v_vent_out_non_ntr_is_n(
	v_leak_is_n *mat.VecDense,
	v_vent_mec_is_n mat.Vector,
) *mat.VecDense {
	__v_vent_out_non_ntr_is_n.AddVec(v_leak_is_n, v_vent_mec_is_n)
	return &__v_vent_out_non_ntr_is_n
}

/*
自然風利用による換気量を求める。

Args:

	operation_mode_is_n: ステップ n からステップ n+1 における室 i の運転モード, [i, 1]
	v_vent_ntr_set_is: 室 i の自然風利用時の換気量, m3/s

Returns:

	ステップ n からステップ n+1 における室 i の自然風利用による換気量, m3/s, [i, 1]

Notes:

	式(2.26)
*/
func get_v_vent_ntr_is_n(operation_mode_is_n []OperationMode, v_vent_ntr_set_is *mat.VecDense) *mat.VecDense {
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

var __f_wsv_js_n_pls mat.VecDense

/*
係数 f_WSV を求める。

Args:

	f_cvl_js_n_pls: ステップ n+1 における係数 f_CVL, degree C, [j, 1]
	f_ax_js_js: 係数 f_AX, -, [j, j]

Returns:

	ステップ n+1 の係数 f_WSV, degree C, [j, 1]

Notes:

	式(2.27)
*/
func get_f_wsv_js_n_pls(
	f_cvl_js_n_pls *mat.VecDense,
	f_ax_js_js_inv *mat.Dense,
	f_ax_js_js *mat.LU,
) *mat.VecDense {
	// 総じて逆行列を用いた場合のほうが速い
	__f_wsv_js_n_pls.MulVec(f_ax_js_js_inv, f_cvl_js_n_pls)

	//f_ax_js_js.SolveVecTo(&__f_wsv_js_n_pls, false, f_cvl_js_n_pls)

	return &__f_wsv_js_n_pls
}

var __f_cvl_js_n_pls *mat.VecDense

/*
係数 f_CVL を求める。

Args:

	theta_dsh_s_a_js_ms_n_pls: ステップ n+1 における境界 j の項別公比法の指数項 m の吸熱応答の項別成分, degree C, [j, m]
	theta_dsh_s_t_js_ms_n_pls: ステップ n+1 における境界 j の項別公比法の指数項 m の貫流応答の項別成分, degree C, [j, m]

Returns:

	ステップ n+1 における係数 f_CVL, degree C, [j, 1]

Notes:

	式(2.28)
*/
func get_f_cvl_js_n_pls(
	theta_dsh_s_a_js_ms_n_pls *mat.Dense,
	theta_dsh_s_t_js_ms_n_pls *mat.Dense,
) *mat.VecDense {
	// Pythonコード: np.sum(theta_dsh_s_t_js_ms_n_pls + theta_dsh_s_a_js_ms_n_pls, axis=1, keepdims=True)
	r, _ := theta_dsh_s_t_js_ms_n_pls.Dims()
	if __f_cvl_js_n_pls == nil {
		__f_cvl_js_n_pls = mat.NewVecDense(r, nil)
	}
	result := __f_cvl_js_n_pls
	for i := 0; i < r; i++ {
		rowSum := floats.Sum(theta_dsh_s_t_js_ms_n_pls.RawRowView(i)) + floats.Sum(theta_dsh_s_a_js_ms_n_pls.RawRowView(i))
		result.SetVec(i, rowSum)
	}

	return result
}

var __theta_dsh_s_a_js_ms_n_pls *mat.Dense
var __theta_dsh_s_a_js_ms_n_pls__temp mat.Dense

/*
項別公比法の吸熱応答の項別成分を求める。

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
	phi_a1_js_ms *mat.Dense,
	q_s_js_n *mat.VecDense,
	r_js_ms mat.Matrix,
	theta_dsh_srf_a_js_ms_n mat.Matrix,
) *mat.Dense {
	if __theta_dsh_s_a_js_ms_n_pls == nil {
		r, c := phi_a1_js_ms.Dims()
		__theta_dsh_s_a_js_ms_n_pls = mat.NewDense(r, c, nil)
	}

	// result = phi_a1_js_ms * q_s_js_n
	result := __theta_dsh_s_a_js_ms_n_pls
	__ScaleRowsTo(result, q_s_js_n, phi_a1_js_ms)

	// temp = r_js_ms * theta_dsh_srf_a_js_ms_n
	temp := &__theta_dsh_s_a_js_ms_n_pls__temp
	temp.MulElem(r_js_ms, theta_dsh_srf_a_js_ms_n)

	// result = result + temp
	result.Add(result, temp)

	return result
}

var __theta_dsh_s_t_js_ms_n_pls *mat.Dense
var __theta_dsh_s_t_js_ms_n_pls__temp mat.Dense

/*
項別公比法の貫流応答の項別成分を求める。

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
	phi_t1_js_ms *mat.Dense,
	r_js_ms *mat.Dense,
	theta_dsh_srf_t_js_ms_n *mat.Dense,
	theta_rear_js_n *mat.VecDense,
) *mat.Dense {
	if __theta_dsh_s_t_js_ms_n_pls == nil {
		r, c := phi_t1_js_ms.Dims()
		__theta_dsh_s_t_js_ms_n_pls = mat.NewDense(r, c, nil)
	}

	result := __theta_dsh_s_t_js_ms_n_pls

	// result = theta_rear_js_n * phi_t1_js_ms
	__ScaleRowsTo(result, theta_rear_js_n, phi_t1_js_ms)

	// temp = r_js_ms * theta_dsh_srf_t_js_ms_n
	temp := &__theta_dsh_s_t_js_ms_n_pls__temp
	temp.MulElem(r_js_ms, theta_dsh_srf_t_js_ms_n)

	// result = result + temp
	result.Add(result, temp)

	return result
}

/*
人体発熱を求める。

Args:

	n_hum_is_n: ステップ n からステップ n+1 における室 i の在室人数, -, [i, 1]
	q_hum_psn_is_n: ステップ n からステップ n+1 における室 i の1人あたりの人体発熱, W, [i, 1]

Returns:

	ステップ n からステップ n+1 における室 i の人体発熱, W, [i, 1]

Notes:

	式(2.31)
*/
func get_q_hum_is_n(n_hum_is_n mat.Vector, q_hum_psn_is_n *mat.VecDense) *mat.VecDense {
	var result mat.VecDense
	result.MulElemVec(n_hum_is_n, q_hum_psn_is_n)
	return &result
}

var __theta_s_rear_js_n__result1 mat.VecDense
var __theta_s_rear_js_n__result2 mat.VecDense
var __theta_s_rear_js_n__result3 mat.VecDense

/*
裏面温度を求める。

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
	k_s_er_js_js *mat.Dense,
	theta_er_js_n *mat.VecDense,
	k_s_eo_js *mat.VecDense,
	theta_eo_js_n mat.Vector,
	k_s_r_js_is *mat.Dense,
	theta_r_is_n *mat.VecDense,
) *mat.VecDense {

	result1 := &__theta_s_rear_js_n__result1
	result2 := &__theta_s_rear_js_n__result2
	result3 := &__theta_s_rear_js_n__result3

	// np.dot(k_s_er_js_js, theta_er_js_n)
	result1.MulVec(k_s_er_js_js, theta_er_js_n)

	// k_s_eo_js*theta_eo_js_n
	result2.MulElemVec(k_s_eo_js, theta_eo_js_n)

	// np.dot(k_s_r_js_is, theta_r_is_n)
	result3.MulVec(k_s_r_js_is, theta_r_is_n)

	//np.dot(k_s_er_js_js, theta_er_js_n) + k_s_eo_js * theta_eo_js_n + np.dot(k_s_r_js_is, theta_r_is_n)
	result1.AddVec(result1, result2) // result1 <= result1 + result2
	result3.AddVec(result3, result1) // result3 <= result3 + result1

	return result3
}

// ----------------------------------------------------------------------------------
// 5 事前計算
// ----------------------------------------------------------------------------------

/*
係数 f_WSC を求める。

Args:

	f_ax_js_js_revert: 係数 f_{AX}^-1, -, [j, j]
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
係数 f_WSR を求める。

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
係数 f_CRX を求める。

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
	h_s_c_js *mat.VecDense,
	h_s_r_js *mat.VecDense,
	k_ei_js_js mat.Matrix,
	phi_a0_js *mat.VecDense,
	phi_t0_js *mat.VecDense,
	q_s_sol_js_ns mat.Matrix,
	k_eo_js *mat.VecDense,
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
係数 f_FIA を求める。

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
	h_s_c_js *mat.VecDense,
	h_s_r_js *mat.VecDense,
	k_ei_js_js mat.Matrix,
	p_js_is mat.Matrix,
	phi_a0_js *mat.VecDense,
	phi_t0_js *mat.VecDense,
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
係数 f_AX を求める。

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
	h_s_c_js *mat.VecDense,
	h_s_r_js *mat.VecDense,
	k_ei_js_js mat.Matrix,
	p_js_is mat.Matrix,
	phi_a0_js *mat.VecDense,
	phi_t0_js *mat.VecDense,
) (*mat.LU, *mat.Dense) {
	// 1.0 + phi_a0_js * (h_s_r_js + p_js_is)
	temp1 := make([]float64, phi_a0_js.Len())
	for i := 0; i < phi_a0_js.Len(); i++ {
		temp1[i] = 1.0 + phi_a0_js.AtVec(i)*(h_s_c_js.AtVec(i)+h_s_r_js.AtVec(i))
	}

	// diag of temp1
	temp2 := NewDiagAsDenseFromFloat64(temp1)

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

	// 高速化のために、ここでLU分解と逆行列の計算を行う
	var lu mat.LU
	lu.Factorize(&result)

	var invert mat.Dense
	invert.Inverse(&result)

	return &lu, &invert
}

/*
機械換気量（全般換気量と局所換気量の合計値）を求める。

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
	// 式(在室者の形態係数:1b)
	// 室 i の在室者に対する境界jの形態係数, [i, j]
	f_mrt_hum_is_js := get_f_mrt_hum_js(bs.p_is_js, bs.a_s_js, bs.is_floor_js)

	// 式(室内の境界の形態係数および放射熱伝達率:1)
	// 室 i の微小球に対する境界 j の重み係数, -, [i, j]
	f_mrt_is_js := get_f_mrt_is_js(bs.a_s_js, bs.h_s_r_js, bs.p_is_js)

	// 式(4.7)
	// ステップ n からステップ n+1 における室 i の機械換気量（全般換気量と局所換気量の合計値）, m3/s, [i, 1]
	v_vent_mec_is_ns := get_v_vent_mec_is_ns(mvs.v_vent_mec_general_is, scd.v_mec_vent_local_is_ns)

	// 式(日射吸収量:1)
	// ステップ n における室 i に設置された備品等による透過日射吸収熱量, W, [i, n+1]
	q_sol_frt_is_ns := get_q_sol_frt_is_ns(bs.q_trs_sol_is_ns)

	// 式(日射吸収量:2)
	// ステップ n における境界 j の透過日射吸収熱量, W/m2, [j, n]
	q_s_sol_js_ns := get_q_s_sol_js_ns(
		bs.p_is_js,
		bs.a_s_js,
		bs.p_s_sol_abs_js,
		bs.p_js_is,
		bs.q_trs_sol_is_ns,
	)

	// 式(4.5)
	// 係数 f_AX, -, [j, j]
	f_ax_js_js, f_ax_js_js_inv := get_f_ax_js_is(
		f_mrt_is_js,
		bs.h_s_c_js,
		bs.h_s_r_js,
		bs.k_ei_js_js,
		bs.p_js_is,
		bs.phi_a0_js,
		bs.phi_t0_js,
	)

	// 式(4.4)
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

	// 式(4.3)
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

	// 式(4.1)
	// 係数 f_WSR, -, [j, i]
	f_wsr_js_is := get_f_wsr_js_is(f_ax_js_js, f_fia_js_is)

	// 式(4.2)
	// 係数 f_{WSC, n}, degree C, [j, n]
	f_wsc_js_ns := get_f_wsc_js_ns(f_ax_js_js, f_crx_js_ns)

	// 式(2.21) (2.22)
	// ステップnにおける室iの在室者表面における対流熱伝達率の総合熱伝達率に対する比, -, [i, 1]
	// ステップ n における室 i の在室者表面における放射熱伝達率の総合熱伝達率に対する比, -, [i, 1]
	// NOTE: 計算仕様では毎時計算になっているが、ここでは事前計算している
	k_c_is_n, k_r_is_n := op.get_k_is()

	// 式(2.20)
	// ステップn+1における室iの係数 XOT, [i, i]
	// NOTE: 計算仕様では毎時計算扱いになっているが、ここでは事前計算している
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
		f_ax_js_js_inv:    f_ax_js_js_inv,
		f_wsc_js_ns:       f_wsc_js_ns,
		k_r_is_n:          k_r_is_n,
		k_c_is_n:          k_c_is_n,
		f_xot_is_is_n_pls: f_xot_is_is_n_pls,
	}

	return &pre_calc_parameters
}
