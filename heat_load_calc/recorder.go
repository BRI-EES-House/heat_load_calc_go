package heat_load_calc

type Recorder struct {
	YEAR                 int
	_itv                 Interval
	_id_rm_is            []int
	_id_bdry_js          []int
	_n_step_i            int
	_n_step_a            int
	theta_i_ns           [][]float64
	theta_o_ns           [][]float64
	x_o_ns               [][]float64
	theta_r_is_ns        [][]float64
	rh_r_is_ns           [][]float64
	x_r_is_ns            [][]float64
	theta_mrt_hum_is_ns  [][]float64
	theta_ot             []float64
	q_trs_sol_is_ns      [][]float64
	theta_frt_is_ns      [][]float64
	q_sol_frt_is_ns      [][]float64
	x_frt_is_ns          [][]float64
	pmv_is_ns            [][]float64
	ppd_is_ns            [][]float64
	theta_s_js_ns        [][]float64
	theta_ei_js_ns       [][]float64
	theta_rear_js_ns     [][]float64
	h_s_r_js_ns          [][]float64
	q_r_js_ns            [][]float64
	h_s_c_js_ns          [][]float64
	q_c_js_ns            [][]float64
	q_i_sol_s_ns_js      [][]float64
	q_s_js_ns            [][]float64
	operation_mode_is_ns [][]float64
	ac_demand_is_ns      [][]float64
	h_hum_c_is_ns        [][]float64
	h_hum_r_is_ns        [][]float64
	q_gen_is_ns          [][]float64
	x_gen_is_ns          [][]float64
	q_hum_is_ns          [][]float64
	x_hum_is_ns          [][]float64
	l_cs_is_ns           [][]float64
	l_rs_is_ns           [][]float64
	l_cl_is_ns           [][]float64
	q_frt_is_ns          [][]float64
	q_l_frt_is_ns        [][]float64
	v_reak_is_ns         [][]float64
	v_ntrl_is_ns         [][]float64
	v_hum_is_ns          [][]float64
	clo_is_ns            [][]float64
	_output_list_room_a  []struct {
		field1 string
		field2 string
	}
	_output_list_room_i []struct {
		field1 string
		field2 string
	}
	_output_list_boudary_i []struct {
		field1 string
		field2 string
	}
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

	r.theta_o_ns = make([][]float64, r._n_step_i)
	r.x_o_ns = make([][]float64, r._n_step_i)
	r.theta_r_is_ns = make([][]float64, n_rm)
	r.rh_r_is_ns = make([][]float64, n_rm)
	r.x_r_is_ns = make([][]float64, n_rm)
	r.theta_mrt_hum_is_ns = make([][]float64, n_rm)
	r.theta_ot = make([]float64, r._n_step_i)
	r.q_trs_sol_is_ns = make([][]float64, n_rm)
	r.theta_frt_is_ns = make([][]float64, n_rm)
	r.q_sol_frt_is_ns = make([][]float64, n_rm)
	r.x_frt_is_ns = make([][]float64, n_rm)
	r.pmv_is_ns = make([][]float64, n_rm)
	r.ppd_is_ns = make([][]float64, n_rm)

	r.theta_s_js_ns = make([][]float64, n_boundries)
	r.theta_ei_js_ns = make([][]float64, n_boundries)
	r.theta_rear_js_ns = make([][]float64, n_boundries)
	r.h_s_r_js_ns = make([][]float64, n_boundries)
	r.q_r_js_ns = make([][]float64, n_boundries)
	r.h_s_c_js_ns = make([][]float64, n_boundries)
	r.q_c_js_ns = make([][]float64, n_boundries)
	r.q_i_sol_s_ns_js = make([][]float64, n_boundries)
	r.q_s_js_ns = make([][]float64, n_boundries)

	r.operation_mode_is_ns = make([][]float64, n_rm)
	r.ac_demand_is_ns = make([][]float64, n_rm)
	r.h_hum_c_is_ns = make([][]float64, n_rm)
	r.h_hum_r_is_ns = make([][]float64, n_rm)
	r.q_gen_is_ns = make([][]float64, n_rm)
	r.x_gen_is_ns = make([][]float64, n_rm)
	r.q_hum_is_ns = make([][]float64, n_rm)
	r.x_hum_is_ns = make([][]float64, n_rm)
	r.l_cs_is_ns = make([][]float64, n_rm)
	r.l_cl_is_ns = make([][]float64, n_rm)
	r.q_frt_is_ns = make([][]float64, n_rm)
	r.q_l_frt_is_ns = make([][]float64, n_rm)
	r.v_reak_is_ns = make([][]float64, n_rm)
	r.v_ntrl_is_ns = make([][]float64, n_rm)
	r.v_hum_is_ns = make([][]float64, n_rm)
	r.clo_is_ns = make([][]float64, n_rm)

	// r._output_list_room_a = []struct{ column_name, field_name string }{
	// 	{"operation_mode_is_ns", "ac_operate"},
	// 	{"ac_demand_is_ns", "occupancy"},
	// 	{"h_hum_c_is_ns", "hc_hum"},
	// 	{"h_hum_r_is_ns", "hr_hum"},
	// 	{"q_gen_is_ns", "q_s_except_hum"},
	// 	{"x_gen_is_ns", "q_l_except_hum"},
	// 	{"q_hum_is_ns", "q_hum_s"},
	// 	{"x_hum_is_ns", "q_hum_l"},
	// 	{"l_cs_is_ns", "l_s_c"},
	// 	{"l_rs_is_ns", "l_s_r"},
	// 	{"l_cl_is_ns", "l_l_c"},
	// 	{"q_frt_is_ns", "q_s_fun"},
	// 	{"q_l_frt_is_ns", "q_l_fun"},
	// 	{"v_reak_is_ns", "v_reak"},
	// 	{"v_ntrl_is_ns", "v_ntrl"},
	// 	{"v_hum_is_ns", "v_hum"},
	// 	{"clo_is_ns", "clo"},
	// }
	// r._output_list_room_i = []struct{ column_name, field_name string }{
	// 	{"theta_r_is_ns", "t_r"},
	// 	{"rh_r_is_ns", "rh_r"},
	// 	{"x_r_is_ns", "x_r"},
	// 	{"theta_mrt_hum_is_ns", "mrt"},
	// 	{"theta_ot", "ot"},
	// 	{"q_trs_sol_is_ns", "q_sol_t"},
	// 	{"theta_frt_is_ns", "t_fun"},
	// 	{"q_sol_frt_is_ns", "q_s_sol_fun"},
	// 	{"x_frt_is_ns", "x_fun"},
	// 	{"pmv_is_ns", "pmv"},
	// 	{"ppd_is_ns", "ppd"},
	// }

	// r._output_list_boudary_i = []struct{ column_name, field_name string }{
	// 	{"theta_s_js_ns", "t_s"},
	// 	{"theta_ei_js_ns", "t_e"},
	// 	{"theta_rear_js_ns", "t_b"},
	// 	{"h_s_r_js_ns", "hir_s"},
	// 	{"q_r_js_ns", "qir_s"},
	// 	{"h_s_c_js_ns", "hic_s"},
	// 	{"q_c_js_ns", "qic_s"},
	// 	{"q_i_sol_s_ns_js", "qisol_s"},
	// 	{"q_s_js_ns", "qiall_s"},
	// }

	return &r
}

func (r *Recorder) pre_recording(ss *PreCalcParameters) {
	// 注意：用意された1年分のデータと実行期間が異なる場合があるためデータスライスする必要がある。

	// ---瞬時値---

	// // ステップ n における外気温度, ℃, [n+1]
	// r.theta_o_ns = ss.theta_o_ns[0:r._n_step_i]

	// // ステップ n における外気絶対湿度, kg/kg(DA), [n+1]
	// r.x_o_ns = ss.x_o_ns[0:r._n_step_i]

	// r.q_trs_sol_is_ns = make([][]float64, len(r._id_rm_is))
	// r.q_sol_frt_is_ns = make([][]float64, len(r._id_rm_is))
	// for i, id_rm := range r._id_rm_is {
	// 	// ステップ n における室 i の窓の透過日射熱取得, W, [i, n+1]
	// 	r.q_trs_sol_is_ns[i] = ss.q_trs_sol_is_ns[i][0:r._n_step_i]

	// 	// ステップ n における室 i に設置された備品等による透過日射吸収熱量, W, [i, n+1]
	// 	r.q_sol_frt_is_ns = ss.q_sol_frt_is_ns[i][0:r._n_step_i]
	// }

	// r.q_i_sol_s_ns_js = make([][]float64, len(r._id_bdry_js))
	// for i, id_js := range r._id_bdry_js {
	// 	// ステップ n の境界 j の表面日射熱流, W, [j, n+1]
	// 	r.q_i_sol_s_ns_js = ss.q_s_sol_js_ns[i][0:r._n_step_i] * ss.a_s_js[i]
	// }

	// r.h_s_c_js_ns = make([][]float64, len(r._id_bdry_js))
	// r.h_s_r_js_ns = make([][]float64, len(r._id_bdry_js))
	// for i, id_js := range r._id_bdry_js {
	// 	r.h_s_c_js_ns[i] = make([]float64, r._n_step_i)
	// 	r.h_s_r_js_ns[i] = make([]float64, r._n_step_i)

	// 	for n := 0; n < r._n_step_i; n++ {
	// 		// ステップ n の境界 j の表面対流熱伝達率, W/m2K, [j, n+1]
	// 		r.h_s_c_js_ns[i][n] = ss.h_s_c_js[n]

	// 		// ステップ n の境界 j の表面放射熱伝達率, W/m2K, [j, n+1]
	// 		r.h_s_r_js_ns[i][n] = ss.h_s_r_js[n]
	// 	}
	// }

	// // ---平均値・積算値---

	// r.ac_demand_is_ns = make([][]float64, len(r._id_rm_is))
	// r.q_gen_is_ns = make([][]float64, len(r._id_rm_is))
	// r.x_gen_is_ns = make([][]float64, len(r._id_rm_is))
	// for i, id_rm := range r._id_rm_is {
	// 	// ステップ n の室 i における当該時刻の空調需要, [i, n]
	// 	r.ac_demand_is_ns[i] = ss.ac_demand_is_ns[i][0:r._n_step_a]

	// 	// ステップnの室iにおける人体発熱を除く内部発熱, W, [i, 8760*4]
	// 	r.q_gen_is_ns[i] = ss.q_gen_is_ns[i][0:r._n_step_a]

	// 	// ステップ n の室 i における人体発湿を除く内部発湿, kg/s, [i, n]
	// 	r.x_gen_is_ns[i] = ss.x_gen_is_ns[i][0:r._n_step_a]
	// }
}

func (r *Recorder) post_recording(ss *PreCalcParameters) {
	// // ---瞬時値---

	// // ステップ n の室 i における飽和水蒸気圧, Pa, [i, n+1]
	// p_vs_is_ns := mat.DenseCopyOf(r.theta_r_is_ns)
	// p_vs_is_ns.Apply(func(_, _ int, v float64) float64 { return get_p_vs(v) })

	// // ステップ n における室 i の水蒸気圧, Pa, [i, n+1]
	// p_v_is_ns := mat.DenseCopyOf(r.x_r_is_ns)
	// p_v_is_ns.Apply(func(_, _ int, v float64) float64 { return get_p_v_r_is_n(v) })

	// // ステップnの室iにおける相対湿度, %, [i, n+1]
	// r.rh_r_is_ns = mat.DenseCopyOf(r.theta_r_is_ns)
	// r.rh_r_is_ns.Apply(func(i, j int, _ float64) float64 {
	// 	return get_h(p_v_is_ns.At(i, j, p_vs_is_ns.At(i, j)))
	// })

	// // ステップnの境界jにおける表面熱流（壁体吸熱を正とする）のうち放射成分, W, [j, n]
	// var temp1, temp2 mat.Dense
	// temp1.Mul(ss.p_js_is, ss.f_mrt_is_js)
	// temp1.Mul(&temp, r.theta_s_js_ns)
	// temp1.Sub(&temp, r.theta_s_js_ns)
	// temp2.MulElem(ss.h_s_r_js, ss.a_s_js)
	// temp2.MulElem(&temp2, &temp1)
	// r.q_r_js_ns = &temp2
	// //ss.h_s_r_js * ss.a_s_js * (np.dot(np.dot(ss.p_js_is, ss.f_mrt_is_js), r.theta_s_js_ns) - r.theta_s_js_ns)

	// // ステップnの境界jにおける表面熱流（壁体吸熱を正とする）のうち対流成分, W, [j, n+1]
	// var temp3, temp4 mat.Dense
	// temp3.Mul(ss.p_js_is, r.theta_r_is_ns)
	// temp3.Sub(&temp3, r.theta_s_js_ns)
	// temp4.MulElem(ss.h_s_c_js, ss.a_s_js)
	// temp4.MulElem(&temp4, &temp3)
	// r.q_c_js_ns = &temp4
	// // ss.h_s_c_js * ss.a_s_js * (np.dot(ss.p_js_is, r.theta_r_is_ns) - r.theta_s_js_ns)

	// // ---平均値・瞬時値---

	// // ステップnの室iにおける家具取得熱量, W, [i, n]
	// // ステップ n+1 の温度を用いてステップ n からステップ n+1 の平均的な熱流を求めている（後退差分）
	// r.q_frt_is_ns = np.delete(ss.g_sh_frt_is*(r.theta_r_is_ns-r.theta_frt_is_ns), 0, 1)

	// // ステップ n の室 i の家具等から空気への水分流, kg/s, [i, n]
	// // ステップ n+1 の湿度を用いてステップ n からステップ n+1 の平均的な水分流を求めている（後退差分）
	// r.q_l_frt_is_ns = np.delete(ss.g_lh_frt_is*(r.x_r_is_ns-r.x_frt_is_ns), 0, 1)

	// // ステップ n+1 のPMVを計算するのに、ステップ n からステップ n+1 のClo値を用いる。
	// // 現在、Clo値の配列数が1つ多いバグがあるため、適切な長さになるようにスライスしている。
	// // TODO: 本来であれば、助走期間における、n=-1 の時の値を用いないといけないが、とりあえず、配列最後の値を先頭に持ってきて代用している。
	// // NOTE: ***************** 以下のコメントは対応するGoへの移植が必要****************
	// //clo_pls = np.append(r.clo_is_ns[:, -1:], self.clo_is_ns, axis=1)[:, 0:self._n_step_i]
	// // ステップ n+1 のPMVを計算するのに、ステップ n からステップ n+1 の人体周りの風速を用いる。
	// // TODO: 本来であれば、助走期間における、n=-1 の時の値を用いないといけないが、とりあえず、配列最後の値を先頭に持ってきて代用している。
	// // NOTE: ***************** 以下のコメントは対応するGoへの移植が必要****************
	// // v_hum_pls = np.append(self.v_hum_is_ns[:, -1:], self.v_hum_is_ns, axis=1)

	// // ---瞬時値---

	// // ステップ n の室 i におけるPMV実現値, [i, n+1]
	// self.pmv_is_ns = get_pmv_is_n(
	// 	p_v_is_ns,
	// 	self.theta_r_is_ns,
	// 	self.theta_mrt_hum_is_ns,
	// 	clo_pls,
	// 	v_hum_pls,
	// 	ss.met_is,
	// )

	// // ステップ n の室 i におけるPPD実現値, [i, n+1]
	// self.ppd_is_ns = pmv.get_ppd_is_n(self.pmv_is_ns)
}

func (self *Recorder) recording(n int) {
	// 瞬時値の書き込み

	// NOTE: ***************** 以下のコメントは対応するGoへの移植が必要****************
	// if n >= -1 {

	// 	// 瞬時値出力のステップ番号
	// 	n_i = n + 1

	// 	// 次の時刻に引き渡す値
	// 	self.theta_r_is_ns[:, n_i] = kwargs["theta_r_is_n_pls"].flatten()
	// 	self.theta_mrt_hum_is_ns[:, n_i] = kwargs["theta_mrt_hum_is_n_pls"].flatten()
	// 	self.x_r_is_ns[:, n_i] = kwargs["x_r_is_n_pls"].flatten()
	// 	self.theta_frt_is_ns[:, n_i] = kwargs["theta_frt_is_n_pls"].flatten()
	// 	self.x_frt_is_ns[:, n_i] = kwargs["x_frt_is_n_pls"].flatten()
	// 	self.theta_ei_js_ns[:, n_i] = kwargs["theta_ei_js_n_pls"].flatten()
	// 	self.q_s_js_ns[:, n_i] = kwargs["q_s_js_n_pls"].flatten()

	// 	// 次の時刻に引き渡さない値
	// 	self.theta_ot[:, n_i] = kwargs["theta_ot_is_n_pls"].flatten()
	// 	self.theta_s_js_ns[:, n_i] = kwargs["theta_s_js_n_pls"].flatten()
	// 	self.theta_rear_js_ns[:, n_i] = kwargs["theta_rear_js_n"].flatten()
	// }

	// // 平均値・積算値の書き込み

	// if n >= 0 {

	// 	// 平均値出力のステップ番号
	// 	n_a = n

	// 	// 次の時刻に引き渡す値
	// 	r.operation_mode_is_ns[:, n_a] = kwargs["operation_mode_is_n"].flatten()

	// 	// 次の時刻に引き渡さない値
	// 	// 積算値
	// 	self.l_cs_is_ns[:, n_a] = kwargs["l_cs_is_n"].flatten()
	// 	self.l_rs_is_ns[:, n_a] = kwargs["l_rs_is_n"].flatten()
	// 	self.l_cl_is_ns[:, n_a] = kwargs["l_cl_is_n"].flatten()
	// 	// 平均値
	// 	self.h_hum_c_is_ns[:, n_a] = kwargs["h_hum_c_is_n"].flatten()
	// 	self.h_hum_r_is_ns[:, n_a] = kwargs["h_hum_r_is_n"].flatten()
	// 	self.q_hum_is_ns[:, n_a] = kwargs["q_hum_is_n"].flatten()
	// 	self.x_hum_is_ns[:, n_a] = kwargs["x_hum_is_n"].flatten()
	// 	self.v_reak_is_ns[:, n_a] = kwargs["v_leak_is_n"].flatten()
	// 	self.v_ntrl_is_ns[:, n_a] = kwargs["v_vent_ntr_is_n"].flatten()
	// 	self.v_hum_is_ns[:, n_a] = kwargs["v_hum_is_n"].flatten()
	// 	self.clo_is_ns[:, n_a] = kwargs["clo_is_n"].flatten()
	// }
}

func (self *Recorder) export_pd() {
	// NOTE: ***************** 以下のコメントは対応するGoへの移植が必要****************

	// // データインデックス（「瞬時値・平均値用」・「積算値用（開始時刻）」・「積算値用（終了時刻）」）を作成する。
	// date_index_a_end, date_index_a_start, date_index_i = self._get_date_index()

	// // dataframe を作成（瞬時値・平均値用）
	// df_a1 = pd.DataFrame(self._get_flat_data_a().T, self.get_header_a(), index=[date_index_a_start, date_index_a_end])

	// // 列入れ替え用の新しいヘッダーを作成
	// new_columns_a = list(itertools.chain.from_iterable(
	// 	[[self._get_room_header_name(id=i, name=column[1]) for column in self._output_list_room_a] for i in r._id_rm_is]
	// ))

	// // 列の入れ替え
	// df_a2 = df_a1.reindex(columns=new_columns_a)

	// // dataframe を作成（積算値用）
	// df_i1 = pd.DataFrame(data=self._get_flat_data_i().T, columns=self.get_header_i(), index=date_index_i)

	// // 列入れ替え用の新しいヘッダーを作成
	// new_columns_i = ["out_temp", "out_abs_humid"] + list(itertools.chain.from_iterable(
	// 	[[self._get_room_header_name(id=id, name=column[1]) for column in r._output_list_room_i] for id in r._id_rm_is]
	// )) + list(itertools.chain.from_iterable(
	// 	[[self._get_boundary_name(id=id, name=column[1]) for column in r._output_list_boundary_i] for id in r._id_bdry_js]
	// ))

	// // 列の入れ替え
	// df_i2 = df_i1.reindex(columns=new_columns_i)

	// return df_i2, df_a2
}

// func (self *Recorder) get_header_i() {

// 	return ["out_temp", "out_abs_humid"] \
// 		+ list(itertools.chain.from_iterable(
// 			[self._get_room_header_names(name=column[1]) for column in r._output_list_room_i]))\
// 		+ list(itertools.chain.from_iterable(
// 			[self._get_boundary_names(name=column[1]) for column in r._output_list_boundary_i]))
// }

// func (self *Recorder) get_header_a() {

// 	return list(itertools.chain.from_iterable(
// 		[r._get_room_header_names(name=column[1]) for column in r._output_list_room_a]))
// }

// func (self *Recorder) _get_flat_data_i() {
// 	/*[列（項目数）✕行（時系列）]のデータを作成する。
// 	Returns:

// 	*/

// 	// 出力リストに従って1つずつ記録された2次元のデータを縦に並べていき（この時点で3次元になる）、concatenate でフラット化する。
// 	// 先頭に外気温度と外気湿度の2つのデータを並べてある。他のデータが2次元データのため、
// 	// 外気温度と外気湿度のデータもあえて 1 ✕ n の2次元データにしてから統合してある。
// 	return np.concatenate(
// 		[[r.theta_o_ns], [r.x_o_ns]]
// 		+ [self.__dict__[column[0]] for column in r._output_list_room_i]
// 		+ [self.__dict__[column[0]] for column in r._output_list_boundary_i]
// 	)
// }

// func (self *Recorder) _get_flat_data_a() {
// 	/*[列（項目数）✕行（時系列）]のデータを作成する。

// 	Returns:

// 	*/

// 	// 出力リストに従って1つずつ記録された2次元のデータを縦に並べていき（この時点で3次元になる）、concatenate でフラット化する。
// 	return np.concatenate([r.__dict__[column[0]] for column in r._output_list_room_a])
// }

// func (self *Recorder) _get_date_index() {
// 	/*データインデックスを作成する。

// 	Returns:

// 	*/

// 	// pandas 用の時間間隔 freq 引数
// 	freq = self._itv.get_pandas_freq()

// 	// date time index 作成（瞬時値・平均値）
// 	date_index_i = pd.date_range(start="1/1/" + self.YEAR, periods=self._n_step_i, freq=freq, name="start_time")

// 	// date time index 作成（積算値）（start と end の2種類作成する）
// 	date_index_a_start = pd.date_range(start="1/1/" + self.YEAR, periods=self._n_step_a, freq=freq)
// 	date_index_a_end = date_index_a_start + dt.timedelta(minutes=15)
// 	date_index_a_start.name = "start_time"
// 	date_index_a_end.name = "end_time"

// 	return date_index_a_end, date_index_a_start, date_index_i
// }

// func _get_room_header_name(cls, id: int, name: str) {
// 	/*room 用のヘッダー名称を取得する。

// 	Args:
// 		id: room のID
// 		name: 出力項目名称

// 	Returns:
// 		ヘッダー名称
// 	*/

// 	return "rm" + str(id) + "_" + name
// }

// func (self *Recorder) _get_room_header_names(self, name: str) {
// 	/*room 用のヘッダー名称を室の数分取得する。

// 	Args:
// 		name: 出力項目名称

// 	Returns:
// 		ヘッダー名称のリスト
// 	*/
// 	return [r._get_room_header_name(id=id, name=name) for id in r._id_rm_is]
// }

// func _get_boundary_name(cls, id: int, name: str) {
// 	/*boundary 用のヘッダ名称を取得する。

// 	Args:
// 		pps: PreCalcParameters クラス
// 		j: boundary の ID
// 		name: 出力項目名称

// 	Returns:

// 	*/

// 	return "b" + str(id) + "_" + name
// }

// func _get_boundary_names(self, name: str) {
// 	/*boundary 用のヘッダ名称を boundary の数だけ取得する。

// 	Args:
// 		name: 出力項目名称

// 	Returns:

// 	*/
// 	return [r._get_boundary_name(id=id, name=name) for id in r._id_bdry_js]
// }
