package heat_load_calc

/*
「透明な開口部」以外の境界又は「透明な開口部」であって日射が当たらない場合のステップnにおける境界jの透過日射量を計算する。
    Args:
        w: Weather Class
    Returns:
        ステップnにおける境界jの透過日射量, W, [N+1]
*/
func get_q_trs_sol_j_ns_for_not(w *Weather) []float64 {
	return make([]float64, w.number_of_data_plus())
}

/*
   Args:
       drct_j: 境界jの Direction Class
       a_s_j: 境界jの面積, m2
       ssp_j: 境界jの SolarShadingPart Class
       wdw_j: 境界jの Window Class
       w: Weather Class
   Returns:
       ステップ n における境界 j の透過日射量, W, [N+1]
*/
func get_q_trs_sol_j_ns_for_transparent_sun_striked(
	drct_j Direction,
	a_s_j float64,
	ssp_j SolarShading,
	wdw_j *Window,
	w *Weather,
) []float64 {
	// ステップnにおける境界jの傾斜面に入射する太陽の入射角, rad, [N+1]
	cos_phi_j_ns := w.cos_phi_j_ns[drct_j]

	// ステップ n における境界 j の傾斜面に入射する日射量のうち直達成分, W/m2 [N+1]
	// ステップ n における境界 j の傾斜面に入射する日射量のうち天空成分, W/m2 [N+1]
	// ステップ n における境界 j の傾斜面に入射する日射量のうち地盤反射成分, W/m2 [N+1]
	// ステップ n における境界 j の傾斜面の夜間放射量, W/m2, [N+1]
	i_s_dn_j_ns, i_s_sky_j_ns, i_s_ref_j_ns, _ := get_i_is_j_ns(w, drct_j)

	// ---日よけの影面積比率

	// ステップnにおける境界jの直達日射に対する日よけの影面積比率, [N+1]
	f_ss_d_j_ns := ssp_j.get_f_ss_dn_j_ns(w.h_sun_ns_plus, w.a_sun_ns_plus)

	// ステップnにおける境界jの天空日射に対する日よけの影面積比率
	f_ss_s_j_ns := ssp_j.get_f_ss_sky_j()

	// ステップnにおける境界jの地面反射日射に対する日よけの影面積比率
	f_ss_r_j_ns := ssp_j.get_f_ss_ref_j()

	// 境界jの窓の天空日射に対する日射透過率, -
	tau_w_s_j := wdw_j.tau_w_s_j()

	// 境界jの窓の地盤反射日射に対する日射透過率, -
	tau_w_r_j := wdw_j.tau_w_r_j()

	q_trs_sol_j_ns := make([]float64, len(cos_phi_j_ns))
	for i := 0; i < len(cos_phi_j_ns); i++ {
		// ステップ n における境界 j の窓の直達日射に対する日射透過率, -, [N+1]
		tau_w_d_j_ns := wdw_j.get_tau_w_d_j_ns(cos_phi_j_ns[i])

		// ステップ n における境界 j の窓の直達日射に対する透過日射量, W/m2, [N+1]
		q_trs_sol_dn_j_ns := tau_w_d_j_ns * (1.0 - f_ss_d_j_ns[i]) * i_s_dn_j_ns[i]

		// ステップ n における境界 j の天空日射に対する透過日射量, W/m2, [N+1]
		q_trs_sol_sky_j_ns := tau_w_s_j * (1.0 - f_ss_s_j_ns) * i_s_sky_j_ns[i]

		// ステップ n における境界 j の地盤反射日射に対する透過日射量, W/m2, [N+1]
		q_trs_sol_ref_j_ns := tau_w_r_j * (1.0 - f_ss_r_j_ns) * i_s_ref_j_ns[i]

		// ステップ n における境界 j の透過日射量, W, [N+1]
		q_trs_sol_j_ns[i] = (q_trs_sol_dn_j_ns + q_trs_sol_sky_j_ns + q_trs_sol_ref_j_ns) * a_s_j
	}

	return q_trs_sol_j_ns
}
