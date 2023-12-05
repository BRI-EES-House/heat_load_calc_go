package heat_load_calc

import (
	"gonum.org/v1/gonum/mat"
)

/*
   間仕切りの相当外気温度を計算する。

   Args:
       w: Weather クラス

   Returns:
       ステップ n におけるの境界 j の相当外気温度, degree C, [N+1]

   Notes:
       本来であれば、間仕切り壁において相当外気温度は定義できない概念ではあるが、
       後々の numpy の行列計算において、相当外気温度を 0.0 にしておく方が都合がよいので、
       ここでは0.0に初期化することにした。
       間仕切り壁に対する外気温度の影響はゼロなので、ここで定める値は何に設定しても用いられることはなく、計算結果に影響しない。

*/
func get_theta_o_eqv_j_ns_for_internal(w *Weather) *mat.VecDense {
	// eq.1
	return mat.NewVecDense(w.number_of_data_plus(), nil)
}

/*
   相当外気温度を計算する。

   Args:
       drct_j: 境界jのDirectionクラス
       a_s_j: 境界jの室外側日射吸収率, -
       eps_r_o_j: 境界jの室外側長波長放射率, -
       r_s_o_j: 境界jの室外側熱伝達抵抗, m2K/W
       ssp_j: 境界jのSolarShadingPartクラス
       w: Weather クラス

   Returns:
       ステップ n における境界 j の相当外気温度, degree C, [N+1]
   Notes:
       eq.2
*/
func get_theta_o_eqv_j_ns_for_external_general_part_and_external_opaque_part(
	drct_j Direction,
	a_s_j float64,
	eps_r_o_j float64,
	r_s_o_j float64,
	ssp_j SolarShading,
	w *Weather,
) *mat.VecDense {
	// ステップnにおける境界jの直達日射に対する日よけの影面積比率, -, [N+
	f_ss_dn_j_ns := ssp_j.get_f_ss_dn_j_ns(w.h_sun_ns_plus, w.a_sun_ns_plus)

	// ステップnにおける境界jの天空日射に対する日よけの影面積比率, -
	f_ss_sky_j_ns := ssp_j.get_f_ss_sky_j()

	// ステップnにおける境界jの地面反射日射に対する日よけの影面積比率, -
	f_ss_ref_j_ns := ssp_j.get_f_ss_ref_j()

	// ステップ n における境界 j の傾斜面に入射する日射量の直達成分, W/m2 [n]
	// ステップ n における境界 j の傾斜面に入射する日射量の天空成分, W/m2 [n]
	// ステップ n における境界 j の傾斜面に入射する日射量の地盤反射成分, W/m2 [n]
	// ステップ n における境界 j の傾斜面の夜間放射量, W/m2, [n]
	i_s_dn_j_ns, i_s_sky_j_ns, i_s_ref_j_ns, r_s_n_j_ns := get_i_is_j_ns(w, drct_j)

	// ステップnにおける境界jの相当外気温度, degree C, [N+1]
	// 一般部位・不透明な開口部の場合、日射・長波長放射を考慮する。
	// eq.2
	theta_o_eqv_j_ns := mat.NewVecDense(w.number_of_data_plus(), nil)
	theta_o_ns_plus := w.theta_o_ns_plus
	for n := 0; n < w.number_of_data_plus(); n++ {
		theta_o_eqv_j_ns.SetVec(
			n,
			theta_o_ns_plus[n]+
				(a_s_j*(i_s_dn_j_ns[n]*(1.0-f_ss_dn_j_ns[n])+
					i_s_sky_j_ns[n]*(1.0-f_ss_sky_j_ns)+
					i_s_ref_j_ns[n]*(1.0-f_ss_ref_j_ns))-
					eps_r_o_j*r_s_n_j_ns[n])*r_s_o_j)
	}

	return theta_o_eqv_j_ns
}

/*
   透明な開口部の場合の相当外気温度を計算する。

   Args:
       w: Weather クラス
       drct_j: Direction クラス
       eps_r_o_j: 室外側長波長放射率, -
       r_s_o_j: 境界jの室外側熱伝達抵抗, m2K/W
       u_j: 境界jの熱貫流率, W/m2K
       ssp_j: SolarShading クラス
       wdw: Window クラス
       w: Weather クラス

   Returns:
       ステップ n における境界 j の相当外気温度, degree C, [N+1]
   Notes:
       eq.3~7
*/
func get_theta_o_eqv_j_ns_for_external_transparent_part(
	drct_j Direction,
	eps_r_o_j float64,
	r_s_o_j float64,
	u_j float64,
	ssp_j SolarShading,
	wdw_j *Window,
	w *Weather,
) *mat.VecDense {
	// ステップ n の境界 j における傾斜面に入射する太陽の入射角, rad, [n]
	cos_phi_j_ns := w.cos_phi_j_ns[drct_j]

	// ステップ n における境界 j の傾斜面に入射する日射量の直達成分, W / m2, [n]
	// ステップ n における境界 j の傾斜面に入射する日射量の天空成分, W / m2, [n]
	// ステップ n における境界 j の傾斜面に入射する日射量の地盤反射成分, W / m2, [n]
	// ステップ n における境界 j の傾斜面の夜間放射量, W/m2, [n]
	i_s_dn_j_ns, i_s_sky_j_ns, i_s_ref_j_ns, r_s_n_j_ns := get_i_is_j_ns(w, drct_j)

	// ---日よけの影面積比率

	// 直達日射に対する日よけの影面積比率, [N+1]
	f_ss_d_j_ns := ssp_j.get_f_ss_dn_j_ns(w.h_sun_ns_plus, w.a_sun_ns_plus)

	// 天空日射に対する日よけの影面積比率
	f_ss_s_j_ns := ssp_j.get_f_ss_sky_j()

	// 地面反射日射に対する日よけの影面積比率
	f_ss_r_j_ns := ssp_j.get_f_ss_ref_j()

	// 境界 ｊ　の開口部の天空日射に対する吸収日射熱取得率, -
	b_w_s_j := wdw_j.b_w_s_j()

	// 境界 ｊ　の開口部の地盤反射日射に対する吸収日射熱取得率, -
	b_w_r_j := wdw_j.b_w_r_j()

	theta_o_ns_plus := w.theta_o_ns_plus
	theta_o_sol_i_j_ns := make([]float64, w.number_of_data_plus())
	for i := 0; i < w.number_of_data_plus(); i++ {
		// ステップ n における境界 ｊ　の開口部の直達日射に対する吸収日射熱取得率, -, [N+1]
		b_w_d_j_ns := wdw_j.get_b_w_d_j_ns(cos_phi_j_ns[i])

		// 直達日射に起因する吸収日射熱取得量, W/m2, [N+1]
		// eq.5
		q_b_dn_j_ns := b_w_d_j_ns * (1.0 - f_ss_d_j_ns[i]) * i_s_dn_j_ns[i]

		// 天空日射に起因する吸収日射熱取得量, W/m2, [N+1]
		// eq.6
		q_b_sky_j_ns := b_w_s_j * (1.0 - f_ss_s_j_ns) * i_s_sky_j_ns[i]

		// 地盤反射日射に起因する吸収日射熱取得量, W/m2, [N+1]
		// eq.7
		q_b_ref_j_ns := b_w_r_j * (1.0 - f_ss_r_j_ns) * i_s_ref_j_ns[i]

		// 吸収日射熱取得量, W/m2, [N+1]
		// eq.4
		q_b_all_j_ns := (q_b_dn_j_ns + q_b_sky_j_ns + q_b_ref_j_ns)

		// 室iの境界jの傾斜面のステップnにおける相当外気温度, degree C, [N+1]
		// 透明な開口部の場合、透過日射はガラス面への透過の項で扱うため、ここでは吸収日射、長波長放射のみ考慮する。
		// eq.3
		theta_o_sol_i_j_ns[i] = theta_o_ns_plus[i] - eps_r_o_j*r_s_n_j_ns[i]*r_s_o_j + q_b_all_j_ns/u_j
	}

	return mat.NewVecDense(w.number_of_data_plus(), theta_o_sol_i_j_ns)
}

/*
   日射があたっていない場合の外皮の相当外気温度を計算する。

   Args:
       w: Weather クラス

   Returns:
       ステップ n における境界 j の相当外気温度, degree C, [N+1]
   Notes:
       eq.8
*/
func get_theta_o_eqv_j_ns_for_external_not_sun_striked(w *Weather) *mat.VecDense {
	theta_o_ns_plus := w.theta_o_ns_plus
	return mat.NewVecDense(len(theta_o_ns_plus), theta_o_ns_plus)
}

/*
	地盤の相当外気温度を計算する。

    Args:
        w: Weather クラス

    Returns:
        ステップ n における室 i の境界 j の傾斜面の相当外気温度, degree C, [N+1]
*/
func get_theta_o_eqv_j_ns_for_ground(w *Weather) *mat.VecDense {
	theta_o_ns_average := w.get_theta_o_ave()
	data := make([]float64, w.number_of_data_plus())

	for i := range data {
		data[i] = theta_o_ns_average
	}

	return mat.NewVecDense(w.number_of_data_plus(), data)
}
