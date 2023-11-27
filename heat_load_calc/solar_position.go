package heat_load_calc

import "math"

/*
ステップnにおける太陽位置を計算する。
*/

/*
   太陽位置を計算する

   Args:
       phi_loc: 緯度, rad
       lambda_loc: 経度, rad
       interval: 生成するデータの時間間隔であり、以下の文字列で指定する。
           1h: 1時間間隔
           30m: 30分間隔
           15m: 15分間隔

   Returns:
       タプル
           (1) 太陽高度, rad [n]
           (2) 太陽方位角, rad [n]
*/
func calc_solar_position(phi_loc float64, lambda_loc float64, interval Interval) ([]float64, []float64) {

	// 標準子午線(meridian), rad
	lambda_loc_mer := _get_lambda_loc_mer()

	// ステップnにおける年通算日（1/1を1とする）, [n]
	d_ns := _get_d_ns(interval)

	// 1968年との年差
	n := _get_n()

	// 平均軌道上の近日点通過日（暦表時による1968年1月1日正午基準の日差）, d
	d_0 := _get_d_0(n)

	// ステップnにおける平均近点離角, rad, [n]
	m_ns := _get_m_ns(d_ns, d_0)

	// ステップnにおける近日点と冬至点の角度, rad, [n]
	epsilon_ns := _get_epsilon_ns(m_ns, n)

	// ステップnにおける真近点離角, rad, [n]
	v_ns := _get_v_ns(m_ns)

	// ステップnにおける均時差, rad, [n]
	e_t_ns := _get_e_t_ns(m_ns, epsilon_ns, v_ns)

	// ステップnにおける赤緯, rad, [n]
	delta_ns := _get_delta_ns(epsilon_ns, v_ns)

	// ステップnにおける標準時, d, [n]
	// 1h: 0, 1.0, .... , 23.0, 0, 1.0, ...23.0
	// 30m: 0, 0.5, 1.0, 1.5, .... , 23.5, 0, 0.5, ...23.5
	// 15m: 0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, .... , 23.75, 0, 0.25, ...23.75
	t_m_ns := _get_t_m_ns(interval)

	// ステップnにおける時角, rad, [n]
	omega_ns := _get_omega_ns(t_m_ns, lambda_loc, lambda_loc_mer, e_t_ns)

	// ステップnにおける太陽高度, rad, [n]
	h_sun_ns := _get_h_sun_ns(phi_loc, omega_ns, delta_ns)

	// 太陽の位置が天頂にないか（天頂にある = False, 天頂にない = True）, [n]
	is_not_zenith_ns := _get_is_not_zenith_ns(h_sun_ns)

	// ステップnにおける太陽の方位角の正弦（太陽が天頂に無い場合のみに定義される）, [n]
	sin_a_sun_ns := _get_sin_a_sun_ns(delta_ns, h_sun_ns, omega_ns, is_not_zenith_ns)

	// ステップnにおける太陽の方位角の余弦（太陽が天頂に無い場合のみに定義される）, [n]
	cos_a_sun_ns := _get_cos_a_sun_ns(delta_ns, h_sun_ns, phi_loc, is_not_zenith_ns)

	// ステップnにおける太陽の方位角（太陽が天頂に無い場合のみに定義される）, rad, [n]
	a_sun_ns := _get_a_sun_ns(cos_a_sun_ns, sin_a_sun_ns, is_not_zenith_ns)

	return h_sun_ns, a_sun_ns
}

/*
   標準子午線を取得する。

   Returns:
       標準子午線における経度, rad

   Notes:
       式(14)
*/
func _get_lambda_loc_mer() float64 {
	// 標準子午線における経度を135°とする。
	return 135.0 * math.Pi / 180.0
}

/*
   ステップnにおける年通算日を取得する 年通算日（1/1を1とする）, d
   Args:
       interval: 生成するデータの時間間隔
   Returns:
       ステップnにおける年通算日, d [n]
   Notes:
       1月1日を1とする。
       返り値は、365✕24✕n_hour の長さの配列
       n_hour: 1時間を何ステップで区切るのか
           1h: 1
           30m: 2
           15m: 4
       出力イメージ （n_hour = 1 の場合）
       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,.....365,365,365,365

       式(13)
*/
func _get_d_ns(interval Interval) []float64 {
	// 1時間を分割するステップ数
	n_hour := interval.get_n_hour()

	d_ns := make([]float64, 365*24*n_hour)

	off := 0
	for i := 0; i < 365; i++ {
		dd := float64(i + 1)
		for j := 0; j < 24*n_hour; j++ {
			d_ns[off] = dd
			off++
		}
	}

	return d_ns
}

/*
   1968年との年差を計算する。

   Returns:
       1968年との年差, year
   Notes:
       式(12)
*/
func _get_n() int {
	// 太陽位置の計算においては1989年で行う。
	y := 1989

	return y - 1968
}

/*
   平均軌道上の近日点通過日を取得する。

   Args:
       n: 1968年との年差, year

   Returns:
       平均軌道上の近日点通過日（暦表時による1968年1月1日正午基準の日差）, d

   Notes:
       式(11)
*/
func _get_d_0(n int) float64 {
	return 3.71 + 0.2596*float64(n) - float64(int((n+3.0)/4.0))
}

/*
   平均近点離角を計算する。
   Args:
       d_ns: 年通算日（1/1を1とする）, d [n]
       d_0: 平均軌道上の近日点通過日（暦表時による1968年1月1日正午基準の日差）, d

   Returns:
       ステップnにおける平均近点離角, rad [n]

   Notes:
       式(10)
*/
func _get_m_ns(d_ns []float64, d_0 float64) []float64 {

	// 近点年（近日点基準の公転周期日数）
	d_ay := 365.2596

	// ステップnにおける平均近点離角, rad, [n]
	m_ns := make([]float64, len(d_ns))
	for i, d := range d_ns {
		m_ns[i] = 2 * math.Pi * (d - d_0) / d_ay
	}

	return m_ns
}

/*
   ステップnにおける近日点と冬至点の角度を計算する。

   Args:
       m_ns: 平均近点離角, rad [n]
       n: 1968年との年差

   Returns:
       ステップnにおける近日点と冬至点の角度, rad [n]

   Notes:
       式(9)
*/
func _get_epsilon_ns(m_ns []float64, n int) []float64 {
	ep_ns := make([]float64, len(m_ns))
	for i, m := range m_ns {
		ep_ns[i] = (12.3901 + 0.0172*(float64(n)+m/(2*math.Pi))) * math.Pi / 180.0
	}
	return ep_ns
}

/*
   ステップnにおける真近点離角を計算する。

   Args:
       m_ns: ステップnにおける平均近点離角, rad [n]

   Returns:
       ステップnにおける真近点離角, rad [n]
   Notes:
       式(8)
*/
func _get_v_ns(m_ns []float64) []float64 {
	v_ns := make([]float64, len(m_ns))
	for i, m := range m_ns {
		v_ns[i] = m + (1.914*math.Sin(m)+0.02*math.Sin(2*m))*math.Pi/180.0
	}
	return v_ns
}

/*

   Args:
       m_ns: ステップnにおける平均近点離角, rad, [n]
       epsilon_ns: ステップnにおける近日点と冬至点の角度, rad, [n]
       v_ns: ステップnにおける真近点離角, rad, [n]
   Returns:
       ステップnにおける均時差, rad, [n]
   Notes:
       式(7)
*/
func _get_e_t_ns(m_ns []float64, epsilon_ns []float64, v_ns []float64) []float64 {
	e_t_ns := make([]float64, len(m_ns))
	for i, m := range m_ns {
		e_t_ns[i] = (m - v_ns[i]) - math.Atan(0.043*math.Sin(2.0*(v_ns[i]+epsilon_ns[i]))/(1.0-0.043*math.Cos(2.0*(v_ns[i]+epsilon_ns[i]))))
	}
	return e_t_ns
}

/*
   ステップnにおける赤緯を計算する。

   Args:
       epsilon_ns: ステップnにおける近日点と冬至点の角度, rad, [n]
       v_ns: ステップnにおける真近点離角, rad, [n]

   Returns:
       ステップnにおける赤緯, rad [n]

   Notes:
       赤緯は -π/2 ～ 0 π/2 の値をとる
       式(6)
*/
func _get_delta_ns(epsilon_ns []float64, v_ns []float64) []float64 {

	// 北半球の冬至の日赤緯, rad
	const delta_0 = -23.4393 * math.Pi / 180.0

	// 赤緯, rad, [n]
	delta_ns := make([]float64, len(epsilon_ns))
	for i, epsilon := range epsilon_ns {
		delta_ns[i] = math.Asin(math.Cos(v_ns[i]+epsilon) * math.Sin(delta_0))
	}

	return delta_ns
}

/*
   ステップnにおける標準時を計算する
   Args:
       interval: 生成するデータの時間間隔
   Returns:
       ステップnにおける標準時, d, [n]
*/
func _get_t_m_ns(interval Interval) []float64 {

	// 1時間を何分割するか
	n_hour := interval.get_n_hour()

	// インターバル時間, h
	int_interval := interval.get_time()

	// 1h: 0, 1.0, .... , 23.0, 0, 1.0, ...23.0
	// 30m: 0, 0.5, 1.0, 1.5, .... , 23.5, 0, 0.5, ...23.5
	// 15m: 0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, .... , 23.75, 0, 0.25, ...23.75
	day := make([]float64, 24*n_hour)
	for i := 0; i < 24*n_hour; i++ {
		day[i] = float64(i) * int_interval
	}

	year := make([]float64, 365*24*n_hour)
	off := 0
	for d := 0; d < 365; d++ {
		for h := 0; h < 24*n_hour; h++ {
			year[off+h] = day[h]
		}
		off += 24 * n_hour
	}

	return year
}

/*
ステップnにおける時角を計算する。

Args:
    t_m_ns: ステップnにおける標準時, h, [n]
    lambda_loc: 経度, rad
    lambda_loc_mer: 標準時の地点の経度, rad
    e_t_ns: ステップnにおける均時差, rad, [n]

Returns:
    ステップnにおける時角, rad, [n]

Notes:
    式(5)
*/
func _get_omega_ns(t_m_ns []float64, lambda_loc float64, lambda_loc_mer float64, e_t_ns []float64) []float64 {
	omega_ns := make([]float64, len(t_m_ns))
	for i, t_m := range t_m_ns {
		omega_ns[i] = ((t_m-12.0)*15.0)*math.Pi/180.0 + (lambda_loc - lambda_loc_mer) + e_t_ns[i]
	}
	return omega_ns
}

/*
ステップnにおける太陽高度を計算する。

Args:
    phi_loc: 経度, rad
    omega_ns: ステップnにおける時角, rad, [n]
    delta_ns: ステップnにおける赤緯, rad, [n]

Returns:
    ステップnにおける太陽高度, rad, [n]

Notes:
    太陽高度はマイナスの値もとり得る。（太陽が沈んでいる場合）
    式(4)
*/
func _get_h_sun_ns(phi_loc float64, omega_ns []float64, delta_ns []float64) []float64 {
	h_sun_ns := make([]float64, len(omega_ns))
	sin_phi_loc := math.Sin(phi_loc)
	cos_phi_loc := math.Cos(phi_loc)
	for i, omega := range omega_ns {
		h_sun_ns[i] = math.Asin(sin_phi_loc*math.Sin(delta_ns[i]) + cos_phi_loc*math.Cos(delta_ns[i])*math.Cos(omega))
	}
	return h_sun_ns
}

/*
Args:
    h_sun_ns: ステップnにおける太陽高度, rad [n]

Returns:
    太陽の位置が天頂にないか（天頂にある = False, 天頂にない = True）, [n]
*/
func _get_is_not_zenith_ns(h_sun_ns []float64) []bool {
	is_not_zenith_ns := make([]bool, len(h_sun_ns))
	for i, h_sun := range h_sun_ns {
		is_not_zenith_ns[i] = h_sun != math.Pi/2
	}
	return is_not_zenith_ns
}

/*

Args:
    delta_ns: ステップnにおける赤緯, rad [n]
    h_sun_ns: ステップnにおける太陽高度, rad [n]
    omega_ns: ステップnにおける時角, rad [n]
    inzs: ステップ n における太陽位置が天頂にあるか否か（True=天頂にない, False=天頂にある）

Returns:
    ステップnにおける太陽の方位角の正弦（太陽が天頂に無い場合のみに定義される）, [n]

Notes:
    式(3)
*/
func _get_sin_a_sun_ns(delta_ns, h_sun_ns, omega_ns []float64, inzs []bool) []float64 {

	sin_a_sun_ns := make([]float64, len(h_sun_ns))
	for i, inz := range inzs {
		if inz {
			// 太陽の方位角の正弦（太陽が天頂に無い場合のみ計算する）
			sin_a_sun_ns[i] = math.Cos(delta_ns[i]) * math.Sin(omega_ns[i]) / math.Cos(h_sun_ns[i])
		} else {
			// 太陽が天頂にある場合は「定義なし = np.nan」とする
			sin_a_sun_ns[i] = math.NaN()
		}
	}

	return sin_a_sun_ns
}

/*

Args:
    delta_ns: ステップnにおける赤緯, rad [n]
    h_sun_ns: ステップnにおける太陽高度, rad [n]
    phi_loc: 緯度, rad
    inzs: ステップ n における太陽位置が天頂にあるか否か（True=天頂にない, False=天頂にある）

Returns:
    ステップnにおける太陽の方位角の余弦（太陽が天頂に無い場合のみに定義される）, [n]

Notes:
    式(2)

*/
func _get_cos_a_sun_ns(delta_ns, h_sun_ns []float64, phi_loc float64, inzs []bool) []float64 {
	cos_a_sun_ns := make([]float64, len(h_sun_ns))
	sin_phi_loc := math.Sin(phi_loc)
	cos_phi_loc := math.Cos(phi_loc)
	for i, inz := range inzs {
		if inz {
			// 太陽の方位角の余弦（太陽が天頂に無い場合のみ計算する）
			cos_a_sun_ns[i] = (math.Sin(h_sun_ns[i])*sin_phi_loc - math.Sin(delta_ns[i])) /
				(math.Cos(h_sun_ns[i]) * cos_phi_loc)
		} else {
			// 太陽が天頂にある場合は「定義なし = np.nan」とする
			cos_a_sun_ns[i] = math.NaN()
		}
	}

	return cos_a_sun_ns
}

/*
Args:
    cos_a_sun_ns: ステップ n における太陽の方位角の余弦（太陽が天頂に無い場合のみ計算する）
    sin_a_sun_ns: ステップ n における太陽の方位角の正弦（太陽が天頂に無い場合のみ計算する）
    inzs: ステップ n における太陽位置が天頂にあるか否か（True=天頂にない, False=天頂にある）

Returns:
    ステップ n における太陽の方位角（太陽が天頂に無い場合のみに定義される）, rad, [n]

Notes:
    式(1)
*/
func _get_a_sun_ns(cos_a_sun_ns, sin_a_sun_ns []float64, inzs []bool) []float64 {
	a_sun_ns := make([]float64, len(cos_a_sun_ns))
	for i, inz := range inzs {
		if inz {
			// 太陽の方位角, rad [n] （太陽が天頂に無い場合のみ計算する）
			// arctan の注意点。
			// arctan2 は、座標上で第1引数をy, 第2引数をxにした際にx軸との角度を求める関数。
			// 従って、単射の通常良く用いられる -π/2 ～ 0 ～ π/2 ではない。
			//   sin_a_s_ns が正 かつ cos_a_s_ns が正 の場合は第1象限（0～π/2）
			//   sin_a_s_ns が正 かつ cos_a_s_ns が負 の場合は第2象限（π/2～π）
			//   sin_a_s_ns が負 かつ cos_a_s_ns が負 の場合は第3象限（-π～-π/2）
			//   sin_a_s_ns が負 かつ cos_a_s_ns が正 の場合は第4象限（-π/2～0）
			a_sun_ns[i] = math.Atan2(sin_a_sun_ns[i], cos_a_sun_ns[i])
		} else {
			// 太陽が天頂にある場合は「定義なし = np.nan」とする
			a_sun_ns[i] = math.NaN()
		}
	}

	return a_sun_ns
}
