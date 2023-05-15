package main

import (
	"errors"
)

// 備品等の物性値を計算する

/*
備品等に関する物性値を取得する。
Args:

	dict_furniture_i: 室 i の備品等に関する入力情報
	v_rm_i: 室 i の容量, m3

Returns:

	備品等に関する物性値
		室 i の備品等の熱容量, J/K
		室 i の空気と備品等間の熱コンダクタンス, W/K
		室 i の備品等の湿気容量, kg/(kg/kgDA)
		室 i の空気と備品等間の湿気コンダクタンス, kg/(s (kg/kgDA))

Notes:

	各値は、特定の値の入力を受け付ける他に、以下の式(1)～(4)により室容積から推定する方法を設定する。
*/
func get_furniture_specs(dict_furniture_i *FurnitureJson, v_rm_i float64) (float64, float64, float64, float64, error) {

	var c_sh_frt_i, g_sh_frt_i, c_lh_frt_i, g_lh_frt_i float64

	if dict_furniture_i.InputMethod == "default" {
		c_sh_frt_i = _get_c_sh_frt_i(v_rm_i)
		g_sh_frt_i = _get_g_sh_frt_i(c_sh_frt_i)
		c_lh_frt_i = _get_c_lh_frt_i(v_rm_i)
		g_lh_frt_i = _get_g_lh_frt_i(c_lh_frt_i)
	} else if dict_furniture_i.InputMethod == "specify" {
		c_sh_frt_i = dict_furniture_i.HeatCapacity
		g_sh_frt_i = dict_furniture_i.HeatCond
		c_lh_frt_i = dict_furniture_i.MoistureCapacity
		g_lh_frt_i = dict_furniture_i.MoistureCond
	} else {
		return 0, 0, 0, 0, errors.New("invalid input method")
	}

	return c_lh_frt_i, c_sh_frt_i, g_lh_frt_i, g_sh_frt_i, nil
}

/*
備品等の熱容量を計算する。
Args:

	v_rm_i: 室iの気積, m3

Returns:

	室 i の備品等の熱容量, J/K

Notes:

	式(4)
*/
func _get_c_sh_frt_i(v_rm_i float64) float64 {
	// 室の気積あたりの備品等の熱容量を表す係数, J/(m3 K)
	const f_c_sh_frt = 12.6 * 1000.0
	return f_c_sh_frt * v_rm_i
}

/*
空気と備品等間の熱コンダクタンスを取得する。
Args:

	c_sh_frt_i: 室 i の備品等の熱容量, J/K

Returns:

	室 i の空気と備品等間の熱コンダクタンス, W/K

Notes:

	式(3)
*/
func _get_g_sh_frt_i(c_sh_frt_i float64) float64 {
	// 備品等の熱容量あたりの空気との間の熱コンダクタンスを表す係数, 1/s
	const f_g_sh_frt = 0.00022
	return f_g_sh_frt * c_sh_frt_i
}

/*
備品等の湿気容量を計算する。
Args:

	v_rm_i: 室iの気積, m3

Returns:

	室 i の備品等の湿気容量, kg/(kg/kg(DA))

Notes:

	式(2)
*/
func _get_c_lh_frt_i(v_rm_i float64) float64 {
	// 室の気積あたりの備品等の湿気容量を表す係数, kg/(m3 kg/kg(DA))
	const f_c_lh_frt = 16.8
	return f_c_lh_frt * v_rm_i
}

/*
空気と備品等間の湿気コンダクタンスを取得する。
Args:

	c_lh_frt_i: 室iの備品等の湿気容量, kg/(kg/kg(DA))

Returns:

	室iの空気と備品等間の湿気コンダクタンス, kg/(s kg/kg(DA))

Notes:

	式(1)
*/
func _get_g_lh_frt_i(c_lh_frt_i float64) float64 {
	// 備品等の湿気容量あたりの空気との間の湿気コンダクタンスを表す係数, 1/s
	const f_g_lh_frt = 0.0018
	return f_g_lh_frt * c_lh_frt_i
}
