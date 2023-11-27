package heat_load_calc

import (
	"math"
)

// フレーム材質
type FlameType int

const (
	FlameTypeRESIN       FlameType = iota // 樹脂製製建具
	FlameTypeWOOD                         // 木製建具
	FlameTypeMIXED_WOOD                   // 木と金属の複合材料製建具
	FlameTypeMIXED_RESIN                  // 樹脂と金属の複合材料製建具
	FlameTypeALUMINUM                     // 金属製建具
)

// ガラスの構成
type GlassType string

// ガラスの構成
const (
	GlassTypeSINGLE   GlassType = "single"   // 単層
	GlassTypeMULTIPLE GlassType = "multiple" // 複層
)

func GlassTypeFromString(str string) GlassType {
	switch str {
	case "single":
		return GlassTypeSINGLE
	case "multiple":
		return GlassTypeMULTIPLE
	default:
		panic("invalid glazing type")
	}
}

// 窓
type Window struct {
	_glass_type    GlassType // ガラスの構成
	_u_w_f_j       float64   // 建具部分の熱損失係数（U値）
	_r_a_w_g_j     float64   // 境界 j の開口部の面積に対するグレージングの面積の比, -
	_u_w_g_j       float64
	_eta_w_g_j     float64 // 境界 j の開口部の日射熱取得率, -
	_u_w_g_s_j     float64
	_r_r_w_g_j     float64
	_rho_w_g_s1f_j float64
	_rho_w_g_s2f_j float64
	_tau_w_g_j     float64
	_tau_w_g_s1_j  float64
	_tau_w_g_s2_j  float64
	_rho_w_g_s1b_j float64
	_tau_w_c_j     float64
	_b_w_c_j       float64
}

/*
   Args:
		u_w_j: 境界 j の窓の熱損失係数, W/m2K
		eta_w_j: 境界 j の窓の日射熱取得率, -
		glass_type: 境界 j の窓のガラス構成
		r_a_w_g_j: 境界 j の窓の面積に対するグレージングの面積の比, -
		flame_type: _description_. funcaults to FlameType.MIXED.
*/
func NewWindow(
	u_w_j float64,
	eta_w_j float64,
	glass_type GlassType,
	r_a_w_g_j float64,
	flame_type FlameType,
) *Window {
	u_w_f_j := _get_u_w_f_j(flame_type)
	r_a_w_g_j = _get_r_a_w_g_j(r_a_w_g_j, flame_type)
	u_w_g_j := _get_u_w_g_j(u_w_j, u_w_f_j, r_a_w_g_j)
	eta_w_g_j := _get_eta_w_g_j(eta_w_j, r_a_w_g_j)
	r_w_o_w, r_w_i_w, r_w_o_s, r_w_i_s := _get_r_w()
	u_w_g_s_j := _get_u_w_g_s_j(u_w_g_j, r_w_o_w, r_w_i_w, r_w_o_s, r_w_i_s)
	r_r_w_g_j := _get_r_r_w_g_j(
		u_w_g_j,
		u_w_g_s_j,
		r_w_o_w,
		r_w_i_w,
		r_w_o_s,
		glass_type,
	)
	rho_w_g_s1f_j := _get_rho_w_g_s1f_j(r_r_w_g_j, eta_w_g_j)
	rho_w_g_s2f_j := _get_rho_w_g_s2f_j(glass_type)
	tau_w_g_j := _get_tau_w_g_j(
		eta_w_g_j,
		r_r_w_g_j,
		rho_w_g_s1f_j,
		rho_w_g_s2f_j,
		glass_type,
	)
	tau_w_g_s1_j := _get_tau_w_g_s1_j(tau_w_g_j, rho_w_g_s2f_j, glass_type)
	tau_w_g_s2_j := _get_tau_w_g_s2_j(tau_w_g_s1_j, glass_type)
	rho_w_g_s1b_j := _get_rho_w_g_s1b_j(tau_w_g_s1_j, glass_type)

	w := &Window{
		_glass_type:    glass_type,
		_u_w_f_j:       u_w_f_j,
		_r_a_w_g_j:     r_a_w_g_j,
		_u_w_g_j:       u_w_g_j,
		_eta_w_g_j:     eta_w_g_j,
		_u_w_g_s_j:     u_w_g_s_j,
		_r_r_w_g_j:     r_r_w_g_j,
		_rho_w_g_s1f_j: rho_w_g_s1f_j,
		_rho_w_g_s2f_j: rho_w_g_s2f_j,
		_tau_w_g_j:     tau_w_g_j,
		_tau_w_g_s1_j:  tau_w_g_s1_j,
		_tau_w_g_s2_j:  tau_w_g_s2_j,
		_rho_w_g_s1b_j: rho_w_g_s1b_j,
	}

	tau_w_c_j, b_w_c_j := w._get_tau_b_w_c_j()

	w._tau_w_c_j = tau_w_c_j
	w._b_w_c_j = b_w_c_j

	return w
}

/*
窓の直達日射に対する日射透過率を計算する。
Args:
	phi_j_ns: ステップnにおける境界jの窓の直達日射の入射角, rad, [N+1]
Returns:
	ステップnにおける境界jの窓の直達日射に対する日射透過率, -, [N+1]
Notes:
	eq.1
*/
func (w *Window) get_tau_w_d_j_ns(cos_phi float64) float64 {
	return w._get_tau_w_j_phis(cos_phi)
}

/*
窓の直達日射に対する吸収日射熱取得率を計算する。
Args:
	phi_j_ns: ステップnにおける境界jの窓の直達日射の入射角, rad, [N+1]
Returns:
	ステップnにおける境界jの窓の直達日射に対する吸収日射熱取得率, [n]
Notes:
	eq.2
*/
func (w *Window) get_b_w_d_j_ns(cos_phi_j_ns float64) float64 {
	return w._get_b_w_j_phis(cos_phi_j_ns)
}

/*
Returns:
	境界jの窓の天空日射に対する日射透過率, -
Notes:
	eq.3
*/
func (w *Window) tau_w_s_j() float64 {
	return w._tau_w_c_j
}

/*
Returns:
	境界jの窓の地面反射日射に対する日射透過率, -
Notes:
	eq.4
*/
func (w *Window) tau_w_r_j() float64 {
	return w._tau_w_c_j
}

/*
Returns:
	境界jの窓の天空日射に対する吸収日射熱取得率, -
Notes:
	eq.5
*/
func (w *Window) b_w_s_j() float64 {
	return w._b_w_c_j
}

/*
窓の地面反射に対する吸収日射熱取得率を取得する。
Returns:
	境界jの窓の地面反射日射に対する吸収日射熱取得率, -
Notes:
	eq.6
*/
func (w *Window) b_w_r_j() float64 {
	return w._b_w_c_j
}

/*
境界jの窓の1/4球の一様拡散日射に対する日射透過率及び吸収日射熱取得率を計算する。
Returns:
	境界jの窓の1/4球の一様拡散日射に対する日射透過率, -
	境界jの窓の1/4球の一様拡散日射に対する吸収日射熱取得率, -
Notes:
	eq.7
	eq.8
*/
func (w *Window) _get_tau_b_w_c_j() (float64, float64) {
	const M = 1000

	var tau_term, b_term float64
	for i := 1; i <= M; i++ {
		ms := float64(i)
		phi_ms := math.Pi / 2.0 * (ms - 1.0/2.0) / M
		cos_phi_ms := math.Cos(phi_ms)
		sin_cos := math.Sin(phi_ms) * cos_phi_ms
		tau := w._get_tau_w_j_phis(cos_phi_ms)
		b := w._get_b_w_j_phis(cos_phi_ms)
		tau_term += tau * sin_cos
		b_term += b * sin_cos
	}

	tau_w_c_j := math.Pi / M * tau_term
	b_w_c_j := math.Pi / M * b_term

	return tau_w_c_j, b_w_c_j
}

/*
入射角Φに対する境界jの窓の日射透過率を取得する。
Args:
	phis: 入射角Φ, rad
Returns:
	入射角Φに対する境界jの窓の日射透過率, -
Notes:
	eq.9
*/
func (w *Window) _get_tau_w_j_phis(cos_phi float64) float64 {
	return w._get_tau_w_g_j_phis(cos_phi) * w._r_a_w_g_j
}

/*
入射角Φに対する境界jの窓の吸収日射熱取得率を取得する。
Args:
	phis: 入射角Φ, rad
Returns:
	入射角Φに対する境界jの窓の吸収日射熱取得率, -
Notes:
	eq.10
*/
func (w *Window) _get_b_w_j_phis(cos_phis float64) float64 {
	return w._get_b_w_g_j_phis(cos_phis) * w._r_a_w_g_j
}

/*
入射角Φに対する境界jの窓のガラス部分の吸収日射熱取得率を取得する。
Args:
	phis: ステップnにおける入射角, rad
Returns:
	入射角Φに対する境界jの窓のガラス部分の吸収日射熱取得率, -
Notes:
	eq.11
*/
func (w *Window) _get_b_w_g_j_phis(cos_phi float64) float64 {
	tau_w_g_j_phis := w._get_tau_w_g_j_phis(cos_phi)
	rho_w_g_j_phis := w._get_rho_w_g_j_phis(cos_phi)
	return (1 - tau_w_g_j_phis - rho_w_g_j_phis) * w._r_r_w_g_j
}

/*
入射角Φに対する境界jの窓のガラス部分の日射透過率を計算する。
Args:
	phis: 入射角Φ, rad
Returns:
	入射角Φに対する境界jの窓のガラス部分の日射透過率, -
Notes:
	eq.12
*/
func (w *Window) _get_tau_w_g_j_phis(cos_phi float64) float64 {
	if w._glass_type == GlassTypeSINGLE {
		tau_w_g_s1_j_phis := w._get_tau_w_g_s1_j_phis(cos_phi)
		return tau_w_g_s1_j_phis
	} else if w._glass_type == GlassTypeMULTIPLE {
		rho_w_g_s1b_j_phis := w._get_rho_w_g_s1b_j_phis(cos_phi)
		rho_w_g_s2f_j_phis := w._get_rho_w_g_s2f_j_phis(cos_phi)
		tau_w_g_s1_j_phis := w._get_tau_w_g_s1_j_phis(cos_phi)
		tau_w_g_s2_j_phis := w._get_tau_w_g_s2_j_phis(cos_phi)
		return tau_w_g_s1_j_phis * tau_w_g_s2_j_phis / (1 - math.Min(rho_w_g_s1b_j_phis*rho_w_g_s2f_j_phis, 0.9999))
	} else {
		panic(w._glass_type)
	}
}

/*
入射角Φに対する境界jの窓のガラス部分の日射反射率を計算する。
Args:
	phis: 入射角, rad
Returns:
	入射角Φに対する境界jの窓のガラス部分の日射反射率, -
Notes:
	eq.13
*/
func (w *Window) _get_rho_w_g_j_phis(cos_phi float64) float64 {
	if w._glass_type == GlassTypeSINGLE {
		rho_w_g_s1f_j_phis := w._get_rho_w_g_s1f_j_phis(cos_phi)
		return rho_w_g_s1f_j_phis
	} else if w._glass_type == GlassTypeMULTIPLE {
		rho_w_g_s1f_j_phis := w._get_rho_w_g_s1f_j_phis(cos_phi)
		rho_w_g_s1b_j_phis := w._get_rho_w_g_s1b_j_phis(cos_phi)
		rho_w_g_s2f_j_phis := w._get_rho_w_g_s2f_j_phis(cos_phi)
		tau_w_g_s1_j_phis := w._get_tau_w_g_s1_j_phis(cos_phi)
		tau_w_g_s2_j_phis := w._get_tau_w_g_s2_j_phis(cos_phi)
		return rho_w_g_s1f_j_phis + tau_w_g_s1_j_phis*tau_w_g_s2_j_phis*rho_w_g_s2f_j_phis/
			(1-math.Min(rho_w_g_s1b_j_phis*rho_w_g_s2f_j_phis, 0.9999))
	} else {
		panic(w._glass_type)
	}
}

/*
入射角Φに対する境界jの窓のガラス部分の室外側から1枚目の板ガラスの透過率を計算する。
Args:
	phis: 入射角Φ, rad
Returns:
	入射角Φに対する境界jの窓のガラス部分の室外側から1枚目の板ガラスの透過率, -
Notes:
	eq.14
*/
func (w *Window) _get_tau_w_g_s1_j_phis(cos_phi float64) float64 {
	return w._tau_w_g_s1_j * _get_tau_n_phi(cos_phi)
}

/*
入射角Φに対する境界jの窓のガラス部分の室外側から2枚目の板ガラスの透過率を計算する。
Args:
	phis: 入射角Φ, rad
Returns:
	入射角Φに対する境界jの窓のガラス部分の室外側から2枚目の板ガラスの透過率, -
Notes:
	この値は複層ガラスのみ計算される。
	eq.15
*/
func (w *Window) _get_tau_w_g_s2_j_phis(cos_phi float64) float64 {
	return w._tau_w_g_s2_j * _get_tau_n_phi(cos_phi)
}

/*
入射角Φに対する境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（正面側）を計算する。
Args:
	phis: 入射角Φ, rad
Returns:
	入射角Φに対する境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（正面側）, -
Notes:
	eq.16
*/
func (w *Window) _get_rho_w_g_s1f_j_phis(cos_phi float64) float64 {
	return w._rho_w_g_s1f_j + (1-w._rho_w_g_s1f_j)*_get_rho_n_phi(cos_phi)
}

/*
入射角Φに対する境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（背面側）を計算する。
Args:
	phis: 入射角Φ, rad
Returns:
	入射角Φに対する境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（背面側）, -
Notes:
	この値は複層ガラスのみ計算される。
	eq.17
*/
func (w *Window) _get_rho_w_g_s1b_j_phis(cos_phi float64) float64 {
	return w._rho_w_g_s1b_j + (1-w._rho_w_g_s1b_j)*_get_rho_n_phi(cos_phi)
}

/*
入射角Φに対する境界jの窓のガラス部分の室外側から2枚目の板ガラスの反射率（正面側）を計算する。
Args:
	phis: 入射角Φ, rad
Returns:
	入射角Φに対する境界jの窓のガラス部分の室外側から2枚目の板ガラスの反射率（正面側）, -
Notes:
	この値は複層ガラスのみ計算される。
	eq.18
*/
func (w *Window) _get_rho_w_g_s2f_j_phis(cos_phi float64) float64 {
	return w._rho_w_g_s2f_j + (1-w._rho_w_g_s2f_j)*_get_rho_n_phi(cos_phi)
}

/*
建具部分の熱損失係数（U値）を取得する。

Returns:
	境界 j の窓の建具部分の熱損失係数（U値）, W/m2K
*/
func (w *Window) u_w_f_j() float64 {
	return w._u_w_f_j
}

/*
窓の面積に対するグレージングの面積の比を取得する。

Returns:
	境界 j の窓の面積に対するグレージングの面積の比, -
*/
func (w *Window) r_a_w_g_j() float64 {
	return w._r_a_w_g_j
}

/*
窓のガラス部分の熱損失係数（U値）を取得する。

Returns:
	境界 j の窓のガラス部分の熱損失係数（U値）, W/m2K
*/
func (w *Window) u_w_g_j() float64 {
	return w._u_w_g_j
}

/*
窓のガラス部分の日射熱取得率を取得する。

Returns:
	境界 j の窓のガラス部分の日射熱取得率, -
*/
func (w *Window) eta_w_g_j() float64 {
	return w._eta_w_g_j
}

/*
窓のガラス部分の日射透過率を取得する。

Returns:
	境界 j の窓のガラス部分の日射透過率, -
*/
func (w *Window) tau_w_g_j() float64 {
	return w._tau_w_g_j
}

/*
規準化透過率を計算する。
    Args:
        phi: 入射角Φ, rad
    Returns:
        入射角Φに対する規準化透過率, -
    Notes:
        eq.19
*/
func _get_tau_n_phi(cos_phi float64) float64 {
	cos_phi_2 := cos_phi * cos_phi
	cos_phi_3 := cos_phi_2 * cos_phi
	cos_phi_4 := cos_phi_2 * cos_phi_2
	cos_phi_5 := cos_phi_2 * cos_phi_3
	return 2.552*cos_phi + 1.364*cos_phi_2 - 11.388*cos_phi_3 +
		13.617*cos_phi_4 - 5.146*cos_phi_5
}

/*
規準化反射率を計算する。
    Args:
        phi: 入射角Φ, rad
    Returns:
        入射角Φに対する規準化反射率, -
    Notes:
        eq.20
*/
func _get_rho_n_phi(cos_phi float64) float64 {
	cos_phi_2 := cos_phi * cos_phi
	cos_phi_3 := cos_phi_2 * cos_phi
	cos_phi_4 := cos_phi_2 * cos_phi_2
	cos_phi_5 := cos_phi_2 * cos_phi_3
	return 1.0 - 5.189*cos_phi + 12.392*cos_phi_2 - 16.593*cos_phi_3 +
		11.851*cos_phi_4 - 3.461*cos_phi_5
}

/*
境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（背面側）を計算する。

    Args:
        tau_w_g_s1_j: 境界jの窓のガラス部分の室外側から1枚目の板ガラスの透過率, -
        glass_type: 境界jの窓のガラス構成
    Returns:
        境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（背面側）, -
    Notes:
        複層ガラスの場合のみ定義される。
        eq.21
*/
func _get_rho_w_g_s1b_j(tau_w_g_s1_j float64, glass_type GlassType) float64 {
	if glass_type == GlassTypeSINGLE {
		return math.NaN()
	} else if glass_type == GlassTypeMULTIPLE {
		return 0.379 * (1 - tau_w_g_s1_j)
	} else {
		panic(glass_type)
	}
}

/*
境界jの窓のガラス部分の室外側から2枚目の板ガラスの透過率を計算する。

    Args:
        tau_w_g_s1_j: 境界jの窓のガラス部分の室外側から1枚目の板ガラスの透過率, -
        glass_type: 境界jの窓のガラス構成
    Returns:
        境界jの窓のガラス部分の室外側から2枚目の板ガラスの透過率, -
    Notes:
        複層ガラスの場合のみ定義される。
        eq.22
*/
func _get_tau_w_g_s2_j(tau_w_g_s1_j float64, glass_type GlassType) float64 {
	if glass_type == GlassTypeSINGLE {
		return math.NaN()
	} else if glass_type == GlassTypeMULTIPLE {
		return tau_w_g_s1_j
	} else {
		panic(glass_type)
	}
}

/*
境界jの窓のガラス部分の室外側から1枚目の板ガラスの透過率を計算する。
    Args:
        tau_w_g_j: 境界jの窓のガラス部分の日射透過率, -
        rho_w_g_s2f_j: 境界jの窓のガラス部分の室外側から2枚目の板ガラスの反射率（正面側）, -
        glass_type: 境界jの窓のガラス構成
    Returns:
        境界jの窓のガラス部分の室外側から1枚目の板ガラスの透過率, -
    Notes:
        eq.23
*/
func _get_tau_w_g_s1_j(tau_w_g_j float64, rho_w_g_s2f_j float64, glass_type GlassType) float64 {
	if glass_type == GlassTypeSINGLE {
		return tau_w_g_j
	} else if glass_type == GlassTypeMULTIPLE {
		temp := 0.379 * rho_w_g_s2f_j * tau_w_g_j
		temp_2 := temp * temp
		return (temp + math.Sqrt(temp_2-4.0*(0.379*rho_w_g_s2f_j-1)*tau_w_g_j)) / 2.0
	} else {
		panic(glass_type)
	}
}

/*
窓のガラス部分の日射透過率を計算する。
    Args:
        eta_w_g_j: 境界jの窓のガラス部分の日射熱取得率, -
        r_r_w_g_j: 境界jの窓のガラス部分の日射吸収量に対する室内側に放出される量の割合, -
        rho_w_g_s1f_j: 境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（正面側）, -
        rho_w_g_s2f_j: 境界jの窓のガラス部分の室外側から2枚目の板ガラスの反射率（正面側）, -
        glass_type: 境界jの窓のガラス構成
    Returns:
        境界jの窓のガラス部分の日射透過率, -
    Notes:
        eq.24
*/
func _get_tau_w_g_j(
	eta_w_g_j float64,
	r_r_w_g_j float64,
	rho_w_g_s1f_j float64,
	rho_w_g_s2f_j float64,
	glass_type GlassType,
) float64 {
	if glass_type == GlassTypeSINGLE {
		return (eta_w_g_j - (1.0-rho_w_g_s1f_j)*r_r_w_g_j) / (1.0 - r_r_w_g_j)
	} else if glass_type == GlassTypeMULTIPLE {
		return (eta_w_g_j - (1.0-rho_w_g_s1f_j)*r_r_w_g_j) / ((1.0 - r_r_w_g_j) - rho_w_g_s2f_j*r_r_w_g_j)
	} else {
		panic(glass_type)
	}
}

/*
境界jの窓のガラス部分の室外側から2枚目の板ガラスの反射率（正面側）を計算する。
    Args:
        glass_type: 境界 j の窓のガラス構成
    Returns:
        境界jの窓のガラス部分の室外側から2枚目の板ガラスの反射率（正面側）
    Notes:
        複層ガラスの場合のみ定義される。
        eq.25
*/
func _get_rho_w_g_s2f_j(glass_type GlassType) float64 {
	if glass_type == GlassTypeSINGLE {
		return math.NaN()
	} else if glass_type == GlassTypeMULTIPLE {
		return 0.077
	} else {
		panic(glass_type)
	}
}

/*
境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（正面側）を計算する。

    Args:
        r_r_w_g_j: 境界jの窓のガラス部分の日射吸収量に対する室内側に放出される量の割合, -
        eta_w_g_j: 境界jの窓のガラス部分の日射熱取得率, -
    Returns:
        境界jの窓のガラス部分の室外側から1枚目の板ガラスの反射率（正面側）, -
    Notes:
        eq.26, 27
*/
func _get_rho_w_g_s1f_j(r_r_w_g_j float64, eta_w_g_j float64) float64 {
	temp := 1.846 * r_r_w_g_j
	temp_2 := temp * temp
	t_j := (-temp + math.Sqrt(temp_2+4*(1-temp)*eta_w_g_j)) /
		(2 * (1 - temp))
	return 0.923*t_j*t_j - 1.846*t_j + 1
}

/*
境界jの窓のガラス部分の日射吸収量に対する室内側に放出される量の割合を計算する。
    Args:
        u_w_g_j: 境界jの窓のガラス部分の熱損失係数（U値）, W/m2K
        u_w_g_s_j: 境界jの窓のガラス部分の熱損失係数（夏期条件）, W/m2K
        r_w_o_w: 窓の室外側表面熱伝達抵抗（冬期条件）, m2K/W
        r_w_i_w: 窓の室内側表面熱伝達抵抗（冬期条件）, m2K/W
        r_w_o_s: 窓の室外側表面熱伝達抵抗（夏期条件）, m2K/W
        glass_type: 境界 j の窓のガラス構成
    Returns:
        境界jの窓のガラス部分の日射吸収量に対する室内側に放出される量の割合, -
    Notes:
        eq.28
*/
func _get_r_r_w_g_j(
	u_w_g_j float64,
	u_w_g_s_j float64,
	r_w_o_w float64,
	r_w_i_w float64,
	r_w_o_s float64,
	glass_type GlassType,
) float64 {
	if glass_type == GlassTypeSINGLE {
		return (1.0/2.0*(1/u_w_g_j-r_w_o_w-r_w_i_w) + r_w_o_s) * u_w_g_s_j
	} else if glass_type == GlassTypeMULTIPLE {
		// 複層ガラスにおける窓の中空層の熱伝達抵抗, m2K/W
		const r_w_air = 0.003
		return (1.0/4.0*(1/u_w_g_j-r_w_o_w-r_w_i_w-r_w_air) + r_w_o_s) * u_w_g_s_j
	} else {
		panic(glass_type)
	}
}

/*
境界jの窓のガラス部分の熱損失係数（夏期条件）を計算する。
    Args:
        u_w_g_j: 境界jの窓のガラス部分の熱損失係数（U値）, W/m2K
        r_w_o_w: 窓の室外側表面熱伝達抵抗（冬期条件）, m2K/W
        r_w_i_w: 窓の室内側表面熱伝達抵抗（冬期条件）, m2K/W
        r_w_o_s: 窓の室外側表面熱伝達抵抗（夏期条件）, m2K/W
        r_w_i_s: 窓の室内側表面熱伝達抵抗（夏期条件）, m2K/W
    Returns:
        境界jの窓のガラス部分の熱損失係数（夏期条件）, W/m2K
    Notes:
        eq.29
*/
func _get_u_w_g_s_j(u_w_g_j float64, r_w_o_w float64, r_w_i_w float64, r_w_o_s float64, r_w_i_s float64) float64 {
	return 1.0 / (1.0/u_w_g_j - r_w_o_w - r_w_i_w + r_w_o_s + r_w_i_s)
}

/*
窓の表面熱伝達抵抗を求める。
    Returns:
        窓の室外側表面熱伝達抵抗（冬季条件）, m2K/W
        窓の室内側表面熱伝達抵抗（冬季条件）, m2K/W
        窓の室外側表面熱伝達抵抗（夏季条件）, m2K/W
        窓の室内側表面熱伝達抵抗（夏季条件）, m2K/W
    Notes:
        table 2
*/
func _get_r_w() (float64, float64, float64, float64) {
	const r_w_o_w = 0.0415
	const r_w_i_w = 0.1228
	const r_w_o_s = 0.0756
	const r_w_i_s = 0.1317
	return r_w_o_w, r_w_i_w, r_w_o_s, r_w_i_s
}

/*
境界jの窓のガラス部分の日射熱取得率を取得する。
    Args:
        eta_w_j: 境界 j の窓の日射熱取得率
        r_a_w_g_j: 境界 j の窓の面積に対するグレージングの面積の比, -
    Returns:
        境界jの窓のガラス部分の日射熱取得率, -
    Notes:
        eq.30
*/
func _get_eta_w_g_j(eta_w_j float64, r_a_w_g_j float64) float64 {
	return eta_w_j / r_a_w_g_j
}

/*
境界jの窓のガラス部分の熱損失係数（U値）を取得する。
    Args:
        u_w_j: 境界jの窓の熱損失係数, W/m2K
        u_w_f_j: 境界jの窓の建具部分の熱損失係数（U値）, W/m2K
        r_a_w_g_j: 境界jの窓の面積に対するグレージングの面積の比, -
    Returns:
        境界jの窓のガラス部分の熱損失係数（U値）, W/m2K
    Notes:
        eq.31
*/
func _get_u_w_g_j(u_w_j float64, u_w_f_j float64, r_a_w_g_j float64) float64 {
	return (u_w_j - u_w_f_j*(1-r_a_w_g_j)) / r_a_w_g_j
}

/*
窓の面積に対するグレージングの面積の比が指定されていない場合に枠（フレーム）材質の種類に応じてデフォルト値を定める。

    Args:
        r_a_w_g_j: 境界jの窓の面積に対するグレージングの面積の比, -
        flame_type: 建具（フレーム）材質の種類
    Returns:
        境界jの窓の面積に対するグレージングの面積の比, -
    Notes:
        table 3
*/
func _get_r_a_w_g_j(r_a_w_g_j float64, flame_type FlameType) float64 {
	if math.IsNaN(r_a_w_g_j) {
		return map[FlameType]float64{
			FlameTypeRESIN:       0.72,
			FlameTypeWOOD:        0.72,
			FlameTypeALUMINUM:    0.8,
			FlameTypeMIXED_WOOD:  0.8,
			FlameTypeMIXED_RESIN: 0.8,
		}[flame_type]
	} else {
		return r_a_w_g_j
	}
}

/*
建具部分の熱損失係数（U値）を取得する。
    Args:
        flame_type: 建具（フレーム）材質の種類
    Returns:
        境界jの窓の建具部分の熱損失係数（U値）, W/m2K
    Notes:
        table 4
*/
func _get_u_w_f_j(flame_type FlameType) float64 {
	return map[FlameType]float64{
		FlameTypeRESIN:       2.2,
		FlameTypeWOOD:        2.2,
		FlameTypeALUMINUM:    6.6,
		FlameTypeMIXED_WOOD:  4.7,
		FlameTypeMIXED_RESIN: 4.7,
	}[flame_type]
}
