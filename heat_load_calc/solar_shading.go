package heat_load_calc

import (
	"math"
)

// 日よけ計算に関するインタフェイス
type SolarShading interface {
	/*
		直達日射に対する日よけの影面積比率を計算する。

		Args:
			h_sun_n: 太陽高度, rad, [N+1]
			a_sun_n: 太陽方位角, rad, [N+1]

		Returns:
			直達日射に対する日除けの影面積比率, [N+1]
	*/
	get_f_ss_dn_j_ns(h_sun_n []float64, a_sun_n []float64) []float64

	/*
		天空放射に対する日よけの影面積比率を計算する。

		Returns:
			天空放射に対する日除けの影面積比率, -
	*/
	get_f_ss_sky_j() float64

	/*
		地面反射に対する日よけの影面積比率を計算する。
		Returns:
			地面反射に対する日よけの影面積比率
	*/
	get_f_ss_ref_j() float64
}

/*
入力ファイルの辞書の"solar_shading_part"を読み込む。

Args:

	ssp_dict: 日除けの仕様に関する辞書
	direction: 方位

Returns:

	SolarShadingPart クラス
*/
func NewSolarShading(ssp_dict *SolarShadingPartJson, direction Direction) SolarShading {
	if ssp_dict.Existence {

		input_method := ssp_dict.InputMethod

		if direction == DirectionTop || direction == DirectionBottom {
			panic("方位が「上方」「下方」の場合に日除けを定義することはできません。")
		}

		// 境界ｊの傾斜面の方位角, rad
		alpha_w_j := direction.alpha_w_j()

		if input_method == "simple" {
			return NewSolarShadingSimple(
				alpha_w_j,
				ssp_dict.Depth,
				ssp_dict.D_h,
				ssp_dict.D_e,
			)
		} else if input_method == "detail" {
			return NewSolarShadingDetail(
				alpha_w_j,
				ssp_dict.X1,
				ssp_dict.X2,
				ssp_dict.X3,
				ssp_dict.Y1,
				ssp_dict.Y2,
				ssp_dict.Y3,
				ssp_dict.Z_x_pls,
				ssp_dict.Z_x_mns,
				ssp_dict.Z_y_pls,
				ssp_dict.Z_y_mns,
			)
		} else {
			panic(input_method)
		}
	} else {
		return NewSolarShadingNot()
	}
}

type SolarShadingSimple struct {
	_alpha_w_j float64 // 開口部 j の方位角, rad
	_l_z_j     float64 // 開口部 j のひさしの出幅, m
	_l_y_h_j   float64 // 開口部 j の高さ, m
	_l_y_e_j   float64 // 開口部 j の上端から日よけまでの長さ, m
}

/*
Args:

	alpha_w_j: 開口部 j の方位角, rad
	l_z_j: 開口部 j のひさしの出幅, m
	l_y_h_j: 開口部 j の高さ, m
	l_y_e_j: 開口部 j の上端から日よけまでの長さ, m
*/
func NewSolarShadingSimple(alpha_w_j float64, l_z_j float64, l_y_h_j float64, l_y_e_j float64) *SolarShadingSimple {
	return &SolarShadingSimple{
		_alpha_w_j: alpha_w_j,
		_l_z_j:     l_z_j,
		_l_y_h_j:   l_y_h_j,
		_l_y_e_j:   l_y_e_j,
	}
}

/*
無限に長い庇がある場合の直達日射に対する日よけの日影面積比率を計算する。

Args:

	h_sun_n: ステップ n における太陽高度, rad, [N+1]
	a_sun_n: ステップ n における太陽方位角, rad, [N+1]

Returns:

	ステップ n における直達日射に対する日除けの日影面積比率, [N+1]

Notes:

	日射が壁にあたらない場合は日影そのものができない。
	その場合は値として 0.0 を返す。
*/
func (self *SolarShadingSimple) get_f_ss_dn_j_ns(h_sun_n, a_sun_n []float64) []float64 {
	n := len(h_sun_n)
	f_ss_d_j_ns := make([]float64, n)

	for i := 0; i < n; i++ {
		h_s_n := math.Max(h_sun_n[i], 0.0)
		var a_s_n float64
		if h_sun_n[i] > 0.0 {
			a_s_n = a_sun_n[i]
		} else {
			a_s_n = 0.0
		}

		cos_a := math.Cos(a_s_n - self._alpha_w_j)
		if cos_a <= 0 {
			cos_a = 1.0
		}

		// ステップ n における境界 j に対する太陽のプロファイル角の正弦, -
		tan_phi_j_n := math.Tan(h_s_n) / cos_a

		// ステップ n における開口部 j に影がかかる長さ（窓上端から下方への長さ）, m
		l_ss_d_y_j_n := self._l_z_j*tan_phi_j_n - self._l_y_e_j

		// 日影面積率の計算 式(79)
		//   マイナスの場合（日陰が窓上端にかからない場合）は 0.0 とする。
		//   1.0を超える場合（日陰が窓下端よりも下になる場合）は 1.0 とする。
		f_ss_d_j_ns[i] = math.Min(math.Max(l_ss_d_y_j_n/self._l_y_h_j, 0.0), 1.0)

		// 日が出ていないときは 0.0 とする。
		if h_sun_n[i] <= 0.0 {
			f_ss_d_j_ns[i] = 0.0
		}

		// 太陽位置が背面にある場合は 0.0 とする。
		if math.Cos(a_s_n-self._alpha_w_j) <= 0.0 {
			f_ss_d_j_ns[i] = 0.0
		}
	}

	return f_ss_d_j_ns
}

/*
無限に長い庇がある場合の天空放射に対する日よけの影面積比率を計算する。

Returns:

	天空放射に対する日除けの影面積比率, -
*/
func (self *SolarShadingSimple) get_f_ss_sky_j() float64 {
	temp := self._l_y_e_j + self._l_y_h_j
	return (((self._l_y_e_j + self._l_y_h_j) + math.Sqrt(self._l_y_e_j*self._l_y_e_j+self._l_z_j*self._l_z_j)) -
		(self._l_y_e_j + math.Sqrt(temp*temp+self._l_z_j*self._l_z_j))) /
		(2.0 * self._l_y_h_j)
}

/*
地面反射に対する日よけの影面積比率を計算する。
Returns:

	地面反射に対する日よけの影面積比率
*/
func (self *SolarShadingSimple) get_f_ss_ref_j() float64 {
	return 0.0
}

type SolarShadingDetail struct {
	w_alpha float64
	x1      float64
	x2      float64
	x3      float64
	y1      float64
	y2      float64
	y3      float64
	z_x_pls float64
	z_x_mns float64
	z_y_pls float64
	z_y_mns float64
}

func NewSolarShadingDetail(alpha_w_j, x1, x2, x3, y1, y2, y3, z_x_pls, z_x_mns, z_y_pls, z_y_mns float64) *SolarShadingDetail {
	return &SolarShadingDetail{
		w_alpha: alpha_w_j,
		x1:      x1,
		x2:      x2,
		x3:      x3,
		y1:      y1,
		y2:      y2,
		y3:      y3,
		z_x_pls: z_x_pls,
		z_x_mns: z_x_mns,
		z_y_pls: z_y_pls,
		z_y_mns: z_y_mns,
	}
}

/*
直達日射に対する日よけの影面積比率を計算する。

Args:

	h_sun_n: 太陽高度, rad, [N+1]
	a_sun_n: 太陽方位角, rad, [N+1]

Returns:

	直達日射に対する日除けの影面積比率, [N+1]
*/
func (self *SolarShadingDetail) get_f_ss_dn_j_ns(h_sun_n, a_sun_n []float64) []float64 {
	panic("Not Implemented Yet")
}

/*
天空放射に対する日よけの影面積比率を計算する。

Returns:

	天空放射に対する日除けの影面積比率, -
*/
func (self *SolarShadingDetail) get_f_ss_sky_j() float64 {
	panic("Not Implemented Yet")
}

/*
地面反射に対する日よけの影面積比率を計算する。
Returns:

	地面反射に対する日よけの影面積比率
*/
func (self *SolarShadingDetail) get_f_ss_ref_j() float64 {
	panic("Not Implemented Yet")
}

type SolarShadingNot struct {
}

func NewSolarShadingNot() *SolarShadingNot {
	return &SolarShadingNot{}
}

/*
直達日射に対する日よけの影面積比率を計算する。

Args:

	h_sun_n: 太陽高度, rad, [N+1]
	a_sun_n: 太陽方位角, rad, [N+1]

Returns:

	直達日射に対する日除けの影面積比率, [N+1]
*/
func (self *SolarShadingNot) get_f_ss_dn_j_ns(h_sun_n []float64, a_sun_n []float64) []float64 {
	return make([]float64, len(h_sun_n))
}

/*
天空放射に対する日よけの影面積比率を計算する。

Returns:

	天空放射に対する日除けの影面積比率, -
*/
func (self *SolarShadingNot) get_f_ss_sky_j() float64 {
	return 0.0
}

/*
地面反射に対する日よけの影面積比率を計算する。
Returns:

	地面反射に対する日よけの影面積比率
*/
func (self *SolarShadingNot) get_f_ss_ref_j() float64 {
	return 0.0
}
