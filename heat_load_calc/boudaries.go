package main

import (
	"errors"
	"fmt"
	"math"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

const InvalidRearSurfaceBoundaryID = 9999

// 境界の種類
type BoundaryType int

const (
	BoundaryTypeInternal                BoundaryType = iota // "internal": 間仕切り
	BoundaryTypeExternalGeneralPart                         // "external_general_part": 外皮_一般部位
	BoundaryTypeExternalTransparentPart                     // "external_transparent_part": 外皮_透明な開口部
	BoundaryTypeExternalOpaquePart                          // "external_opaque_part": 外皮_不透明な開口部
	BoundaryTypeGround                                      // "ground": 地盤
)

func (b BoundaryType) String() string {
	return [...]string{"internal", "external_general_part", "external_transparent_part", "external_opaque_part", "ground"}[b]
}

func BoundaryFromString(str string) (BoundaryType, error) {
	switch str {
	case "internal":
		return BoundaryTypeInternal, nil
	case "external_general_part":
		return BoundaryTypeExternalGeneralPart, nil
	case "external_transparent_part":
		return BoundaryTypeExternalTransparentPart, nil
	case "external_opaque_part":
		return BoundaryTypeExternalOpaquePart, nil
	case "ground":
		return BoundaryTypeGround, nil
	default:
		return 0, errors.New("invalid boundary type")
	}
}

type Boundary struct {

	// ID
	id int

	// 名称
	name string

	// 副名称
	sub_name string

	// 接する室のID
	connected_room_id int

	// 境界の種類
	boundary_type BoundaryType

	// 面積, m2
	area float64

	// 温度差係数
	h_td float64

	// 裏側表面の境界ID
	// internal_wall の場合のみ定義される。
	//rear_surface_boundary_id int

	// 床か否か
	is_floor bool

	// 室内侵入日射吸収の有無
	is_solar_absorbed_inside bool

	// 室外側の日射の有無
	// True 当たる
	// False 当たらない
	// 境界の種類が"external_general_part", "external_transparent_part", "external_opaque_part"の場合に定義される。
	is_sun_striked_outside bool

	// 室内側表面対流熱伝達率, W/m2K
	h_s_c float64

	// 室内側表面放射熱伝達率, W/m2K
	h_s_r float64

	// 相当外気温度, ℃, [8760 * 4]
	theta_o_sol *mat.VecDense

	// 透過日射熱取得, W, [8760*4]
	q_trs_sol []float64

	// 応答係数データクラス
	rf *ResponseFactor

	// 裏面温度に他の境界 j の等価室温が与える影響, [j, j]
	k_ei_js_j []float64

	// 裏面温度に室 i の室温が与える影響, [i]
	k_s_r_j_is []float64

	// 計算で使用する熱貫流率, W/m2K
	simulation_u_value float64
}

type Boundaries struct {
	bss                   []*Boundary // 境界
	p_is_js               *mat.Dense  // 室iと境界jの関係を表す係数（境界jから室iへの変換）, [i, j]
	q_trs_sol_is_ns       mat.Matrix  // ステップ n の室 i における窓の透過日射熱取得, W, [n]
	n_b                   int         // 境界の数
	n_ground              int         // 地盤の数
	id_bs_js              []int       // ID
	name_js               []string    // 名前, [j, 1]
	sub_name_js           []string    // 名前2, [j, 1]
	p_is_is               mat.Matrix  // 室iと境界jの関係を表す係数（境界jから室iへの変換）
	p_js_is               mat.Matrix  // 室iと境界jの関係を表す係数（室iから境界jへの変換）
	is_floor_js           []bool      // 床かどうか, [j, 1]
	is_ground_js          []bool      // 地盤かどうか, [j, 1]
	k_ei_js_js            mat.Matrix  // 境界jの裏面温度に他の境界の等価室温が与える影響, [j, j]
	k_eo_js               mat.Vector  // 温度差係数
	k_s_r_js              mat.Matrix  // 境界 j の裏面温度に室温が与える影響, [j, i]
	p_s_sol_abs_js        mat.Vector  // 境界jの日射吸収の有無, [j, 1]
	h_s_r_js              mat.Vector  // 境界jの室内側表面放射熱伝達率, W/m2K, [j, 1]
	h_s_c_js              mat.Vector  // 境界jの室内側表面対流熱伝達率, W/m2K, [j, 1]
	simulation_u_value_js mat.Vector  // 境界jの計算で使用する熱貫流率, W/m2K, [j, 1]
	a_s_js                mat.Vector  // 境界jの面積, m2, [j, 1]
	phi_a0_js             mat.Vector  // 境界jの吸熱応答係数の初項, m2K/W, [j, 1]
	phi_a1_js_ms          [][]float64 // 境界jの項別公比法における項mの吸熱応答係数の第一項 , m2K/W, [j, 12]
	phi_t0_js             mat.Vector  // 境界jの貫流応答係数の初項, [j, 1]
	phi_t1_js_ms          [][]float64 // 境界jの項別公比法における項mの貫流応答係数の第一項, [j, 12]
	r_js_ms               [][]float64 // 境界jの項別公比法における項mの公比, [j, 12]
	theta_o_eqv_js_ns     *mat.Dense  // ステップ n の境界 j における相当外気温度, ℃, [j, n+1]
	q_trs_sol_js_ns       mat.Matrix  // ステップ n の室 i における窓の透過日射熱取得, W, [n]
}

/*
Args

	id_rm_is 室のID, [i, 1]
	bs_list 境界に関する辞書
	w Weather クラス

Notes

	本来であれば Boundaries クラスにおいて境界に関する入力用辞書から読み込みを境界個別に行う。
	しかし、室内側表面放射熱伝達は室内側の形態係数によって値が決まり、ある室に接する境界の面積の組み合わせで決定されるため、
	境界個別に値を決めることはできない。（すべての境界の情報が必要である。）
	一方で、境界の集約を行うためには、応答係数を Boundary クラス生成時に求める必要があり、
	さらに応答係数の計算には裏面の表面放射・対流熱伝達率の値が必要となるため、
	Boundary クラスを生成する前に、予め室内側表面放射・対流熱伝達率を計算しておき、
	Boundary クラスを生成する時に必要な情報としておく。
*/
func NewBoundaries(id_rm_is []int, bs_list []BoudaryJson, w *Weather) *Boundaries {

	_areas := mat.NewVecDense(len(bs_list), nil)
	_connected_room_ids := make([]int, len(bs_list))
	h_c_js := make([]float64, len(bs_list))
	for j := 0; j < len(bs_list); j++ {
		b := bs_list[j]
		_areas.SetVec(j, b.Area)
		_connected_room_ids[j] = int(b.ConnectedRoomId)

		// 境界jの室内側表面対流熱伝達率, W/m2K, [J, 1]
		h_c_js[j] = b.H_c
	}

	// 境界jの室内側表面放射熱伝達率, W/m2K, [J, 1]
	h_s_r_js := get_h_s_r_js(id_rm_is, _areas, _connected_room_ids)

	// 室の数
	n_rm := len(id_rm_is)

	// 境界 j, [J]
	bss := make([]*Boundary, len(bs_list))
	for j := 0; j < len(bs_list); j++ {
		b := &bs_list[j]
		//log.Printf("境界%d - %s (%s)", j, b["name"].(string), b["boundary_type"].(string))
		bss[j] = _get_boundary(b, h_c_js, h_s_r_js, w, n_rm)
	}

	// 室iと境界jの関係を表す係数（境界jから室iへの変換）, [i, j]
	p_is_js := _get_p_is_js(n_rm, bss)

	// ステップ n の室 i における窓の透過日射熱取得, W, [n]
	q_trs_sol_is_ns := get_q_trs_sol_is_ns(n_rm, bss)

	bs := &Boundaries{
		bss:             bss,
		p_is_js:         p_is_js,
		q_trs_sol_is_ns: q_trs_sol_is_ns,
	}
	self := bs

	// 境界の数
	bs.n_b = len(bs.bss)

	// 地盤の数
	sum := 0
	for _, bs := range self.bss {
		if bs.boundary_type == BoundaryTypeGround {
			sum++
		}
	}
	bs.n_ground = sum

	// ID
	ids := make([]int, len(self.bss))
	for i, bs := range self.bss {
		ids[i] = bs.id
	}
	bs.id_bs_js = ids

	// 名前, [j, 1]
	names := make([]string, len(self.bss))
	for i, bs := range self.bss {
		names[i] = bs.name
	}
	bs.name_js = names

	// 名前2, [j, 1]
	sub_names := make([]string, len(self.bss))
	for i, bs := range self.bss {
		sub_names[i] = bs.sub_name
	}
	bs.sub_name_js = sub_names

	// 室iと境界jの関係を表す係数（境界jから室iへの変換）
	bs.p_is_js = p_is_js

	// 室iと境界jの関係を表す係数（室iから境界jへの変換）
	bs.p_js_is = p_is_js.T()

	// 床かどうか, [j, 1]
	is_floor_js := make([]bool, len(self.bss))
	for i, bs := range self.bss {
		is_floor_js[i] = bs.is_floor
	}
	bs.is_floor_js = is_floor_js

	// 地盤かどうか, [j, 1]
	is_grounds := make([]bool, len(self.bss))
	for i, bs := range self.bss {
		is_grounds[i] = (bs.boundary_type == BoundaryTypeGround)
	}
	bs.is_ground_js = is_grounds

	// 境界jの裏面温度に他の境界の等価室温が与える影響, [j, j]
	n := len(self.bss[0].k_ei_js_j)
	k_ei_js_j := mat.NewDense(self.n_b, n, nil)
	for i, bs := range self.bss {
		k_ei_js_j.SetRow(i, bs.k_ei_js_j)
	}
	bs.k_ei_js_js = k_ei_js_j

	// 温度差係数
	h_td := make([]float64, len(self.bss))
	for i, bs := range self.bss {
		h_td[i] = bs.h_td
	}
	bs.k_eo_js = mat.NewVecDense(len(h_td), h_td)

	// 境界 j の裏面温度に室温が与える影響, [j, i]
	k_s_r_j_is := mat.NewDense(self.n_b, len(self.bss[0].k_s_r_j_is), nil)
	for i, bs := range self.bss {
		k_s_r_j_is.SetRow(i, bs.k_s_r_j_is)
	}
	bs.k_s_r_js = k_s_r_j_is

	// 境界jの日射吸収の有無, [j, 1]
	is_solar_absorbed_inside := make([]float64, len(self.bss))
	for i, bs := range self.bss {
		if bs.is_solar_absorbed_inside {
			is_solar_absorbed_inside[i] = 1.0
		} else {
			is_solar_absorbed_inside[i] = 0.0
		}
	}
	bs.p_s_sol_abs_js = mat.NewVecDense(len(is_solar_absorbed_inside), is_solar_absorbed_inside)

	// 境界jの室内側表面放射熱伝達率, W/m2K, [j, 1]
	h_s_r := make([]float64, len(self.bss))
	for i, bs := range self.bss {
		h_s_r[i] = bs.h_s_r
	}
	bs.h_s_r_js = mat.NewVecDense(len(h_s_r), h_s_r)

	// 境界jの室内側表面対流熱伝達率, W/m2K, [j, 1]
	h_s_c := make([]float64, len(self.bss))
	for i, bs := range self.bss {
		h_s_c[i] = bs.h_s_c
	}
	bs.h_s_c_js = mat.NewVecDense(len(h_s_c), h_s_c)

	// 境界jの室内側表面対流熱伝達率, W/m2K, [j, 1]
	u := make([]float64, len(self.bss))
	for i, bs := range self.bss {
		u[i] = bs.simulation_u_value
	}
	bs.simulation_u_value_js = mat.NewVecDense(len(u), u)

	// 境界jの面積, m2, [j, 1]
	area := make([]float64, len(self.bss))
	for i, bs := range self.bss {
		area[i] = bs.area
	}
	bs.a_s_js = mat.NewVecDense(len(area), area)

	// 境界jの吸熱応答係数の初項, m2K/W, [j, 1]
	rfa0 := make([]float64, len(self.bss))
	for i, bs := range self.bss {
		rfa0[i] = bs.rf.rfa0
	}
	bs.phi_a0_js = mat.NewVecDense(len(rfa0), rfa0)

	// 境界jの項別公比法における項mの吸熱応答係数の第一項 , m2K/W, [j, 12]
	rfa1 := make([][]float64, self.n_b)
	for i := 0; i < self.n_b; i++ {
		rfa1[i] = make([]float64, 12)
		copy(rfa1[i], self.bss[i].rf.rfa1)
	}
	bs.phi_a1_js_ms = rfa1

	// 境界jの貫流応答係数の初項, [j, 1]
	rft0 := make([]float64, len(self.bss))
	for i, bs := range self.bss {
		rft0[i] = bs.rf.rft0
	}
	bs.phi_t0_js = mat.NewVecDense(len(rft0), rft0)

	// 境界jの項別公比法における項mの貫流応答係数の第一項, [j, 12]
	rft1 := make([][]float64, self.n_b)
	for i := 0; i < self.n_b; i++ {
		rft1[i] = make([]float64, 12)
		copy(rft1[i], self.bss[i].rf.rft1)
	}
	bs.phi_t1_js_ms = rft1

	// 境界jの項別公比法における項mの公比, [j, 12]
	row := make([][]float64, self.n_b)
	for i := 0; i < self.n_b; i++ {
		row[i] = make([]float64, 12)
		copy(row[i], self.bss[i].rf.row)
	}
	bs.r_js_ms = row

	// ステップ n の境界 j における相当外気温度, ℃, [j, n+1]
	theta := mat.NewDense(self.n_b, self.bss[0].theta_o_sol.Len(), nil)
	for i := 0; i < self.n_b; i++ {
		theta.SetRow(i, self.bss[i].theta_o_sol.RawVector().Data)
	}
	bs.theta_o_eqv_js_ns = theta

	// ステップ n の室 i における窓の透過日射熱取得, W, [n]
	bs.q_trs_sol_is_ns = q_trs_sol_is_ns

	return bs
}

/*
Args

	b Boundary　の辞書
	h_c_js 境界 j の室内側表面対流熱伝達率, W/m2K, [J, 1]
	h_s_r_js 境界 j の室内側表面放射熱伝達率, W/m2K, [J, 1]
	w Weather クラス
	n_rm 室の数

Returns

	Boundary クラス
*/
func _get_boundary(
	b *BoudaryJson,
	h_c_js []float64,
	h_s_r_js []float64,
	w *Weather,
	n_rm int,
) *Boundary {

	// ID
	// TODO ID が0始まりで1ずつ増え、一意であることのチェックを行うコードを追記する。
	boundary_id := int(b.Id)

	// 名前
	name := b.Name

	// 副名称
	sub_name := b.SubName

	// 接する室のID
	// TODO 指定された room_id が存在するかどうかをチェックするコードを追記する。
	connected_room_id := int(b.ConnectedRoomId)

	// 境界の種類
	boundary_type, err := BoundaryFromString(b.BoundaryType)
	if err != nil {
		panic(err)
	}

	// 面積, m2
	area := b.Area
	if area <= 0.0 {
		panic(fmt.Sprintf("境界(ID=%d)の面積で0以下の値が指定されました。", boundary_id))
	}

	// 温度差係数
	// 境界の種類が"external_general_part", "external_transparent_part", "external_opaque_part"の場合に定義される。
	// TODO 下記、BoundaryType.Ground が指定されているのは間違い？　要チェック。
	var h_td float64
	if boundary_type == BoundaryTypeExternalGeneralPart ||
		boundary_type == BoundaryTypeExternalTransparentPart ||
		boundary_type == BoundaryTypeExternalOpaquePart {
		h_td = b.TempDifCoef
	} else if boundary_type == BoundaryTypeGround {
		h_td = 1.0
	} else {
		h_td = 0.0
	}

	if h_td > 1.0 {
		panic(fmt.Sprintf("境界(ID=%d)の温度差係数で1.0を超える値が指定されました。", boundary_id))
	}
	if h_td < 0.0 {
		panic(fmt.Sprintf("境界(ID=%d)の温度差係数で0.0を下回る値が指定されました。", boundary_id))
	}

	// TODO 指定された rear_surface_boundary_id の rear_surface_boundary_id が自分自信のIDかどうかのチェックが必要かもしれない。
	var rear_surface_boundary_id int
	if boundary_type == BoundaryTypeInternal {
		rear_surface_boundary_id = int(b.RearSurfaceBoundaryId)
	} else {
		// 便宜的にありえない値を入れている
		rear_surface_boundary_id = InvalidRearSurfaceBoundaryID
	}

	// 室内侵入日射吸収の有無 (True吸収する/False吸収しない)
	is_solar_absorbed_inside := b.IsSolarAbsorbedInside

	// 床か否か(True床/False床以外)
	is_floor := b.IsFloor

	// 室内側表面対流熱伝達率, W/m2K
	h_s_c := b.H_c

	// 室内側表面放射熱伝達率, W/m2K
	h_s_r := h_s_r_js[boundary_id]

	// 日射の有無 (True当たる/False 当たらない)
	// 境界の種類が"external_general_part", "external_transparent_part", "external_opaque_part"の場合に定義される。
	var is_sun_striked_outside bool
	if boundary_type == BoundaryTypeExternalGeneralPart ||
		boundary_type == BoundaryTypeExternalTransparentPart ||
		boundary_type == BoundaryTypeExternalOpaquePart {
		is_sun_striked_outside = b.IsSunStrikedOtside
	} else {
		is_sun_striked_outside = false
	}

	var simulation_u_value float64
	var theta_o_eqv_j_ns *mat.VecDense
	var rf *ResponseFactor
	var q_trs_sol []float64

	if boundary_type == BoundaryTypeInternal {

		// 相当外気温度, ℃
		theta_o_eqv_j_ns = get_theta_o_eqv_j_ns_for_internal(w)

		// 透過日射量, W, [N+1]
		q_trs_sol = get_q_trs_sol_j_ns_for_not(w)

		layers := b.Layers
		cs := make([]float64, len(layers))
		rs := make([]float64, len(layers))
		rs_sum := 0.0
		for i := range layers {
			layer := &layers[i]
			cs[i] = _read_cs_j_l(layer)
			rs[i] = _read_rs_j_l(layer)
			rs_sum += rs[i]
		}

		rear_h_c := h_c_js[int(b.RearSurfaceBoundaryId)]
		rear_h_r := h_s_r_js[int(b.RearSurfaceBoundaryId)]

		r_o := 1.0 / (rear_h_c + rear_h_r)

		// 応答係数
		rf = create_for_unsteady_not_ground(cs, rs, r_o)

		// U値
		simulation_u_value = 1.0 / (1.0/(h_s_c+h_s_r) + rs_sum + 1.0/(rear_h_c+rear_h_r))

	} else if boundary_type == BoundaryTypeExternalGeneralPart {

		if is_sun_striked_outside {

			// 方位
			drct_j := DirectionFromString(b.Direction)

			// 日除け
			ssp_j := NewSolarShading(&b.SolarShadingPart, drct_j)

			// 境界jの室外側日射吸収率, -
			a_s_j := _read_a_s_j(b, boundary_id)

			r_s_o_j := _read_r_s_o_j(b, boundary_id)

			eps_r_o_j := _read_eps_r_o_j(b, boundary_id)

			// 相当外気温度, ℃
			theta_o_eqv_j_ns = get_theta_o_eqv_j_ns_for_external_general_part_and_external_opaque_part(
				drct_j, a_s_j, eps_r_o_j, r_s_o_j, ssp_j, w,
			)

		} else {

			// 相当外気温度, ℃
			theta_o_eqv_j_ns = get_theta_o_eqv_j_ns_for_external_not_sun_striked(w)
		}

		// 透過日射量, W, [N+1]
		q_trs_sol = get_q_trs_sol_j_ns_for_not(w)

		layers := b.Layers
		cs := make([]float64, len(layers))
		rs := make([]float64, len(layers))
		rs_sum := 0.0
		for i := range layers {
			layer := &layers[i]
			cs[i] = _read_cs_j_l(layer)
			rs[i] = _read_rs_j_l(layer)
			rs_sum += rs[i]
		}

		r_o := b.OutsideHeatTransferResistance

		rf = create_for_unsteady_not_ground(cs, rs, r_o)

		// U値
		simulation_u_value = 1.0 / (1.0/(h_s_c+h_s_r) + rs_sum + r_o)

	} else if boundary_type == BoundaryTypeExternalTransparentPart {

		u_value_j := _read_u_nominal_j(b, boundary_id)

		if is_sun_striked_outside {

			// 方位
			drct_j := DirectionFromString(b.Direction)

			// 日除け
			ssp_j := NewSolarShading(&b.SolarShadingPart, drct_j)

			// 日射熱取得率
			eta_value := b.EtaValue
			if eta_value <= 0.0 {
				panic(fmt.Sprintf("境界(ID=%d)の日射熱取得率で0.0以下の値が指定されました。", boundary_id))
			}

			// 開口部の面積に対するグレージングの面積の比率
			glass_area_ratio := b.GlassAreaRatio
			if glass_area_ratio < 0.0 {
				panic(fmt.Sprintf("境界(ID=%d)の開口部の面積に対するグレージング面積の比率で0.0未満の値が指定されました。", boundary_id))
			}
			if glass_area_ratio > 1.0 {
				panic(fmt.Sprintf("境界(ID=%d)の開口部の面積に対するグレージング面積の比率で1.0より大の値が指定されました。", boundary_id))
			}

			// グレージングの種類
			glazing_type := GlassTypeFromString(b.IncidentAngleCharacteristics)

			wdw_j := NewWindow(
				u_value_j, eta_value, glazing_type, glass_area_ratio, FlameTypeMIXED_WOOD,
			)

			r_s_o_j := _read_r_s_o_j(b, boundary_id)

			eps_r_o_j := _read_eps_r_o_j(b, boundary_id)

			// 相当外気温度, ℃
			theta_o_eqv_j_ns = get_theta_o_eqv_j_ns_for_external_transparent_part(
				drct_j, eps_r_o_j, r_s_o_j, u_value_j, ssp_j, wdw_j, w,
			)

			// 透過日射量, W, [N+1]
			q_trs_sol = get_q_trs_sol_j_ns_for_transparent_sun_striked(
				drct_j, area, ssp_j, wdw_j, w,
			)

		} else {

			// 相当外気温度, ℃
			theta_o_eqv_j_ns = get_theta_o_eqv_j_ns_for_external_not_sun_striked(w)

			// 透過日射量, W, [N+1]
			q_trs_sol = get_q_trs_sol_j_ns_for_not(w)
		}

		// 室内側熱伝達抵抗, m2K/W
		r_i_nominal := _read_r_i_nominal(b, boundary_id)

		// 応答係数
		rf = create_for_steady(u_value_j, r_i_nominal)

		u_value_nominal := b.UValue

		simulation_u_value = 1.0 / (1.0/u_value_nominal - r_i_nominal + 1.0/(h_s_c+h_s_r))

	} else if boundary_type == BoundaryTypeExternalOpaquePart {

		if is_sun_striked_outside {

			// 方位
			drct_j := DirectionFromString(b.Direction)

			// 日除け
			ssp_j := NewSolarShading(&b.SolarShadingPart, drct_j)

			// 室外側日射吸収率
			a_s_j := _read_a_s_j(b, boundary_id)

			r_s_o_j := _read_r_s_o_j(b, boundary_id)

			eps_r_o_j := _read_eps_r_o_j(b, boundary_id)

			// 相当外気温度, ℃
			theta_o_eqv_j_ns = get_theta_o_eqv_j_ns_for_external_general_part_and_external_opaque_part(
				drct_j, a_s_j, eps_r_o_j, r_s_o_j, ssp_j, w,
			)

		} else {

			// 相当外気温度, ℃
			theta_o_eqv_j_ns = get_theta_o_eqv_j_ns_for_external_not_sun_striked(w)

		}

		// 透過日射量, W, [N+1]
		q_trs_sol = get_q_trs_sol_j_ns_for_not(w)

		// 室内側熱伝達抵抗, m2K/W
		r_i_nominal := _read_r_i_nominal(b, boundary_id)

		u_value_j := _read_u_nominal_j(b, boundary_id)

		rf = create_for_steady(u_value_j, r_i_nominal)

		u_value_nominal := b.UValue
		simulation_u_value = 1.0 / (1.0/u_value_nominal - r_i_nominal + 1.0/(h_s_c+h_s_r))

	} else if boundary_type == BoundaryTypeGround {

		// 相当外気温度, ℃
		theta_o_eqv_j_ns = get_theta_o_eqv_j_ns_for_ground(w)

		// 透過日射量, W, [N+1]
		q_trs_sol = get_q_trs_sol_j_ns_for_not(w)

		layers := b.Layers
		cs := make([]float64, len(layers))
		rs := make([]float64, len(layers))
		rs_sum := 0.0
		for i := range layers {
			layer := &layers[i]
			cs[i] = _read_cs_j_l(layer)
			rs[i] = _read_rs_j_l(layer)
			rs_sum += rs[i]
		}

		// 応答係数
		rf = create_for_unsteady_ground(cs, rs)

		// U値
		simulation_u_value = 1.0 / (1.0/(h_s_c+h_s_r) + rs_sum)

	} else {
		panic("invalid boundary type")
	}

	// Boundary の数
	n_b := len(h_c_js)

	k_ei_js_j := make([]float64, n_b)

	if boundary_type == BoundaryTypeExternalOpaquePart ||
		boundary_type == BoundaryTypeExternalTransparentPart ||
		boundary_type == BoundaryTypeExternalGeneralPart {
		//pass
	} else if boundary_type == BoundaryTypeInternal {

		// 室内壁の場合にk_ei_jsを登録する。
		k_ei_js_j[int(rear_surface_boundary_id)] = 1.0

	} else {

		// 外皮に面していない場合、室内壁ではない場合（地盤の場合が該当）は、k_ei_js に操作は行わない。
		//pass
	}

	k_s_r_j_is := make([]float64, n_rm)

	if boundary_type == BoundaryTypeExternalOpaquePart ||
		boundary_type == BoundaryTypeExternalTransparentPart ||
		boundary_type == BoundaryTypeExternalGeneralPart {

		h := h_td

		// 温度差係数が1.0でない場合はk_ei_jsに値を代入する。
		// id は自分自身の境界IDとし、自分自身の表面の影響は1.0から温度差係数を減じた値になる。
		if h < 1.0 {
			k_s_r_j_is[connected_room_id] = math.Round((1.0-h)*100) / 100
		} else {
			// 温度差係数が1.0の場合は裏面の影響は何もないため k_ei_js に操作は行わない。
			//pass
		}

	} else {

		//pass
	}

	return &Boundary{
		id:                boundary_id,
		name:              name,
		sub_name:          sub_name,
		connected_room_id: connected_room_id,
		boundary_type:     boundary_type,
		area:              area,
		h_td:              h_td,
		//rear_surface_boundary_id: rear_surface_boundary_id,
		is_floor:                 is_floor,
		is_solar_absorbed_inside: is_solar_absorbed_inside,
		is_sun_striked_outside:   is_sun_striked_outside,
		h_s_c:                    h_s_c,
		h_s_r:                    h_s_r,
		simulation_u_value:       simulation_u_value,
		theta_o_sol:              theta_o_eqv_j_ns,
		q_trs_sol:                q_trs_sol,
		rf:                       rf,
		k_ei_js_j:                k_ei_js_j,
		k_s_r_j_is:               k_s_r_j_is,
	}
}

func _get_p_is_js(n_rm int, bss []*Boundary) *mat.Dense {
	// 室iと境界jの関係を表す係数（境界jから室iへの変換）
	// [[p_0_0 ... ... p_0_j]
	//  [ ...  ... ...  ... ]
	//  [p_i_0 ... ... p_i_j]]

	p_is_js := mat.NewDense(n_rm, len(bss), nil)

	for i, bs := range bss {
		p_is_js.Set(bs.connected_room_id, i, 1)
	}

	return p_is_js
}

func get_q_trs_sol_is_ns(n_rm int, bss []*Boundary) *mat.Dense {
	n := len(bss[0].q_trs_sol)
	q_trs_sol_is_ns := mat.NewDense(n_rm, n, nil)
	for i := 0; i < n_rm; i++ {
		data := make([]float64, n)
		for _, bs := range bss {
			if bs.connected_room_id == i {
				floats.Add(data, bs.q_trs_sol)
			}
		}
		q_trs_sol_is_ns.SetRow(i, data)
	}
	return q_trs_sol_is_ns
}

func (self *Boundaries) get_room_id_by_boundary_id(boundary_id int) int {

	bs := self._get_boundary_by_id(boundary_id)

	return bs.connected_room_id
}

func (self *Boundaries) _get_boundary_by_id(boundary_id int) *Boundary {

	// 指定された boundary_id に一致する Boundary を取得する。
	bss := make([]*Boundary, 0, len(self.bss))
	for _, bs := range self.bss {
		if bs.id == boundary_id {
			bss = append(bss, bs)
		}
	}

	// 取得された Boundary は必ず1つのはずなので、「見つからない場合」「複数該当した場合」にはエラーを出す。
	if len(bss) == 0 {
		panic("指定された boundary_id に一致する boundary が見つかりませんでした。")
	}
	if len(bss) > 1 {
		panic("指定された boundary_id に一致する boundary が複数見つかりました。")
	}

	return bss[0]
}

/*
室内側熱伝達抵抗を取得する。
Args:

	b: 境界の辞書
	boundary_id: 境界のID

Returns:

	室内側熱伝達抵抗, m2K/W
*/
func _read_r_i_nominal(b *BoudaryJson, boundary_id int) float64 {
	// 室内側熱伝達抵抗, m2K/W
	r_i := b.InsideHeatTransferResistance

	if r_i <= 0.0 {
		panic(fmt.Sprintf("境界(ID=%d)の室内側熱伝達抵抗で00.0以下の値が指定されました。", boundary_id))
	}

	return r_i
}

func _read_cs_j_l(layer *LayerJson) float64 {
	return layer.ThermalCapacity
}

func _read_rs_j_l(layer *LayerJson) float64 {
	return layer.ThermalResistance
}

/*
境界jの室外側日射吸収率を取得する。

	Args:
	    b: 境界を表す辞書
	Returns:
	    境界jの室外側日射吸収率, -
*/
func _read_a_s_j(b *BoudaryJson, boundary_id int) float64 {

	a_s := b.OutsideSolarAbsorption
	if a_s < 0.0 {
		panic(fmt.Sprintf("境界(ID=%d)の日射吸収率で0.0未満の値が指定されました。", boundary_id))
	}
	if a_s > 1.0 {
		panic(fmt.Sprintf("境界(ID=%d)の日射吸収率で1.0より大の値が指定されました。", boundary_id))
	}

	return a_s
}

func _read_u_nominal_j(b *BoudaryJson, boundary_id int) float64 {

	u_nominal_j := b.UValue

	if u_nominal_j <= 0.0 {
		panic(fmt.Sprintf("境界(ID=%d)の熱貫流率で0.0以下の値が指定されました。", boundary_id))
	}

	return u_nominal_j
}

func _read_r_s_o_j(b *BoudaryJson, boundary_id int) float64 {

	r_surf := b.OutsideHeatTransferResistance
	if r_surf <= 0.0 {
		panic(fmt.Sprintf("境界(ID=%d)の室外側熱伝達抵抗で0.0以下の値が指定されました。", boundary_id))
	}

	return r_surf
}

func _read_eps_r_o_j(b *BoudaryJson, boundary_id int) float64 {
	eps_r := b.OutsideEmissivity

	if eps_r > 1.0 {
		panic(fmt.Sprintf("境界(ID=%d)の室外側長波長放射率で1.0を超える値が指定されました。", boundary_id))
	}
	if eps_r < 0.0 {
		panic(fmt.Sprintf("境界(ID=%d)の室外側長波長放射率で0.0を下回る値が指定されました。", boundary_id))
	}

	return eps_r
}
