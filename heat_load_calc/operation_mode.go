package main

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

type ACMethod string

const (
	PMV             ACMethod = "pmv"             // PMV 制御
	SIMPLE          ACMethod = "simple"          // OT制御と同じ　互換性を保つために残しているがOTを使用することを推奨する
	OT              ACMethod = "ot"              // OT制御
	AIR_TEMPERATURE ACMethod = "air_temperature" // 空気温度制御
)

// 運転状態
type OperationMode int

// 運転状態
const (
	COOLING    OperationMode = iota + 1 // COOLING ： 冷房
	HEATING                             // HEATING : 暖房
	STOP_OPEN                           // STOP_OPEN : 暖房・冷房停止で窓「開」
	STOP_CLOSE                          // STOP_CLOSE : 暖房・冷房停止で窓「閉」
)

type Operation struct {
	_ac_method          ACMethod
	_ac_config          []ACConfigJson
	_lower_target_is_ns *ScheduleData
	_upper_target_is_ns *ScheduleData
	_ac_demand_is_ns    *ScheduleData
	_n_rm               int
}

func NewOperation(
	ac_method ACMethod,
	ac_config []ACConfigJson,
	lower_target_is_ns *ScheduleData,
	upper_target_is_ns *ScheduleData,
	ac_demand_is_ns *ScheduleData,
	n_rm int,
) *Operation {
	return &Operation{
		_ac_method:          ac_method,
		_ac_config:          ac_config,
		_lower_target_is_ns: lower_target_is_ns,
		_upper_target_is_ns: upper_target_is_ns,
		_ac_demand_is_ns:    ac_demand_is_ns,
		_n_rm:               n_rm,
	}
}

func make_operation(
	d *CommonJson,
	ac_setting_is_ns *ScheduleData,
	ac_demand_is_ns *ScheduleData,
	n_rm int,
) *Operation {
	ac_method := ACMethod(d.AcMethod)

	var ac_config []ACConfigJson
	if d.AcConfig != nil {
		ac_config = d.AcConfig
	} else {
		switch ac_method {
		case AIR_TEMPERATURE, SIMPLE, OT:
			ac_config = []ACConfigJson{
				ACConfigJson{
					Mode:  1,
					Lower: 20.0,
					Upper: 27.0,
				},
				ACConfigJson{
					Mode:  2,
					Lower: 20.0,
					Upper: 27.0,
				},
			}
		case PMV:
			ac_config = []ACConfigJson{
				ACConfigJson{
					Mode:  1,
					Lower: -0.5,
					Upper: 0.5,
				},
				ACConfigJson{
					Mode:  2,
					Lower: -0.5,
					Upper: 0.5,
				},
			}
		default:
			panic("invalid ac_method")
		}
	}

	r, c := ac_setting_is_ns.BatchSize, ac_setting_is_ns.Len()

	lower_target_is_ns_data := make([]float64, len(ac_setting_is_ns.Data))
	upper_target_is_ns_data := make([]float64, len(ac_setting_is_ns.Data))
	for i := 0; i < len(ac_setting_is_ns.Data); i++ {
		lower_target_is_ns_data[i] = math.NaN()
		upper_target_is_ns_data[i] = math.NaN()
	}

	to_int := func(v interface{}) int {
		switch v.(type) {
		case float64:
			return int(v.(float64))
		case int:
			return v.(int)
		default:
			panic("invalid type")
		}
	}

	for _, conf := range ac_config {
		mode := to_int(conf.Mode)
		lower := conf.Lower
		upper := conf.Upper

		off := 0
		for j := 0; j < c; j++ {
			for i := 0; i < r; i++ {
				batch := ac_setting_is_ns.Get(j)
				if int(batch.AtVec(i)) == mode {
					lower_target_is_ns_data[off] = lower
					upper_target_is_ns_data[off] = upper
				}
				off++
			}
		}
	}

	lower_target_is_ns := &ScheduleData{
		Data:      lower_target_is_ns_data,
		BatchSize: r,
	}

	upper_target_is_ns := &ScheduleData{
		Data:      upper_target_is_ns_data,
		BatchSize: r,
	}

	return NewOperation(
		ac_method,
		ac_config,
		lower_target_is_ns,
		upper_target_is_ns,
		ac_demand_is_ns,
		n_rm,
	)
}

func (o *Operation) ac_method() ACMethod {
	return o._ac_method
}

func (o *Operation) ac_config() []ACConfigJson {
	return o._ac_config
}

/*
		運転モードを決定する。

	    Args:
	        ac_demand_is_n: ステップnにおける室iの空調需要の有無, 0.0~1.0, [i, 1]
	        operation_mode_is_n_mns: ステップn-1における室iの運転状態, [i, 1]
	        pmv_heavy_is_n: ステップnにおける室iの厚着時のPMV, [i, 1]
	        pmv_middle_is_n: ステップnにおける室iの中間着時のPMV, [i, 1]
	        pmv_light_is_n: ステップnにおける室iの薄着時のPMV, [i, 1]

	    Returns:
	        ステップnの室iにおける運転状態, [i, 1]
*/
func (self *Operation) get_operation_mode_is_n(
	n int,
	nn int,
	is_radiative_heating_is []bool,
	is_radiative_cooling_is []bool,
	met_is []float64,
	theta_r_ot_ntr_non_nv_is_n_pls *mat.VecDense,
	theta_r_ot_ntr_nv_is_n_pls *mat.VecDense,
	theta_r_ntr_non_nv_is_n_pls []float64,
	theta_r_ntr_nv_is_n_pls []float64,
	theta_mrt_hum_ntr_non_nv_is_n_pls []float64,
	theta_mrt_hum_ntr_nv_is_n_pls []float64,
	x_r_ntr_non_nv_is_n_pls []float64,
	x_r_ntr_nv_is_n_pls []float64,
) ([]OperationMode, bool) {

	upper_target_is_n := self._upper_target_is_ns.Get(nn)
	lower_target_is_n := self._lower_target_is_ns.Get(nn)
	ac_demand_is_n := self._ac_demand_is_ns.Get(nn)

	var x_cooling_is_n_pls, x_window_open_is_n_pls, x_heating_is_n_pls []float64

	switch self._ac_method {
	case AIR_TEMPERATURE, SIMPLE, OT:
		x_cooling_is_n_pls, x_window_open_is_n_pls, x_heating_is_n_pls = _get_operation_mode_simple_is_n(
			theta_r_ot_ntr_non_nv_is_n_pls.RawVector().Data,
			theta_r_ot_ntr_nv_is_n_pls.RawVector().Data,
		)

	case PMV:
		x_cooling_is_n_pls, x_window_open_is_n_pls, x_heating_is_n_pls = _get_operation_mode_pmv_is_n(
			is_radiative_cooling_is,
			is_radiative_heating_is,
			"constant",
			met_is,
			theta_r_ntr_non_nv_is_n_pls,
			theta_r_ntr_nv_is_n_pls,
			theta_mrt_hum_ntr_non_nv_is_n_pls,
			theta_mrt_hum_ntr_nv_is_n_pls,
			x_r_ntr_non_nv_is_n_pls,
			x_r_ntr_nv_is_n_pls,
		)
	default:
		panic("invalid ac method")
	}

	all_stop := true
	v := make([]OperationMode, self._n_rm)
	for i := 0; i < self._n_rm; i++ {
		is_op := ac_demand_is_n.AtVec(i) > 0.0
		if is_op {
			if x_cooling_is_n_pls[i] > upper_target_is_n.AtVec(i) && x_window_open_is_n_pls[i] > upper_target_is_n.AtVec(i) {
				v[i] = COOLING
				all_stop = false
			} else if x_cooling_is_n_pls[i] > upper_target_is_n.AtVec(i) && x_window_open_is_n_pls[i] <= upper_target_is_n.AtVec(i) {
				v[i] = STOP_OPEN
			} else if x_heating_is_n_pls[i] < lower_target_is_n.AtVec(i) {
				v[i] = HEATING
				all_stop = false
			} else {
				v[i] = STOP_CLOSE
			}
		} else {
			v[i] = STOP_CLOSE
		}
	}

	return v, all_stop
}

func (self *Operation) get_theta_target_is_n(
	operation_mode_is_n []OperationMode,
	theta_r_is_n []float64,
	theta_mrt_hum_is_n []float64,
	x_r_ntr_is_n_pls []float64,
	n int,
	nn int,
	is_radiative_heating_is []bool,
	is_radiative_cooling_is []bool,
	met_is []float64,
) (*mat.VecDense, *mat.VecDense, []float64, []float64) {
	lower_target_is_n := self._lower_target_is_ns.Get(nn)
	upper_target_is_n := self._upper_target_is_ns.Get(nn)

	// ステップnの室iにおけるClo値, [i, 1]
	clo_is_n := get_clo_is_ns(operation_mode_is_n)

	// ステップnにおける室iの在室者周りの風速, m/s, [i, 1]
	v_hum_is_n := get_v_hum_is_n(
		operation_mode_is_n,
		is_radiative_heating_is,
		is_radiative_cooling_is,
	)

	// (1) ステップ n における室 i の在室者周りの対流熱伝達率, W/m2K, [i, 1]
	// (2) ステップ n における室 i の在室者周りの放射熱伝達率, W/m2K, [i, 1]
	// (3) ステップ n における室 i の在室者周りの総合熱伝達率, W/m2K, [i, 1]
	h_hum_c_is_n, h_hum_r_is_n, h_hum_is_n := get_h_hum(
		theta_mrt_hum_is_n,
		theta_r_is_n,
		clo_is_n,
		v_hum_is_n,
		"constant",
		met_is,
	)

	switch self.ac_method() {
	case AIR_TEMPERATURE, SIMPLE, OT:
		return lower_target_is_n, upper_target_is_n, h_hum_c_is_n, h_hum_r_is_n
	case PMV:

		// ステップnにおける室iの水蒸気圧, Pa, [i, 1]
		p_v_r_is_n := get_p_v_r_is_n(x_r_ntr_is_n_pls)

		theta_lower_target_is_n, theta_upper_target_is_n := _get_theta_target(
			operation_mode_is_n,
			p_v_r_is_n,
			met_is,
			lower_target_is_n,
			upper_target_is_n,
			h_hum_is_n,
			clo_is_n,
		)

		return theta_lower_target_is_n, theta_upper_target_is_n, h_hum_c_is_n, h_hum_r_is_n

	default:
		panic("invalid ac method")
	}
}

/*
Returns:

	ステップ n における室 i の人体表面の対流熱伝達率が総合熱伝達率に占める割合, -, [i, 1]
	ステップ n における室 i の人体表面の放射熱伝達率が総合熱伝達率に占める割合, -, [i, 1]
*/
func (o *Operation) get_k_is() (*mat.VecDense, *mat.VecDense) {
	k_c_is := mat.NewVecDense(o._n_rm, nil)
	k_r_is := mat.NewVecDense(o._n_rm, nil)
	switch o.ac_method() {
	case AIR_TEMPERATURE:
		for i := 0; i < o._n_rm; i++ {
			k_c_is.SetVec(i, 1.0)
			k_r_is.SetVec(i, 0.0)
		}
	case OT, PMV, SIMPLE:
		for i := 0; i < o._n_rm; i++ {
			k_c_is.SetVec(i, 0.5)
			k_r_is.SetVec(i, 0.5)
		}
	default:
		panic("invalid ac_method")
	}

	return k_c_is, k_r_is
}

func _get_theta_target_simple(
	p_v_r_is_n mat.Vector,
	operation_mode_is_n []OperationMode,
	theta_r_is_n mat.Vector,
	theta_mrt_hum_is_n mat.Vector,
	theta_lower_target_is_ns mat.Matrix,
	theta_upper_target_is_ns mat.Matrix,
	n int) (
	mat.Vector,
	mat.Vector,
	mat.Vector,
	mat.Vector,
	mat.Vector,
	mat.Vector,
) {

	// ステップnの室iにおけるClo値, [i, 1]
	l := theta_r_is_n.Len()
	clo_is_n := mat.NewVecDense(l, nil)
	for i := 0; i < l; i++ {
		clo_is_n.SetVec(i, clo_light)
	}

	// ステップnにおける室iの在室者周りの風速, m/s, [i, 1]
	v_hum_is_n := mat.NewVecDense(l, nil)
	h_hum_c_is_n := mat.NewVecDense(l, nil)
	h_hum_r_is_n := mat.NewVecDense(l, nil)
	for i := 0; i < l; i++ {
		h_hum_c_is_n.SetVec(i, 1.0)
		h_hum_r_is_n.SetVec(i, 1.0)
	}

	if n < 0 {
		_, c := theta_lower_target_is_ns.Dims()
		n += c
	}
	theta_lower_target_is_n := theta_lower_target_is_ns.(mat.ColViewer).ColView(n)
	theta_upper_target_is_n := theta_upper_target_is_ns.(mat.ColViewer).ColView(n)

	return theta_lower_target_is_n, theta_upper_target_is_n, h_hum_c_is_n, h_hum_r_is_n, v_hum_is_n, clo_is_n
}

func _get_theta_target(
	operation_mode_is_n []OperationMode,
	p_v_r_is_n []float64,
	met_is []float64,
	lower_target_is_n mat.Vector,
	upper_target_is_n mat.Vector,
	h_hum_is_n []float64,
	clo_is_n []float64,
) (
	*mat.VecDense,
	*mat.VecDense,
) {
	theta_lower_target_is_n := mat.NewVecDense(len(operation_mode_is_n), nil)
	theta_upper_target_is_n := mat.NewVecDense(len(operation_mode_is_n), nil)
	for i, om := range operation_mode_is_n {
		if om == HEATING {
			v := get_theta_ot_target(
				clo_is_n[i],
				p_v_r_is_n[i],
				h_hum_is_n[i],
				met_is[i],
				lower_target_is_n.AtVec(i),
			)
			theta_lower_target_is_n.SetVec(i, v)
		} else if om == COOLING {
			v := get_theta_ot_target(
				clo_is_n[i],
				p_v_r_is_n[i],
				h_hum_is_n[i],
				met_is[i],
				lower_target_is_n.AtVec(i),
			)
			theta_upper_target_is_n.SetVec(i, v)
		} else {
			//PASS
		}
	}

	return theta_lower_target_is_n, theta_upper_target_is_n
}

/*
		在室者周りの風速を求める。

	    Args:
	        operation_mode_is:
	        is_radiative_cooling_is:
	        is_radiative_heating_is:

	    Returns:
	        ステップnにおける室iの在室者周りの風速, m/s, [i, 1]
*/
func get_v_hum_is_n(
	operation_mode_is []OperationMode,
	is_radiative_cooling_is []bool,
	is_radiative_heating_is []bool,
) []float64 {
	// 在室者周りの風速はデフォルトで 0.0 m/s とおく
	v_hum_is_n := make([]float64, len(operation_mode_is))

	for i := 0; i < len(v_hum_is_n); i++ {
		if operation_mode_is[i] == HEATING {
			if is_radiative_heating_is[i] {
				// 対流暖房時の風速を 0.2 m/s とする
				v_hum_is_n[i] = 0.2
			} else {
				// 放射暖房時の風速を 0.0 m/s とする
				v_hum_is_n[i] = 0.0
			}
		} else if operation_mode_is[i] == COOLING {
			if is_radiative_cooling_is[i] {
				// 対流冷房時の風速を 0.2 m/s とする
				v_hum_is_n[i] = 0.2
			} else {
				// 放射冷房時の風速を 0.0 m/s とする
				v_hum_is_n[i] = 0.0
			}
		}

		// 暖冷房をせずに窓を開けている時の風速を 0.1 m/s とする
		// 対流暖房・冷房時と窓を開けている時は同時には起こらないことを期待しているが
		// もし同時にTrueの場合は窓を開けている時の風速が優先される（上書きわれる）
		if operation_mode_is[i] == STOP_OPEN {
			v_hum_is_n[i] = 0.1
		}

		// 上記に当てはまらない場合の風速は 0.0 m/s のままである。
	}

	return v_hum_is_n
}

/*
運転モードに応じた在室者のClo値を決定する。

Args

	operation_mode_is_n ステップnにおける室iの運転状態, [i, 1]

Returns

	ステップnにおける室iの在室者のClo値, [i, 1]
*/
func get_clo_is_ns(operation_mode_is_n []OperationMode) []float64 {

	// ステップnにおける室iの在室者のClo値, [i, 1]
	clo_is_n := make([]float64, len(operation_mode_is_n))

	// 運転方法に応じてclo値の設定を決定する。
	for i, v := range operation_mode_is_n {
		switch v {
		case HEATING:
			// 暖房運転時の場合は厚着とする。
			clo_is_n[i] = clo_heavy
		case COOLING:
			// 冷房運転時の場合は薄着とする。
			clo_is_n[i] = clo_light
		case STOP_OPEN, STOP_CLOSE:
			// 暖冷房停止時は、窓の開閉によらず中間着とする。
			clo_is_n[i] = clo_middle
		}
	}

	return clo_is_n
}

func _get_operation_mode_simple_is_n(
	theta_r_ot_ntr_non_nv_is_n_pls []float64,
	theta_r_ot_ntr_nv_is_n_pls []float64,
) ([]float64, []float64, []float64) {

	x_cooling_is_n_pls := theta_r_ot_ntr_nv_is_n_pls
	x_window_open_is_n_pls := theta_r_ot_ntr_non_nv_is_n_pls
	x_heating_is_n_pls := theta_r_ot_ntr_nv_is_n_pls

	return x_cooling_is_n_pls, x_window_open_is_n_pls, x_heating_is_n_pls
}

func _get_operation_mode_pmv_is_n(
	is_radiative_cooling_is []bool,
	is_radiative_heating_is []bool,
	method string,
	met_is []float64,
	theta_r_ntr_non_nv_is_n_pls []float64,
	theta_r_ntr_nv_is_n_pls []float64,
	theta_mrt_hum_ntr_non_nv_is_n_pls []float64,
	theta_mrt_hum_ntr_nv_is_n_pls []float64,
	x_r_ntr_non_nv_is_n_pls []float64,
	x_r_ntr_nv_is_n_pls []float64,
) ([]float64, []float64, []float64) {
	// ステップnにおける室iの水蒸気圧, Pa, [i, 1]
	p_v_r_ntr_non_nv_is_n_pls := get_p_v_r_is_n(x_r_ntr_non_nv_is_n_pls)
	p_v_r_ntr_nv_is_n_pls := get_p_v_r_is_n(x_r_ntr_nv_is_n_pls)

	// 薄着時のClo値
	clo_light_slice := make([]float64, len(is_radiative_cooling_is))
	for i := 0; i < len(is_radiative_cooling_is); i++ {
		clo_light_slice[i] = clo_light
	}

	// 厚着時のClo値
	clo_heavy_slice := make([]float64, len(is_radiative_cooling_is))
	for i := 0; i < len(is_radiative_cooling_is); i++ {
		clo_heavy_slice[i] = clo_heavy
	}

	////// 冷房判定用（窓開け時）のPMV計算

	// 窓を開けている時の風速を 0.1 m/s とする
	v_hum_window_open_is_n := make([]float64, len(is_radiative_cooling_is))
	for i := 0; i < len(is_radiative_cooling_is); i++ {
		v_hum_window_open_is_n[i] = 0.1
	}

	// 冷房判定用（窓開け時）のPMV
	pmv_window_open_is_n := get_pmv_is_n(
		p_v_r_ntr_nv_is_n_pls,
		theta_r_ntr_nv_is_n_pls,
		theta_mrt_hum_ntr_nv_is_n_pls,
		clo_light_slice,
		v_hum_window_open_is_n,
		met_is,
		method,
	)

	// 冷房判定用（窓閉め時）のPMV計算

	// 冷房時の風速を対流冷房時0.2m/s・放射冷房時0.0m/sに設定する。
	v_hum_cooling_is_n := make([]float64, len(is_radiative_cooling_is))
	for i := 0; i < len(is_radiative_cooling_is); i++ {
		if is_radiative_cooling_is[i] {
			v_hum_cooling_is_n[i] = 0.0
		} else {
			v_hum_cooling_is_n[i] = 0.2
		}
	}

	// 冷房判定用（窓閉め時）のPMV
	pmv_cooling_is_n := get_pmv_is_n(
		p_v_r_ntr_non_nv_is_n_pls,
		theta_r_ntr_non_nv_is_n_pls,
		theta_mrt_hum_ntr_non_nv_is_n_pls,
		clo_light_slice,
		v_hum_cooling_is_n,
		met_is,
		method,
	)

	// 暖房判定用のPMV計算

	// 暖房時の風速を対流暖房時0.2m/s・放射暖房時0.0m/sに設定する。
	v_hum_heating_is_n := make([]float64, len(is_radiative_heating_is))
	for i := 0; i < len(is_radiative_heating_is); i++ {
		if is_radiative_heating_is[i] {
			v_hum_heating_is_n[i] = 0.0
		} else {
			v_hum_heating_is_n[i] = 0.2
		}
	}

	// 暖房判定用のPMV
	pmv_heating_is_n := get_pmv_is_n(
		p_v_r_ntr_non_nv_is_n_pls,
		theta_r_ntr_non_nv_is_n_pls,
		theta_mrt_hum_ntr_non_nv_is_n_pls,
		clo_heavy_slice,
		v_hum_heating_is_n,
		met_is,
		method,
	)

	x_cooling_is_n_pls := pmv_cooling_is_n
	x_window_open_is_n_pls := pmv_window_open_is_n
	x_heating_is_n_pls := pmv_heating_is_n

	return x_cooling_is_n_pls, x_window_open_is_n_pls, x_heating_is_n_pls
}

/*
運転モードを決定する。

	Args:
	    ac_demand_i_n: ステップnにおける室iの空調需要の有無, 0.0～1.0
	    operation_mode_i_n_mns: ステップn-1における室iの運転状態
	    pmv_heavy_i_n: ステップnにおける室iの厚着時のPMV
	    pmv_middle_i_n: ステップnにおける室iの中間着時のPMV
	    pmv_light_i_n: ステップnにおける室iの薄着時のPMV

	Returns:
	    ステップnにおける室iの運転状態
*/
func get_operation_mode_i_n(
	ac_demand_i_n float64,
	operation_mode_i_n_mns OperationMode,
	pmv_heavy_i_n float64,
	pmv_middle_i_n float64,
	pmv_light_i_n float64,
) OperationMode {

	if ac_demand_i_n > 0.0 { // 空調需要がある場合
		switch operation_mode_i_n_mns {
		case HEATING:
			if pmv_heavy_i_n <= 0.7 {
				return HEATING
			} else if pmv_light_i_n >= 0.7 {
				return COOLING
			} else {
				return STOP_CLOSE
			}
		case COOLING:
			if pmv_light_i_n >= -0.7 {
				return COOLING
			} else if pmv_heavy_i_n <= -0.7 {
				return HEATING
			} else {
				return STOP_CLOSE
			}
		case STOP_OPEN:
			if pmv_light_i_n >= 0.7 {
				return COOLING
			} else if pmv_heavy_i_n <= -0.7 {
				return HEATING
			} else if pmv_middle_i_n <= 0.0 {
				return STOP_CLOSE
			} else {
				return STOP_OPEN
			}
		case STOP_CLOSE:
			if pmv_light_i_n >= 0.7 {
				return COOLING
			} else if pmv_heavy_i_n <= -0.7 {
				return HEATING
			} else if pmv_middle_i_n >= 0.0 {
				return STOP_OPEN
			} else {
				return STOP_CLOSE
			}
		default:
			panic("Invalid operation_mode_i_n_mns")
		}
	} else { // 空調需要がない場合（窓閉鎖、空調停止）
		return STOP_CLOSE
	}
}
