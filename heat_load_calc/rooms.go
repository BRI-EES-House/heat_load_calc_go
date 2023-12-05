package heat_load_calc

import (
	"gonum.org/v1/gonum/mat"
)

type Room struct {
	id             int     // id
	name           string  // 名称
	sub_name       string  // 副名称
	v              float64 // 気積, m3
	a_f            float64 // 床面積, m2
	c_sh_frt       float64 // 備品等の熱容量, J/K
	g_sh_frt       float64 // 空気と備品等間の熱コンダクタンス, W/K
	c_lh_frt       float64 // 備品等の湿気容量, kg/(kg/kgDA)
	g_lh_frt       float64 // 空気と備品等間の湿気コンダクタンス, kg/(kg/kgDA)
	v_vent_ntr_set float64 // 自然風利用時の換気量, m3/s
}

type Rooms struct {
	n_rm              int           // 室の数
	rms               []Room        // 室, [I]
	id_rm_is          []int         // 空間のID, [I]
	name_rm_is        []string      // 室 i の名前, [I]
	sub_name_rm_is    []string      // 室 i の名前2, [I]
	a_f_rm_is         *mat.VecDense // 室 i の面積, m2, [I]
	v_rm_is           *mat.VecDense // 室 i の容積, m3, [I]
	c_sh_frt_is       *mat.VecDense // 室 i の備品等の熱容量, J/K, [I]
	g_sh_frt_is       *mat.VecDense // 室 i の空気と備品等間の熱コンダクタンス, W/K, [I]
	c_lh_frt_is       *mat.VecDense // 室 i の備品等の湿気容量, kg/(kg/kgDA), [I]
	g_lh_frt_is       *mat.VecDense // 室 i の空気と備品等間の湿気コンダクタンス, kg/(s (kg/kgDA)), [I]
	v_vent_ntr_set_is []float64     // 室 i の自然風利用時の換気量, m3/s, [I]
	met_is            []float64     // 室 i の在室者のMet値, [I]
}

func NewRooms(ds []RoomJson) (*Rooms, error) {
	n_rm := len(ds)
	rms := make([]Room, n_rm)
	id_rm_is := make([]int, n_rm)
	name_rm_is := make([]string, n_rm)
	sub_name_rm_is := make([]string, n_rm)
	v_rm_is := make([]float64, n_rm)
	floor_area_rm_is := make([]float64, n_rm)
	c_sh_frt_rm_is := make([]float64, n_rm)
	g_sh_frt_rm_is := make([]float64, n_rm)
	c_lh_frt_rm_is := make([]float64, n_rm)
	g_lh_frt_rm_is := make([]float64, n_rm)
	v_vent_ntr_set_rm_is := make([]float64, n_rm)
	met_is := make([]float64, n_rm)

	for i, d := range ds {
		rm, err := _get_rm(d)
		if err != nil {
			return nil, err
		}
		rms[i] = rm
		id_rm_is[i] = rm.id
		name_rm_is[i] = rm.name
		sub_name_rm_is[i] = rm.sub_name
		v_rm_is[i] = rm.v
		floor_area_rm_is[i] = rm.a_f
		c_sh_frt_rm_is[i] = rm.c_sh_frt
		g_sh_frt_rm_is[i] = rm.g_sh_frt
		c_lh_frt_rm_is[i] = rm.c_lh_frt
		g_lh_frt_rm_is[i] = rm.g_lh_frt
		v_vent_ntr_set_rm_is[i] = rm.v_vent_ntr_set
		met_is[i] = 1.0
	}

	return &Rooms{
		n_rm:              n_rm,
		rms:               rms,
		id_rm_is:          id_rm_is,
		name_rm_is:        name_rm_is,
		sub_name_rm_is:    sub_name_rm_is,
		a_f_rm_is:         mat.NewVecDense(n_rm, floor_area_rm_is),
		v_rm_is:           mat.NewVecDense(n_rm, v_rm_is),
		c_sh_frt_is:       mat.NewVecDense(n_rm, c_sh_frt_rm_is),
		g_sh_frt_is:       mat.NewVecDense(n_rm, g_sh_frt_rm_is),
		c_lh_frt_is:       mat.NewVecDense(n_rm, c_lh_frt_rm_is),
		g_lh_frt_is:       mat.NewVecDense(n_rm, g_lh_frt_rm_is),
		v_vent_ntr_set_is: v_vent_ntr_set_rm_is,
		met_is:            met_is,
	}, nil
}

func _get_rm(d RoomJson) (Room, error) {
	v_rm_i := d.Volume

	c_lh_frt, c_sh_frt, g_lh_frt, g_sh_frt, err := get_furniture_specs(
		&d.Furniture,
		v_rm_i,
	)
	if err != nil {
		return Room{}, err
	}

	// v_vent_ntr_set については m3/h から m3/s の単位変換を行う。
	return Room{
		id:             d.Id,
		name:           d.Name,
		sub_name:       d.SubName,
		a_f:            d.FloorArea,
		v:              v_rm_i,
		c_sh_frt:       c_sh_frt,
		g_sh_frt:       g_sh_frt,
		c_lh_frt:       c_lh_frt,
		g_lh_frt:       g_lh_frt,
		v_vent_ntr_set: d.Ventilation.Natural / 3600.0,
	}, nil
}
