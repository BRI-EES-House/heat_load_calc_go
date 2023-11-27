package heat_load_calc

import (
	"gonum.org/v1/gonum/mat"
)

type VentilationType string

const (
	VentilationTypeTYPE1 VentilationType = "type1"
	VentilationTypeTYPE2 VentilationType = "type2"
	VentilationTypeTYPE3 VentilationType = "type3"
	NATURAL_LOOP         VentilationType = "natural_loop"
)

// 機械換気の要素
type MechanicalVentilation struct {
	id        int             // ID
	root_type VentilationType // 換気経路のタイプ
	volume    float64         // 換気量, m3/h
	root      []int           // 換気のルート
}

// 機械換気の全体
type MechanicalVentilations struct {
	_mechanical_ventilations []*MechanicalVentilation // 機械換気の要素
	_n_rm                    int                      // 室数
	v_vent_int_is_is         *mat.Dense               // 室間換気量, m3/h
	v_vent_mec_general_is    []float64                // 室ごとの機械換気量, m3/h
}

func NewMechanicalVentilations(vs []MechanicalVentilationJson, n_rm int) *MechanicalVentilations {
	mvs := &MechanicalVentilations{_n_rm: n_rm}
	for i := range vs {
		input := &vs[i]
		mv := &MechanicalVentilation{
			id:        int(input.Id),                   // ID
			root_type: VentilationType(input.RootType), // 換気経路のタイプ
			volume:    input.Volume,                    // 換気量, m3/h
			root:      input.Root,                      // 換気のルート
		}
		mvs._mechanical_ventilations = append(mvs._mechanical_ventilations, mv)
	}

	mvs.v_vent_mec_general_is = mvs.get_v_vent_mec_general_is()
	mvs.v_vent_int_is_is = mvs.get_v_vent_int_is_is()

	return mvs
}

/*
Calculate the mechanical ventilation (general) of room i.

Returns:

	mechanical ventilation (general) of room i, m3/s

Note:

	eq. 1
*/
func (mvs *MechanicalVentilations) get_v_vent_mec_general_is() []float64 {
	v1 := make([]float64, mvs._n_rm)

	for _, v := range mvs._mechanical_ventilations {
		switch v.root_type {

		// In case that the mechanical ventialtion tyepe is type1, type2 or type3,
		// the outdoor air flows into the room with top ID.
		// Therefore, add the ventilation amount to the matrix of the room ID.
		case VentilationTypeTYPE1, VentilationTypeTYPE2, VentilationTypeTYPE3:
			// the list of the mechanical ventilation loot
			r := v.root

			// Add the ventilation amount to the room which ID equals to the top ID of mechanical ventilation list.
			// Convert the unit from m3/h to m3/s.
			v1[r[0]] += v.volume / 3600
		case NATURAL_LOOP:
			// PASS
		default:
			panic(v.root_type)
		}
	}

	return v1
}

/*
Calculate the tamount of the air moved from the room i* to the room i.

Returns:

	the amount of the air moved from the room i* to tohe room i, [i, i], m3/s
*/
func (mvs *MechanicalVentilations) get_v_vent_int_is_is() *mat.Dense {
	v2 := mat.NewDense(mvs._n_rm, mvs._n_rm, nil)

	for _, v := range mvs._mechanical_ventilations {
		switch v.root_type {
		case VentilationTypeTYPE1, VentilationTypeTYPE2, VentilationTypeTYPE3:
			r := v.root
			for i := 1; i < len(r); i++ {
				v2.Set(r[i], r[i-1], v2.At(r[i], r[i-1])+v.volume/3600)
				v2.Set(r[i], r[i], v2.At(r[i], r[i])-v.volume/3600)
			}
		case NATURAL_LOOP:
			r := v.root
			for i := 1; i < len(r); i++ {
				var i_upstream int
				if i == 0 {
					i_upstream = len(r) - 1
				} else {
					i_upstream = i - 1
				}
				v2.Set(r[i], r[i_upstream], v2.At(r[i], r[i_upstream])+v.volume/3600)
				v2.Set(r[i], r[i], v2.At(r[i], r[i])-v.volume/3600)
			}
		default:
			panic(v.root_type)
		}
	}

	return v2
}
