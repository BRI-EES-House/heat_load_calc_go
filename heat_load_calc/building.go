package main

import (
	"math"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/stat"
)

// 建物の階数（共同住宅の場合は住戸の階数）
type Story int

// 建物の階数（共同住宅の場合は住戸の階数）
const (
	StoryOne Story = 1 // 1階
	StoryTwo Story = 2 // 2階（2階以上の階数の場合も2階とする。）
)

func (s Story) String() string {
	return [...]string{"", "ONE", "TWO"}[s]
}

func StoryFromString(s string) Story {
	return map[string]Story{
		"ONE": StoryOne,
		"TWO": StoryTwo,
	}[s]
}

//---------------------------------------------------------------------------------------------------//

// 室内圧力
type InsidePressure int

// 室内圧力
const (
	InsidePressurePositive InsidePressure = iota // 正圧
	InsidePressureNegative                       // 負圧
	InsidePressureBalanced                       // ゼロバランス
)

func (ip InsidePressure) String() string {
	return [...]string{"positive", "negative", "balanced"}[ip]
}

func InsidePressureFromString(s string) InsidePressure {
	return map[string]InsidePressure{
		"positive": InsidePressurePositive,
		"negative": InsidePressureNegative,
		"balanced": InsidePressureBalanced,
	}[s]
}

//---------------------------------------------------------------------------------------------------//

// 構造を表す列挙型
type Structure int

// 構造を表す列挙型
const (
	StructureRC     Structure = iota // RC
	StructureSRC                     // SRC
	StructureWooden                  //木造
	StructureSteel                   //鉄骨
)

func (s Structure) String() string {
	return [...]string{"rc", "src", "wooden", "steel"}[s]
}

func StructureFromString(s string) Structure {
	return map[string]Structure{
		"rc":     StructureRC,
		"src":    StructureSRC,
		"wooden": StructureWooden,
		"steel":  StructureSteel,
	}[s]
}

//---------------------------------------------------------------------------------------------------//

type Building struct {
	infiltration_method string
	story               Story
	c_value             float64
	inside_pressure     InsidePressure
}

func NewBuilding(
	infiltration_method string,
	story Story,
	c_value float64,
	inside_pressure InsidePressure,
) *Building {
	return &Building{
		infiltration_method: infiltration_method,
		story:               story,
		c_value:             c_value,
		inside_pressure:     inside_pressure,
	}
}

func CreateBuilding(d *BuildingJson) *Building {
	ifl := d.Infiltration
	ifl_method := ifl.Method

	var story Story
	var c_value float64
	var inside_pressure InsidePressure

	if ifl_method == "balance_residential" {
		// 建物の階数
		story = Story(int(ifl.Story))

		// C値
		switch ifl.CValueEstimate {
		case "specify":
			c_value = ifl.CValue
		case "calculate":
			ua_value := ifl.UAValue
			structure := StructureFromString(ifl.Structure)
			c_value = _estimate_c_value(ua_value, structure)
		default:
			panic(ifl.CValueEstimate)
		}

		// 換気の種類
		inside_pressure = InsidePressureFromString(ifl.InsidePressure)
	} else {
		panic(ifl_method)
	}

	return NewBuilding(
		ifl_method,
		story,
		c_value,
		inside_pressure,
	)
}

/*
Calculate the leakage air volume
This calculation is approx. expression based on the elaborate results obtained by solving for pressure balance
Args:

	theta_r_is_n: room temperature of room i in step n, degree C, [i,1]
	theta_o_n: outdoor temperature at step n, degree C
	v_room_is: room volume of room i, m3, [i,1]

Returns:

	leakage air volume of rooms at step n, m3/s, [i,1]
*/
func (self *Building) get_v_leak_is_n(
	theta_r_is_n []float64,
	theta_o_n float64,
	v_rm_is []float64,
) []float64 {

	// average air temperature at step n which is weghted by room volumes, degree C
	theta_average_r_n := _get_theta_average_r_n(theta_r_is_n, v_rm_is)

	// temperature difference between room and outdoor at step n, K
	delta_theta_n := _get_delta_theta_n(theta_average_r_n, theta_o_n)

	// ventilation rate of air leakage at step n, 1/h
	n_leak_n := _get_n_leak_n(
		self.c_value,
		self.story,
		self.inside_pressure,
		delta_theta_n,
	)

	// leakage air volume of rooms at step n, m3/s, [i, 1]
	v_leak_is_n := _get_v_leak_is_n(n_leak_n, v_rm_is)

	return v_leak_is_n
}

/*
Args

	ua_value: UA値, W/m2 K
	struct: 構造

Returns:

	C値, cm2/m2
*/
func _estimate_c_value(uaValue float64, structure Structure) float64 {
	a := map[Structure]float64{
		StructureRC:     4.16, // RC造
		StructureSRC:    4.16, // SRC造
		StructureWooden: 8.28, // 木造
		StructureSteel:  8.28, // 鉄骨造
	}[structure]

	return a * uaValue
}

/*
		calculate leakage air volume of rooms at step n
	    Args:
	        n_leak_n: ventilation rate of air leakage at step n, 1/h
	        v_rm_is: room volume of rooms, m3, [i, 1]
	    Returns:
	        air leakage volume of rooms at step n, m3/s, [i, 1]
	    Note:
	        eq.2
*/
func _get_v_leak_is_n(n_leak_n float64, v_rm_is []float64) []float64 {
	v_leak_is_n := make([]float64, len(v_rm_is))

	floats.ScaleTo(v_leak_is_n, n_leak_n/3600, v_rm_is)

	return v_leak_is_n
}

/*
		Calculate the leakage air volume
	    This calculation is approx. expression based on the elaborate results obtained by solving for pressure balance
	    Args:
	        c_value: equivalent leakage area (C value), cm2/m2
	        story: story
	        inside_pressure: inside pressure against outdoor pressure
	            'negative': negative pressure
	            'positive': positive pressure
	            'balanced': balanced
	    Returns:
	        air leakage volume at step n, m3/s, [i,1]
	    Note:
	        eq.3
*/
func _get_n_leak_n(
	c_value float64,
	story Story,
	inside_pressure InsidePressure,
	delta_theta_n float64,
) float64 {
	var a float64 // 係数aの計算, 回/(h (cm2/m2 K^0.5))

	switch story {
	case StoryOne:
		// 1階建ての時の係数
		a = 0.022
	case StoryTwo:
		// 2階建ての時の係数
		a = 0.020
	default:
		panic(story)
	}

	// 係数bの計算, 回/h
	// 階数と換気方式の組み合わせで決定する
	var b float64
	switch inside_pressure {
	case InsidePressureBalanced:
		switch story {
		case StoryOne:
			b = 0.00
		case StoryTwo:
			b = 0.0
		default:
			panic(story)
		}
	case InsidePressurePositive:
		switch story {
		case StoryOne:
			b = 0.26
		case StoryTwo:
			b = 0.14
		default:
			panic(story)
		}
	case InsidePressureNegative:
		switch story {
		case StoryOne:
			b = 0.28
		case StoryTwo:
			b = 0.13
		default:
			panic(story)
		}
	default:
		panic(inside_pressure)
	}

	// 換気回数の計算
	// Note: 切片bの符号は-が正解（報告書は間違っている）
	n_leak_n := math.Max(a*(c_value*math.Sqrt(delta_theta_n))-b, 0.0)

	return n_leak_n
}

/*
		Calculate the temperature difference between room and outdoor.

	    Args:
	        theta_average_r_n: averate room temperature at step n, degree C
	        theta_o_n: outdoor temperature at step n, degree C

	    Returns:
	        temperature difference between room and outdoor at step n, K

	    Notes:
	        eq.4
*/
func _get_delta_theta_n(theta_average_r_n float64, theta_o_n float64) float64 {

	delta_theta_n := math.Abs(theta_average_r_n - theta_o_n)

	return delta_theta_n
}

/*
		Calculate the average air temperature at step n which is weghted by room volumes.

	    Args:
	        theta_r_is_n: room temperature of room i in step n, degree C, [i, 1]
	        v_rm_is: room volume of room i, m3, [i, 1]

	    Returns:
	        average air temperature at step n, degree C

	    Note:
	        eq.5
*/
func _get_theta_average_r_n(theta_r_is_n []float64, v_rm_is []float64) float64 {
	return stat.Mean(theta_r_is_n, v_rm_is)
}
