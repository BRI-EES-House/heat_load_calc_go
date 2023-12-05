package heat_load_calc

// **** 建物全般のパラメータ ****
// https://hc-energy.readthedocs.io/ja/latest/contents/03_10_eval_building.html

import (
	"math"

	"gonum.org/v1/gonum/mat"
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
	C値を求める。
	Args:
		ua_value: UA値, W/m2 K
		struct: 構造
	Returns:
		C値, cm2/m2
	Notes:
		式(1)
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

var __v_leak_is_n mat.VecDense

/*
	室 i のすきま風量を求める。

	Args:

		theta_r_is_n: ステップ n における室 i の温度, degree C, [I]
		theta_o_n: ステップ n における外気温, degree C
		v_rm_is: 室 i の容積, m3, [I]

	Returns:
		ステップ n における室 i のすきま風量, m3/s, [I]

	Notes:
		式(2)
*/
func (self *Building) get_v_leak_is_n(
	theta_r_is_n *mat.VecDense,
	theta_o_n float64,
	v_rm_is *mat.VecDense,
) *mat.VecDense {

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
	__v_leak_is_n.ScaleVec(n_leak_n/3600, v_rm_is)

	return &__v_leak_is_n
}

/*
	すきま風による住宅全体の換気回数を求める。

	Args:
		c_value: 相当隙間面積(C値), cm2/m2
		story: 建物の階数
		inside_pressure: inside pressure against outdoor pressure
			'negative': negative pressure
			'positive': positive pressure
			'balanced': balanced
	Returns:
		ステップnのすきま風による住宅全体の換気回数, 1/h, [i,1]
	Note:
		式(3)
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
	室内外温度差を求める。

	Args:
		theta_average_r_n: ステップnにおける平均室温, degree C
		theta_o_n: ステップnにおける外気温, degree C

	Returns:
		ステップnにおける室内外温度差, K

	Notes:
		式(4)
*/
func _get_delta_theta_n(theta_average_r_n float64, theta_o_n float64) float64 {

	delta_theta_n := math.Abs(theta_average_r_n - theta_o_n)

	return delta_theta_n
}

/*
	平均室温を求める。

	Args:
		theta_r_is_n: ステップnにおける室iの温度, degree C, [i, 1]
		v_rm_is: 室iの容積, m3, [i,1]

	Returns:
		ステップnにおける平均室温, degree C

	Note:
		式(5)
*/
func _get_theta_average_r_n(theta_r_is_n mat.Vector, v_rm_is mat.Vector) float64 {
	sumValues := mat.Dot(theta_r_is_n, v_rm_is)
	sumWeights := mat.Sum(v_rm_is)
	return sumValues / sumWeights
}
