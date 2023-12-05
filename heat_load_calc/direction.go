package heat_load_calc

import (
	"math"
)

// 方位
type Direction string

const (
	DirectionS      Direction = "s"
	DirectionSW     Direction = "sw"
	DirectionW      Direction = "w"
	DirectionNW     Direction = "nw"
	DirectionN      Direction = "n"
	DirectionNE     Direction = "ne"
	DirectionE      Direction = "e"
	DirectionSE     Direction = "se"
	DirectionTop    Direction = "top"
	DirectionBottom Direction = "bottom"
)

func DirectionFromString(str string) Direction {
	switch str {
	case "s":
		return DirectionS
	case "sw":
		return DirectionSW
	case "w":
		return DirectionW
	case "nw":
		return DirectionNW
	case "n":
		return DirectionN
	case "ne":
		return DirectionNE
	case "e":
		return DirectionE
	case "se":
		return DirectionSE
	case "top":
		return DirectionTop
	case "bottom":
		return DirectionBottom
	default:
		panic("invalid direction")
	}
}

func (d Direction) alpha_w_j() float64 {
	/*境界 j の傾斜面の方位角を取得する。

	Returns:
		alpha_w_j: 境界 j の傾斜面の方位角, rad
	*/

	if d == DirectionTop || d == DirectionBottom {
		panic("方位が上面・下面が定義されているにもかかわらず、方位角を取得しようとしました。")
	}

	switch d {
	case DirectionS:
		return math.Pi * 0.0 / 180.0
	case DirectionSW:
		return math.Pi * 45.0 / 180.0
	case DirectionW:
		return math.Pi * 90.0 / 180.0
	case DirectionNW:
		return math.Pi * 135.0 / 180.0
	case DirectionN:
		return math.Pi * 180.0 / 180.0
	case DirectionNE:
		return math.Pi * -135.0 / 180.0
	case DirectionE:
		return math.Pi * -90.0 / 180.0
	case DirectionSE:
		return math.Pi * -45.0 / 180.0
	default:
		panic("invalid direction")
	}
}

func (d Direction) beta_w_j() float64 {
	/*境界 j の傾斜面の傾斜角を取得する。

	Returns:
		境界 j の傾斜面の傾斜角, rad
	*/

	if d == DirectionS || d == DirectionSW || d == DirectionW || d == DirectionNW || d == DirectionN || d == DirectionNE || d == DirectionE || d == DirectionSE {
		return math.Pi * 90.0 / 180.0
	} else if d == DirectionTop {
		return math.Pi * 0.0 / 180.0
	} else if d == DirectionBottom {
		return math.Pi * 180.0 / 180.0
	} else {
		panic("invalid direction")
	}
}
