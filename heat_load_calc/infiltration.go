package heat_load_calc

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
   室温と外気温度から隙間風を計算する関数を作成する。
   Args:
       v_rm_is: 空間 i の気積, m3, [i, 1]
       building: Building クラス
   Returns:
       室温と外気温度から隙間風を計算する関数
   Notes:
       作成される関数の引数と戻り値は以下のとおり。
           引数:
               theta_r_is_n: 時刻nの室温, degree C, [i,1]
               theta_o_n: 時刻n+1の外気温度, degree C
           戻り値:
               すきま風量, m3/s, [i,1]
*/
func make_get_infiltration_function(v_rm_is mat.Vector, building *Building) func(theta_r_is_n mat.Vector, theta_o_n float64) mat.Vector {
	return func(theta_r_is_n mat.Vector, theta_o_n float64) mat.Vector {
		return _get_infiltration_residential(
			theta_r_is_n,
			theta_o_n,
			building.c_value,
			v_rm_is,
			building.story,
			building.inside_pressure,
		)
	}
}

/*
	すきま風量を取得する関数（住宅用、圧力バランスを解いた近似式バージョン）
    住宅を１つの空間に見立てて予め圧力バランスを解き、
    Args:
        theta_r_is_n: 時刻nの室温, degree C, [i,1]
        theta_o_n: 時刻n+1の外気温度, degree C
        c_value: 相当隙間面積, cm2/m2
        v_room_is: 室iの容積, m3, [i,1]
        story: 階
        inside_pressure: 室内側の圧力
            'negative': 負圧
            'positive': 正圧
            'balanced': バランス圧力
    Returns:
        すきま風量, m3/s, [i,1]
*/
func _get_infiltration_residential(
	theta_r_is_n mat.Vector,
	theta_o_n float64,
	c_value float64,
	v_room_is mat.Vector,
	story Story,
	inside_pressure InsidePressure,
) mat.Vector {
	// 室気積加重平均室温theta_r_nの計算, degree C, float
	theta_average_r_n := average(theta_r_is_n, v_room_is)

	// 室内外温度差の計算, degree C, float
	deltaTheta := math.Abs(theta_average_r_n - theta_o_n)

	// 係数aの計算, 回/(h (cm2/m2 K^0.5))
	a := map[Story]float64{
		StoryOne: 0.022, // 1階建ての時の係数
		StoryTwo: 0.020, // 2階建ての時の係数
	}[story]

	// 係数bの計算, 回/h
	// 階数と換気方式の組み合わせで決定する
	b := map[InsidePressure]map[Story]float64{
		InsidePressureBalanced: {
			StoryOne: 0.00,
			StoryTwo: 0.00,
		},
		InsidePressurePositive: {
			StoryOne: 0.26,
			StoryTwo: 0.14,
		},
		InsidePressureNegative: {
			StoryOne: 0.28,
			StoryTwo: 0.13,
		},
	}[inside_pressure]

	// 換気回数の計算
	// Note: 切片bの符号は-が正解（報告書は間違っている）
	infiltration_rate := math.Max(a*(c_value*math.Sqrt(deltaTheta))-b[story], 0)

	// すきま風量の計算
	infiltration := make([]float64, v_room_is.Len())
	for i := 0; i < v_room_is.Len(); i++ {
		infiltration[i] = infiltration_rate * v_room_is.AtVec(i) / 3600.0
	}

	return mat.NewVecDense(len(infiltration), infiltration)
}

func average(values mat.Vector, weights mat.Vector) float64 {
	sumValues := 0.0
	sumWeights := 0.0

	for i := 0; i < values.Len(); i++ {
		sumValues += values.AtVec(i) * weights.AtVec(i)
		sumWeights += weights.AtVec(i)
	}

	return sumValues / sumWeights
}
