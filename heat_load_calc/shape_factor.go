package main

// ***** 室内の境界の形態係数および放射熱伝達率 *****
// See: https://hc-energy.readthedocs.io/ja/latest/contents/03_14_eval_shape_factor.html

import (
	"fmt"
	"log"
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
室 i の微小球に対する境界 j の重み係数を求める

    Args:
        a_s_js: 境界 j の面積, m2, [j, 1]
        h_s_r_js: 境界 j の室内側放射熱伝達率, W/m2K, [j, 1]
        p_is_js: 室 i と境界 j の関係を表す係数（境界から室への変換）, [i, j]

    Returns:
        室 i の微小球に対する境界 j の重み係数, -, [i, j]

    Notes:
        式(1)
*/
func get_f_mrt_is_js(a_s_js, h_s_r_js mat.Vector, p_is_js mat.Matrix) *mat.Dense {
	var ah mat.VecDense
	ah.MulElemVec(a_s_js, h_s_r_js)
	_ah := ah.RawVector().Data

	var term1 mat.Dense
	term1.Apply(func(i, j int, v float64) float64 {
		return v * _ah[j]
	}, p_is_js)

	var term2 mat.VecDense
	term2.MulVec(p_is_js, &ah)
	_term2 := term2.RawVector().Data

	term1.Apply(func(i, j int, v float64) float64 {
		return v / _term2[i]
	}, &term1)

	return &term1
}

/*
室iに接する境界jの放射熱伝達率を求める

Notes:
	式(2)
*/
func get_h_s_r_js(id_rm_is []int, a_s_js *mat.VecDense, connected_room_id_js []int) []float64 {
	h_s_r_js := make([]float64, a_s_js.Len())

	for i := 0; i < len(id_rm_is); i++ {

		a_s_js_filtered := make([]float64, 0)

		for j, room_id := range connected_room_id_js {
			is_connected := room_id == id_rm_is[i]
			if is_connected {
				a_s_js_filtered = append(a_s_js_filtered, a_s_js.AtVec(j))
			}
		}

		temp := _calc_h_s_r_i_js(mat.NewVecDense(len(a_s_js_filtered), a_s_js_filtered))

		jj := 0
		for j, room_id := range connected_room_id_js {
			is_connected := room_id == id_rm_is[i]
			if is_connected {
				h_s_r_js[j] = temp[jj]
				jj++
			}
		}
	}

	return h_s_r_js
}

/*
	放射熱伝達率（室単位で計算する）

    Args:
        a_s_i_js: 境界jの面積, m2, [j, 1]
    Returns:
        放射熱伝達率, W/m2K, [j, 1]
*/
func _calc_h_s_r_i_js(a_s_i_js *mat.VecDense) []float64 {
	// 面積比, [j]
	r_a_i_js := mat.NewVecDense(a_s_i_js.Len(), nil)
	sum_a := sum(a_s_i_js.RawVector().Data)
	r_a_i_js.ScaleVec(1.0/sum_a, a_s_i_js)

	// 非線形方程式L(f̅)=0の解, float
	f_ver := _get_f_ver(r_a_i_js)

	// 放射伝熱計算で使用する微小球に対する部位の形態係数, -, [j]
	f_i_js := _get_f_i_k(f_ver, r_a_i_js)

	// 総和のチェック
	if math.Abs(sum(f_i_js)-1.0) > 1.0e-3 {
		log.Printf("形態係数の合計値が不正 TotalFF=%.10f", sum(f_i_js))
	}

	h_s_r_i_js := _get_h_s_r_i_js(f_i_js)

	return h_s_r_i_js
}

func sum(arr []float64) float64 {
	sum := 0.0
	for _, val := range arr {
		sum += val
	}
	return sum
}

/*
	室 i に接する境界 j の放射熱伝達率を計算する。

    Args:
        f_i_js: 室 i の微小球からみた境界 j への形態係数

    Returns:
        室 i に接する境界 j の放射熱伝達率, W/m2K, [j]
*/
func _get_h_s_r_i_js(f_i_js []float64) []float64 {

	// 境界間の放射熱伝達率を決定する際、平均放射温度を20℃固定値であるとして計算する。
	const theta_mrt_K3 = (20.0 + 273.15) * (20.0 + 273.15) * (20.0 + 273.15)

	h_s_r_i_js := make([]float64, len(f_i_js))
	for j := range f_i_js {
		h_s_r_i_js[j] = eps / (1.0 - eps*f_i_js[j]) * 4.0 * sgm * theta_mrt_K3
	}
	return h_s_r_i_js
}

/*
   非線形方程式L(f̅)=0の解
   Args:
       r_a_i_k: 面積比, [j]
   Returns:
       非線形方程式L(f̅)=0の解
   Notes:
       式(5)
*/
func _get_f_ver(r_a_i_k *mat.VecDense) float64 {

	f := func(f_ver float64) float64 {
		var sum float64
		for _, f_i := range _get_f_i_k(f_ver, r_a_i_k) {
			sum += f_i
		}
		return sum - 1.0
	}

	// Bisection method
	findRoot := func(a float64, b float64, tol float64, maxIter int) (float64, error) {
		if f(a)*f(b) >= 0 {
			return 0, fmt.Errorf("no root found in the interval [%f, %f]", a, b)
		}

		var c float64
		for i := 0; i < maxIter; i++ {
			c = (a + b) / 2

			if f(c) == 0 || (b-a)/2 < tol {
				return c, nil
			}

			if f(c)*f(a) < 0 {
				b = c
			} else {
				a = c
			}
		}
		return 0, fmt.Errorf("failed to find root within %d iterations", maxIter)
	}

	//NOTE: fsolveはFortranのMINPACK相当だが、Goでは相当の実装が見当たらない
	//ref: https://github.com/fortran-lang/minpack

	a := 0.0
	b := 5.0
	tol := 1e-6
	maxIter := 100

	root, err := findRoot(a, b, tol, maxIter)
	if err != nil {
		panic(err)
	}

	return root
}

/*
   空間内の微小球からみた面iへの形態係数を計算する。
   Args:
       f_ver_i: 非線形方程式L(f_ver_i)の解, -
       r_a_i_k: 同一方位となる表面のグループ k の面積が室 i 内の表面積の総和に占める比, -, [k]
   Returns:
       室 i の微小球から同一方位となる表面のグループ k への形態係数
   Notes:
       式(4)
*/
func _get_f_i_k(f_ver_i float64, r_a_i_k *mat.VecDense) []float64 {
	result := make([]float64, r_a_i_k.Len())
	for i, r := range r_a_i_k.RawVector().Data {
		term := 1 - 4*r/f_ver_i
		if term < 0 {
			result[i] = 0.5 * (1.0 + math.Sqrt(math.Abs(term)))
		} else {
			result[i] = 0.5 * (1.0 - math.Sqrt(math.Abs(term)))
		}
	}
	return result
}
