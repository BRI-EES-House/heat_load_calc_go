package heat_load_calc

import "math"

type EarthenFloor struct {
	psi   float64   //線熱貫流率, W/m K
	rfa0  float64   //吸熱応答係数の初項, -
	rft0  float64   //貫流応答係数の初項, -
	rfa1  []float64 //指数項別吸熱応答係数, -
	rft1  []float64 //指数項別貫流応答係数, -
	nRoot int       //根の数
}

/*
	土間床外周部の貫流応答特性のパラメータを返す
	Returns:
		土間床外周部のパラメータ
*/
func get_rf_parameters_t() []float64 {
	return []float64{
		-0.0374209855594571,
		0.0941209125618409,
		-0.136869498520281,
		-0.493811564418963,
		0.0541164032804825,
		-0.415684492437867,
		0.198293349862905,
		-0.987302200526257,
		1.61323086601284,
		-0.320547647481877,
	}
}

/*
	土間床外周部の吸熱応答特性のパラメータを返す
	Returns:
		土間床外周部のパラメータ
*/
func get_rf_parameters_a() []float64 {
	return []float64{
		-0.000375432069330334,
		0.00585684299109825,
		-0.054548410348266,
		-0.613447717287527,
		0.212610151690178,
		-0.922186268181334,
		0.790978785668312,
		-1.26629291125237,
		1.26609505908759,
		-4.0149975814435,
	}
}

/*
	土間床外周部の応答係数の初項を計算する
	Args:
		parameters: 応答特性のパラメータ
		alpha_m: 固定根数列

	Returns:
		土間床外周部の応答係数の初項
*/
func rf_initial_term(parameters []float64, alpha_m []float64) float64 {
	sum := 0.0
	for i := 0; i < len(parameters); i++ {
		sum += parameters[i] / (alpha_m[i] * 900) * (1.0 - math.Exp(-alpha_m[i]*900))
	}
	return 1.0 + sum
}

/*
	土間床外周部の指数項別応答係数を計算する
	Args:
		parameters: 応答特性のパラメータ
		alpha_m: 固定根数列

	Returns:
		土間床外周部の指数項別応答
*/
func calc_rf_exponential(parameters []float64, alpha_m []float64) []float64 {
	rf_e := make([]float64, len(alpha_m))
	for i := 0; i < len(alpha_m); i++ {
		term := 1.0 - math.Exp(-alpha_m[i]*900.0)
		rf_e[i] = -parameters[i] / (alpha_m[i] * 900.0) * term * term
	}
	return rf_e
}

func NewEarthFloor(psi float64) *EarthenFloor {

	var ef EarthenFloor

	// 線熱貫流率
	ef.psi = psi

	// 貫流応答、吸熱応答のパラメータ取得
	parameters_t := get_rf_parameters_t()
	parameters_a := get_rf_parameters_a()

	// 土壌の根を取得
	alpha_m := get_alpha_m(true)

	// 根の数
	ef.nRoot = len(alpha_m)

	// 吸熱応答の初項
	ef.rfa0 = psi * rf_initial_term(parameters_a, alpha_m)

	// 貫流応答の初項
	ef.rft0 = psi * rf_initial_term(parameters_t, alpha_m)

	// 指数項別吸熱応答係数
	var rfa1, rft1 []float64
	rfa1 = calc_rf_exponential(parameters_a, alpha_m)
	rft1 = calc_rf_exponential(parameters_t, alpha_m)
	for i := 0; i < len(rfa1); i++ {
		rfa1[i] *= psi
	}
	for i := 0; i < len(rft1); i++ {
		rft1[i] *= psi
	}
	ef.rfa1 = rfa1
	ef.rft1 = rft1

	return &ef
}
