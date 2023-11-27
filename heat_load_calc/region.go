package heat_load_calc

import (
	"math"
)

// 地域の区分
type Region string

// 地域の区分の定数
const (
	Region1 Region = "1"
	Region2 Region = "2"
	Region3 Region = "3"
	Region4 Region = "4"
	Region5 Region = "5"
	Region6 Region = "6"
	Region7 Region = "7"
	Region8 Region = "8"
)

/*
地域の区分から緯度、経度を設定する

Returns:
	以下のタプル
		(1) 緯度, rad
		(2) 経度, rad
*/
func (r Region) get_phi_loc_and_lambda_loc() (phi_loc float64, lambda_loc float64) {
	var latitude, longitude float64

	switch r {
	case Region1:
		latitude, longitude = 43.82, 143.91 // 1地域（北見）
	case Region2:
		latitude, longitude = 43.21, 141.79 // 2地域（岩見沢）
	case Region3:
		latitude, longitude = 39.70, 141.17 // 3地域（盛岡）
	case Region4:
		latitude, longitude = 36.66, 138.20 // 4地域（長野）
	case Region5:
		latitude, longitude = 36.55, 139.87 // 5地域（宇都宮）
	case Region6:
		latitude, longitude = 34.66, 133.92 // 6地域（岡山）
	case Region7:
		latitude, longitude = 31.94, 131.42 // 7地域（宮崎）
	case Region8:
		latitude, longitude = 26.21, 127.685 // 8地域（那覇）
	}

	const to_rad = math.Pi / 180
	phi_loc = latitude * to_rad
	lambda_loc = longitude * to_rad

	return
}
