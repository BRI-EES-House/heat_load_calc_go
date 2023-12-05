package heat_load_calc

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/gocarina/gocsv"
	"gonum.org/v1/gonum/mat"
)

type Weather struct {
	_a_sun_ns       []float64               // ステップ n における太陽高度, rad, [N]
	_h_sun_ns       []float64               // ステップ n における太陽方位角, rad, [N]
	_i_dn_ns        []float64               // ステップ n における法線面直達日射量, W/m2, [n]
	_i_sky_ns       []float64               // ステップ n における水平面天空日射量, W/m2, [N]
	_r_n_ns         []float64               // ステップ n における夜間放射量, W/m2, [N]
	_theta_o_ns     []float64               // ステップ n における外気温度, degree C, [N]
	_x_o_ns         []float64               // ステップ n における外気絶対湿度, kg/kg(DA), [N]
	_itv            Interval                // 時間間隔
	a_sun_ns_plus   []float64               // ステップ n における太陽高度, rad, [N+1] (先頭用途と最後尾の要素の値は同じ)
	h_sun_ns_plus   []float64               // ステップ n における太陽方位角, rad, [N+1] (先頭用途と最後尾の要素の値は同じ)
	i_dn_ns_plus    []float64               // ステップ n における法線面直達日射量, W/m2, [N+1] (先頭用途と最後尾の要素の値は同じ)
	i_sky_ns_plus   []float64               // ステップ n における水平面天空日射量, W/m2, [N+1] (先頭用途と最後尾の要素の値は同じ)
	r_n_ns_plus     []float64               // ステップ n における夜間放射量, W/m2, [N+1] (先頭用途と最後尾の要素の値は同じ)
	theta_o_ns_plus []float64               // ステップ n における外気温度, degree C, [N+1] (先頭用途と最後尾の要素の値は同じ)
	x_o_ns_plus     *mat.VecDense           // ステップ n における外気温度, degree C, [N+1] (先頭用途と最後尾の要素の値は同じ)
	cos_phi_j_ns    map[Direction][]float64 // ステップ n における方位 dir の傾斜面に入射する太陽の入射角の余弦, -,  [N+1] (先頭用途と最後尾の要素の値は同じ)
}

/*
Args
	a_sun_ns 外気温度, degree C, [n]
	h_sun_ns 外気絶対湿度, kg/kg(DA), [n]
	i_dn_ns 法線面直達日射量, W/m2, [n]
	i_sky_ns 水平面天空日射量, W/m2, [n]
	r_n_ns 夜間放射量, W/m2, [n]
	theta_o_ns 太陽高度, rad, [n]
	x_o_ns 太陽方位角, rad, [n]
	itv 時間間隔
*/
func NewWeather(
	a_sun_ns, h_sun_ns, i_dn_ns, i_sky_ns, r_n_ns, theta_o_ns, x_o_ns []float64,
	itv Interval,
) *Weather {
	w := &Weather{
		_a_sun_ns:       a_sun_ns,
		_h_sun_ns:       h_sun_ns,
		_i_dn_ns:        i_dn_ns,
		_i_sky_ns:       i_sky_ns,
		_r_n_ns:         r_n_ns,
		_theta_o_ns:     theta_o_ns,
		_x_o_ns:         x_o_ns,
		_itv:            itv,
		a_sun_ns_plus:   _add_index_0_data_to_end_slice(a_sun_ns),
		h_sun_ns_plus:   _add_index_0_data_to_end_slice(h_sun_ns),
		i_dn_ns_plus:    _add_index_0_data_to_end_slice(i_dn_ns),
		i_sky_ns_plus:   _add_index_0_data_to_end_slice(i_sky_ns),
		r_n_ns_plus:     _add_index_0_data_to_end_slice(r_n_ns),
		theta_o_ns_plus: _add_index_0_data_to_end_slice(theta_o_ns),
		x_o_ns_plus:     _add_index_0_data_to_end(x_o_ns),
	}

	w.cos_phi_j_ns = make(map[Direction][]float64)
	for _, d := range []Direction{DirectionS, DirectionSW, DirectionW, DirectionN, DirectionNE, DirectionE, DirectionSE, DirectionTop, DirectionBottom} {
		phi_j_ns := get_phi_j_ns(w.h_sun_ns_plus, w.a_sun_ns_plus, d)

		cos_phi := make([]float64, len(phi_j_ns))
		for i, phi_j := range phi_j_ns {
			cos_phi[i] = math.Cos(phi_j)
		}
		w.cos_phi_j_ns[d] = cos_phi
	}

	return w
}

func make_weather(method string, itv Interval, file_path string, region Region) *Weather {
	if method == "file" {
		log.Printf("Load weather data from `%s`", file_path)
		return make_from_pd(file_path, itv)
	} else if method == "ees" {
		log.Printf("make weather data based on the EES region")
		return _make_weather_ees(region, itv)
	} else {
		panic(method)
	}
}

func (self *Weather) get_weather_csv() string {
	var sb strings.Builder

	// Header
	sb.WriteString(",temperature,absolute humidity,normal direct solar radiation,horizontal sky solar radiation,outward radiation,sun altitude,sun azimuth\n")

	// インターバル(分)
	freq := self._itv.get_pandas_freq()

	// 開始日時と終了日時を定義
	start := time.Date(1989, 1, 1, 0, 0, 0, 0, time.Local)
	end := time.Date(1989, 12, 31, 23, 59, 59, 0, time.Local)

	// 現在の日時を開始日時にセット
	current := start

	// 終了日時に到達するまでループ
	n := 0
	for current.Before(end) {
		sb.WriteString(fmt.Sprintf("%s,%g,%g,%g,%g,%g,%g,%g\n",
			current.Format("2006-01-02 15:04:05"),
			self._theta_o_ns[n],
			self._x_o_ns[n],
			self._i_dn_ns[n],
			self._i_sky_ns[n],
			self._r_n_ns[n],
			self._h_sun_ns[n],
			self._a_sun_ns[n],
		))

		// 30分を加算
		current = current.Add(time.Duration(freq * float64(time.Minute)))
		n++
	}

	return sb.String()
}

// データの数を取得する。
func (self *Weather) number_of_data() int {
	return self._itv.get_n_hour() * 8760
}

/*
データの数に1を足したものを取得する。
例えば、1時間間隔の場合、データの数は8760なので、返す値は8761
15分間隔の場合、データの数は35040なので、返す値は35041
Returns
	データの数に1を足したもの
*/
func (self *Weather) number_of_data_plus() int {
	return self.number_of_data() + 1
}

/*
外気温度の年間平均値を取得する。
Returns
	外気温度の年間平均値, degree C
*/
func (self *Weather) get_theta_o_ave() float64 {
	var avg float64
	for _, v := range self._theta_o_ns {
		avg += v
	}
	avg /= float64(len(self._theta_o_ns))
	return avg
}

type WeatherDataRow struct {
	longitude                      string  `csv:"longitude"`
	latitude                       string  `csv:"latitude"`
	temperature                    float64 `csv:"temperature"`
	absolute_humidity              float64 `csv:"absolute_humidity"`
	normal_direct_solar_radiation  float64 `csv:"normal_direct_solar_radiation"`
	horizontal_sky_solar_radiation float64 `csv:"horizontal_sky_solar_radiation"`
	outward_radiation              float64 `csv:"outward_radiation"`
	sun_altitude                   float64 `csv:"sun_altitude"`
	sun_azimuth                    float64 `csv:"sun_azimuth"`
}

/*
気象データを読み込む。

Args
	file_path 気象データのファイルのパス
	itv Interval 列挙体
Returns
	OutdoorCondition クラス
*/
func make_from_pd(file_path string, itv Interval) *Weather {

	// file is exist
	if _, err := os.Stat(file_path); os.IsNotExist(err) {
		panic(fmt.Sprintf("Error File %s is not exist.", file_path))
	}

	// Open the CSV file
	file, err := os.Open(file_path)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create a slice of MyStruct to store the CSV data
	var pp []*WeatherDataRow

	// Unmarshal the CSV data into the slice of MyStruct
	if err := gocsv.UnmarshalFile(file, &pp); err != nil {
		panic(err)
	}

	if len(pp) != 8760 {
		panic("Error Row length of the file should be 8760.")
	}

	latitude, err := strconv.ParseFloat(pp[0].longitude, 64)
	if err != nil {
		panic(err)
	}
	longitude, err := strconv.ParseFloat(pp[0].latitude, 64)
	if err != nil {
		panic(err)
	}

	phi_loc, lambda_loc := math.Pi/180*latitude, math.Pi/180*longitude

	// 太陽位置
	//   (1) ステップ n における太陽高度, rad, [n]
	//   (2) ステップ n における太陽方位角, rad, [n]
	h_sun_ns, a_sun_ns := calc_solar_position(phi_loc, lambda_loc, itv)

	f := func(getc func(row *WeatherDataRow) float64) []float64 {
		var ret []float64 = make([]float64, len(pp))
		for i := range pp {
			ret[i] = getc(pp[i])
		}
		return ret
	}

	// 外気温度, degree C
	theta_o_ns := _interpolate(
		f(func(row *WeatherDataRow) float64 {
			return row.temperature
		}),
		itv,
		false,
	)

	// 外気絶対湿度, kg/kg(DA)
	// g/kgDA から kg/kgDA へ単位変換を行う。
	x_o_ns := _interpolate(
		f(func(row *WeatherDataRow) float64 {
			return row.absolute_humidity
		}),
		itv,
		false,
	)
	for i, v := range x_o_ns {
		x_o_ns[i] = v / 1000.0
	}

	// 法線面直達日射量, W/m2
	i_dn_ns := _interpolate(
		f(func(row *WeatherDataRow) float64 {
			return row.normal_direct_solar_radiation
		}),
		itv,
		false,
	)

	// 水平面天空日射量, W/m2
	i_sky_ns := _interpolate(
		f(func(row *WeatherDataRow) float64 {
			return row.horizontal_sky_solar_radiation
		}),
		itv,
		false,
	)

	// 夜間放射量, W/m2
	r_n_ns := _interpolate(
		f(func(row *WeatherDataRow) float64 {
			return row.outward_radiation
		}),
		itv,
		false,
	)

	return NewWeather(
		a_sun_ns,
		h_sun_ns,
		i_dn_ns,
		i_sky_ns,
		r_n_ns,
		theta_o_ns,
		x_o_ns,
		itv,
	)
}

/*
リストの最後に一番最初のデータを追加する。

Args
	d リスト

Returns
	追加されたリスト
*/
func _add_index_0_data_to_end(d []float64) *mat.VecDense {
	ret := make([]float64, len(d)+1)
	copy(ret, d)
	ret[len(d)] = ret[0]
	return mat.NewVecDense(len(ret), ret)
}

/*
リストの最後に一番最初のデータを追加する。

Args
	d リスト

Returns
	追加されたリスト
*/
func _add_index_0_data_to_end_slice(d []float64) []float64 {
	ret := make([]float64, len(d)+1)
	copy(ret, d)
	ret[len(d)] = ret[0]
	return ret
}

/*
気象データを作成する。
Args
	rgn 地域の区分
	itv Interval 列挙体
Returns
	OutdoorCondition クラス
*/
func _make_weather_ees(rgn Region, itv Interval) *Weather {
	// 気象データの読み込み
	//   (1)ステップnにおける外気温度, degree C, [n]
	//   (2)ステップnにおける法線面直達日射量, W/m2, [n]
	//   (3)ステップnにおける水平面天空日射量, W/m2, [n]
	//   (4)ステップnにおける夜間放射量, W/m2, [n]
	//   (5)ステップnにおける外気絶対湿度, kg/kgDA, [n]
	// インターバルごとの要素数について
	//   interval = "1h" -> n = 8760
	//   interval = "30m" -> n = 8760 * 2
	//   interval = "15m" -> n = 8760 * 4
	theta_o_ns, i_dn_ns, i_sky_ns, r_n_ns, x_o_ns := _load(rgn, itv)

	// 緯度, rad & 経度, rad
	phi_loc, lambda_loc := rgn.get_phi_loc_and_lambda_loc()

	// 太陽位置
	//   (1) ステップ n における太陽高度, rad, [n]
	//   (2) ステップ n における太陽方位角, rad, [n]
	h_sun_ns, a_sun_ns := calc_solar_position(phi_loc, lambda_loc, itv)

	return NewWeather(
		a_sun_ns,
		h_sun_ns,
		i_dn_ns,
		i_sky_ns,
		r_n_ns,
		roundFloat64Slice(theta_o_ns, 3),
		roundFloat64Slice(x_o_ns, 6),
		itv,
	)
}

func roundToDecimalPlaces(value float64, decimalPlaces int) float64 {
	shift := math.Pow10(decimalPlaces)
	return math.Round(value*shift) / shift
}

func roundFloat64Slice(slice []float64, decimalPlaces int) []float64 {
	roundedSlice := make([]float64, len(slice))
	for i, value := range slice {
		roundedSlice[i] = roundToDecimalPlaces(value, decimalPlaces)
	}
	return roundedSlice
}

/*
地域の区分に応じて気象データを読み込み、指定された時間間隔で必要に応じて補間を行いデータを作成する。

Args
	rgn 地域の区分
	itv Interval 列挙体

Returns
	以下の5項目
		(1) ステップnにおける外気温度, degree C [n]
		(2) ステップnにおける法線面直達日射量, W/m2 [n]
		(3) ステップnにおける水平面天空日射量, W/m2 [n]
		(4) ステップnにおける夜間放射量, W/m2 [n]
		(5) ステップnにおける外気絶対湿度, g/kgDA [n]

Notes
	interval = "1h" -> n = 8760
	interval = "30m" -> n = 8760 * 2
	interval = "15m" -> n = 8760 * 4
*/
func _load(rgn Region, itv Interval) ([]float64, []float64, []float64, []float64, []float64) {

	// 地域の区分に応じたファイル名の取得
	weather_data_filename := _get_filename(rgn)

	// ファイル読み込み
	path_and_filename := filepath.Join("expanded_amedas/", weather_data_filename)

	// Open the CSV file
	file, err := os.Open(path_and_filename)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	// Create a CSV reader and read the file
	reader := csv.NewReader(file)

	// Skip the first two header lines
	_, _ = reader.Read()
	_, _ = reader.Read()

	// Read all records from the CSV file
	records, err := reader.ReadAll()
	if err != nil {
		log.Fatal(err)
	}

	// データの整形
	usecols := []int{2, 3, 4, 5, 6} // 外気温度, 法線面直達日射量, 水平面天空日射量, 夜間放射量, 外気絶対湿度
	weather := make([][]float64, len(usecols))
	for i, col := range usecols {
		weather[i] = make([]float64, len(records))
		for j, record := range records {
			weather[i][j], _ = strconv.ParseFloat(record[col], 64)
		}
	}

	// ステップnにおける外気温度, degree C
	theta_o_ns := _interpolate(weather[0], itv, true)

	// ステップnにおける法線面直達日射量, W/m2
	i_dn_ns := _interpolate(weather[1], itv, true)

	// ステップnにおける水平面天空日射量, W/m2
	i_sky_ns := _interpolate(weather[2], itv, true)

	// ステップnにおける夜間放射量, W/m2
	r_n_ns := _interpolate(weather[3], itv, true)

	// ステップnにおける外気絶対湿度, kg/kgDA
	// g/kgDA から kg/kgDA へ単位変換を行う。
	x_o_ns := _interpolate(weather[4], itv, true)
	for i := range x_o_ns {
		x_o_ns[i] /= 1000.0
	}

	return theta_o_ns, i_dn_ns, i_sky_ns, r_n_ns, x_o_ns
}

/*
1時間ごとの8760データを指定された間隔のデータに補間する。
"1h" 1時間間隔の場合、 n = 8760
"30m" 30分間隔の場合、 n = 8760 * 2 = 17520
"15m" 15分間隔の場合、 n = 8760 * 4 = 35040

Args
	weather_data 1時間ごとの気象データ [8760]
	interval 生成するデータの時間間隔
	rolling rolling するか否か。データが1時始まりの場合は最終行の 12/31 2400 のデータを 1/1 000 に持ってくるため、この値は True にすること。

Returns
	指定する時間間隔に補間された気象データ [n]
*/
func _interpolate(weather_data []float64, interval Interval, rolling bool) []float64 {
	if interval == IntervalH1 {

		if rolling {
			// 拡張アメダスのデータが1月1日の1時から始まっているため1時間ずらして0時始まりのデータに修正する。
			return roll(weather_data, 1)
		} else {
			return weather_data
		}
	} else {
		// 補間比率の係数
		alpha := map[Interval][]float64{
			IntervalM30: {1.0, 0.5},
			IntervalM15: {1.0, 0.75, 0.5, 0.25},
		}[interval]

		// 補間元データ1, 補間元データ2
		var data1, data2 []float64
		if rolling {
			// 拡張アメダスのデータが1月1日の1時から始まっているため1時間ずらして0時始まりのデータに修正する。
			data1 = roll(weather_data, 1) // 0時=24時のため、1回分前のデータを参照
			data2 = weather_data
		} else {
			data1 = weather_data
			data2 = roll(weather_data, -1)
		}

		ndata := len(data1)                  // 補間元のデータ数
		nalpha := len(alpha)                 // 補間比率の係数の数
		n := len(data1) * nalpha             // 補間後のデータ数
		data_interp_1d := make([]float64, n) // 補間後のデータ
		off := 0
		for i := 0; i < ndata; i++ {
			for j := 0; j < nalpha; j++ {
				data_interp_1d[off] = alpha[j]*data1[i] + (1.0-alpha[j])*data2[i]
				off++
			}
		}

		return data_interp_1d
	}
}

func roll(slice []float64, shift int) []float64 {
	length := len(slice)
	shift %= length
	if shift < 0 {
		shift += length
	}
	result := append(slice[length-shift:], slice[:length-shift]...)
	return result
}

/*
地域の区分に応じたファイル名を取得する。

Args
	rgn 地域の区分

Returns
	地域の区分に応じたファイル名（CSVファイル）（拡張子も含む）
*/
func _get_filename(rgn Region) string {
	weather_data_filename := map[Region]string{
		Region1: "01_kitami.csv",     // 1地域（北見）
		Region2: "02_iwamizawa.csv",  // 2地域（岩見沢）
		Region3: "03_morioka.csv",    // 3地域（盛岡）
		Region4: "04_nagano.csv",     // 4地域（長野）
		Region5: "05_utsunomiya.csv", // 5地域（宇都宮）
		Region6: "06_okayama.csv",    // 6地域（岡山）
		Region7: "07_miyazaki.csv",   // 7地域（宮崎）
		Region8: "08_naha.csv",       // 8地域（那覇）
	}[rgn]

	return weather_data_filename
}
