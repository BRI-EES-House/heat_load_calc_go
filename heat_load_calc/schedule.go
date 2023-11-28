package heat_load_calc

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"path/filepath"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

/*
居住人数の指定方法
1～4人を指定するか、auto の場合は床面積から居住人数を指定する方法を選択する。
*/
type NumberOfOccupants int

// 居住人数の指定方法
const (
	One NumberOfOccupants = iota + 1
	Two
	Three
	Four
	Auto
)

// 居住人数の指定方法を文字列で返す。
func (n NumberOfOccupants) String() string {
	return [...]string{"1", "2", "3", "4", "auto"}[n-1]
}

type Schedule struct {
	q_gen_is_ns            *ScheduleData // ステップ　n　の室　i　における内部発熱, W, [i, n]
	x_gen_is_ns            *ScheduleData // ステップ　n　の室　i　における人体発湿を除く内部発湿, kg/s, [i, n]
	v_mec_vent_local_is_ns *ScheduleData // ステップ　n　の室　i　における局所換気量, m3/s, [i, n]
	n_hum_is_ns            *ScheduleData // ステップ　n　の室　i　における在室人数, [i, n]
	ac_demand_is_ns        *ScheduleData // ステップ　n　の室　i　における空調需要, [i, n]
	ac_setting_is_ns       *ScheduleData // ステップ n の室 i における空調モード, [i, n]
}

/*
Args:

	q_gen_is_ns: ステップ　n　の室　i　における内部発熱, W, [i, n]
	x_gen_is_ns: ステップ　n　の室　i　における人体発湿を除く内部発湿, kg/s, [i, n]
	v_mec_vent_local_is_ns: ステップ　n　の室　i　における局所換気量, m3/s, [i, n]
	n_hum_is_ns: ステップ　n　の室　i　における在室人数, [i, n]
	ac_demand_is_ns: ステップ　n　の室　i　における空調需要, [i, n]
	ac_setting_is_ns: ステップ n の室 i における空調モード, [i, n]
*/
func NewSchedule(
	q_gen_is_ns *ScheduleData,
	x_gen_is_ns *ScheduleData,
	v_mec_vent_local_is_ns *ScheduleData,
	n_hum_is_ns *ScheduleData,
	ac_demand_is_ns *ScheduleData,
	ac_setting_is_ns *ScheduleData,
) *Schedule {
	return &Schedule{
		q_gen_is_ns:            q_gen_is_ns,
		x_gen_is_ns:            x_gen_is_ns,
		v_mec_vent_local_is_ns: v_mec_vent_local_is_ns,
		n_hum_is_ns:            n_hum_is_ns,
		ac_demand_is_ns:        ac_demand_is_ns,
		ac_setting_is_ns:       ac_setting_is_ns,
	}
}

/*
Schedule クラスを生成する。

Args:

	number_of_occupants: 居住人数の指定方法
	s_name_is: 室 i のスケジュールの名称, [i]
	a_floor_is: 室 i の床面積, m2, [i]

Returns:

	Schedule クラス
*/
func get_schedule(number_of_occupants NumberOfOccupants, s_name_is []string, a_floor_is []float64) *Schedule {

	// 居住人数の指定モード
	noo := number_of_occupants

	// 居住人数
	n_p := _get_n_p(noo, a_floor_is)

	// ステップ n の室 i における局所換気量, m3/s, [i, n]
	// jsonファイルでは、 m3/h で示されているため、単位換算(m3/h -> m3/s)を行っている。
	v_mec_vent_local_is_ns := _get_schedules(s_name_is, noo, n_p, "local_vent_amount", true, false)
	floats.Scale(1.0/3600.0, v_mec_vent_local_is_ns.Data)

	// ステップ n の室 i における機器発熱, W, [i, n]
	q_gen_app_is_ns := _get_schedules(s_name_is, noo, n_p, "heat_generation_appliances", true, false)

	// ステップ n の室 i における調理発熱, W, [i, n]
	q_gen_ckg_is_ns := _get_schedules(s_name_is, noo, n_p, "heat_generation_cooking", true, false)

	// ステップ n の室 i における調理発湿, kg/s, [i, n]
	// jsonファイルでは、g/h で示されているため、単位換算(g/h->kg/s)を行っている。
	x_gen_ckg_is_ns := _get_schedules(s_name_is, noo, n_p, "vapor_generation_cooking", true, false)
	floats.Scale(1.0/1000.0/3600.0, x_gen_ckg_is_ns.Data)

	// ステップ n の室 i における照明発熱, W/m2, [i, n]
	// 単位面積あたりで示されていることに注意
	q_gen_lght_is_ns := _get_schedules(s_name_is, noo, n_p, "heat_generation_lighting", true, false)

	// ステップ n の室 i における在室人数, [i, n]
	// 居住人数で按分しているため、整数ではなく小数であることに注意
	n_hum_is_ns := _get_schedules(s_name_is, noo, n_p, "number_of_people", true, false)

	// ステップ n の室 i における空調割合, [i, n]
	ac_demand_is_ns := _get_schedules(s_name_is, noo, n_p, "is_temp_limit_set", true, true)

	// ステップ n の室 i における空調モード, [i, n]
	ac_setting_is_ns := _get_schedules(s_name_is, noo, n_p, "is_temp_limit_set", false, false)

	// ステップ n の室 i における人体発熱を除く内部発熱, W, [i, n]
	n_rm := len(s_name_is)
	n_step := len(q_gen_app_is_ns.Data) / n_rm
	q_gen_is_ns := make([]float64, n_rm*n_step)
	off := 0
	for j := 0; j < n_step; j++ {
		_q_gen_app_is_ns := q_gen_app_is_ns.Get(j)
		_q_gen_ckg_is_ns := q_gen_ckg_is_ns.Get(j)
		_q_gen_lght_is_ns := q_gen_lght_is_ns.Get(j)
		for i := 0; i < n_rm; i++ {
			q_gen_is_ns[off] = _q_gen_app_is_ns.AtVec(i) + _q_gen_ckg_is_ns.AtVec(i) + _q_gen_lght_is_ns.AtVec(i)*a_floor_is[i]
			off++
		}
	}

	// ステップ n の室 i における人体発湿を除く内部発湿, kg/s, [i, n]
	x_gen_is_ns := x_gen_ckg_is_ns

	return NewSchedule(
		&ScheduleData{Data: q_gen_is_ns, BatchSize: n_rm},
		x_gen_is_ns,
		v_mec_vent_local_is_ns,
		n_hum_is_ns,
		ac_demand_is_ns,
		ac_setting_is_ns,
	)
}

func (self *Schedule) save_schedule(output_data_dir string) {
	// ステップnの室iにおける局所換気量, m3/s, [i, 8760*4]
	mid_data_local_vent_path := filepath.Join(output_data_dir, "mid_data_local_vent.csv")
	log.Printf("Save v_mec_vent_local_is_ns to `%s`", mid_data_local_vent_path)
	self.v_mec_vent_local_is_ns.Save(mid_data_local_vent_path)

	// ステップnの室iにおける内部発熱, W, [8760*4]
	mid_data_heat_generation_path := filepath.Join(output_data_dir, "mid_data_heat_generation.csv")
	log.Printf("Save q_gen_is_ns to `%s`", mid_data_heat_generation_path)
	self.q_gen_is_ns.Save(mid_data_heat_generation_path)

	// ステップnの室iにおける人体発湿を除く内部発湿, kg/s, [8760*4]
	mid_data_moisture_generation_path := filepath.Join(output_data_dir, "mid_data_moisture_generation.csv")
	log.Printf("Save x_gen_is_ns to `%s`", mid_data_moisture_generation_path)
	self.x_gen_is_ns.Save(mid_data_moisture_generation_path)

	// ステップnの室iにおける在室人数, [8760*4]
	mid_data_occupants_path := filepath.Join(output_data_dir, "mid_data_occupants.csv")
	log.Printf("Save n_hum_is_ns to `%s`", mid_data_occupants_path)
	self.n_hum_is_ns.Save(mid_data_occupants_path)

	// ステップnの室iにおける空調需要, [8760*4]
	mid_data_ac_demand_path := filepath.Join(output_data_dir, "mid_data_ac_demand.csv")
	log.Printf("Save ac_demand_is_ns to `%s`", mid_data_ac_demand_path)
	self.ac_demand_is_ns.Save(mid_data_ac_demand_path)
}

/*
365日分のカレンダーを取得する。

Returns:

	365日分のカレンダー
*/
func _load_calendar() []string {
	calender_dict, err := _load_schedule("calendar")
	if err != nil {
		panic(err)
	}

	calobj := calender_dict["calendar"].([]interface{})
	cal := make([]string, len(calobj))
	for i, v := range calobj {
		cal[i] = v.(string)
	}
	return cal
}

// スケジュールを読み込む
func _load_schedule(filename string) (map[string]interface{}, error) {
	filePath := filepath.Join("schedule", filename+".json")
	content, err := ioutil.ReadFile(filePath)
	if err != nil {
		return nil, err
	}

	var jsonData map[string]interface{}
	err = json.Unmarshal(content, &jsonData)
	if err != nil {
		return nil, err
	}

	return jsonData, nil
}

type ScheduleData struct {
	Data      []float64
	BatchSize int
}

func (self *ScheduleData) Get(i int) *mat.VecDense {
	return mat.NewVecDense(self.BatchSize, self.Data[i*self.BatchSize:(i+1)*self.BatchSize])
}

func (self *ScheduleData) Len() int {
	return len(self.Data) / self.BatchSize
}

func (self *ScheduleData) Save(fileName string) error {
	f, err := os.Create(fileName)
	if err != nil {
		log.Fatal(err)
		return err
	}
	defer f.Close()

	r, c := self.Len(), self.BatchSize

	off := 0
	for i := 0; i < r; i++ {
		f.WriteString(_format_float64(self.Data[off]))
		off++
		for j := 1; j < c; j++ {
			f.WriteString(",")
			f.WriteString(_format_float64(self.Data[off]))
			off++
		}
		f.WriteString("\n")
	}

	return nil
}

func _format_float64(f float64) string {
	if f == 1.0 {
		return "1.0"
	} else if f == 0.0 {
		return "0.0"
	}
	return fmt.Sprintf("%g", f)
}

func _get_schedules(
	s_name_is []string,
	noo NumberOfOccupants,
	n_p float64,
	schedule_type string,
	is_proportionable bool,
	is_zero_one bool,
) *ScheduleData {
	data := make([][]float64, len(s_name_is))
	for i, schedule_name_i := range s_name_is {
		data[i] = _get_schedule(
			schedule_name_i,
			noo,
			n_p,
			schedule_type,
			is_proportionable,
			is_zero_one,
		)
	}

	n_rm := len(s_name_is)
	n_step := len(data[0])

	off := 0
	flatdata := make([]float64, n_rm*n_step)
	for i := 0; i < n_step; i++ {
		for j := 0; j < n_rm; j++ {
			flatdata[off] = data[j][i]
			off++
		}
	}

	return &ScheduleData{
		Data:      flatdata,
		BatchSize: n_rm,
	}
}

/*
スケジュールを取得する。
Args:

	schedule_name_i: 室 i のスケジュールの名称
	noo: 居住人数の指定方法（NumberOfOccupants 列挙体）
	n_p: 居住人数
	schedule_type: どのようなスケジュールを扱うのか？　以下から指定する。
		"local_vent_amount"
		"heat_generation_appliances"
		"vapor_generation_cooking"
		"heat_generation_cooking"
		"heat_generation_lighting"
		"number_of_people"
		"is_temp_limit_set"
	is_proportionable: 按分可能かどうか
		按分可能な場合は居住人数により按分が行われる
		按分可能でない場合は2つの数字のうち大きい方の値が採用される
		按分作業が発生しない場合（schedule_type が const の場合または schedule_type が number でかつ居住人数が auto ではない場合）、本パラメータは無視される。
		これが適用されないのは唯一、ac_setting を想定している。
	is_zero_one: 数字データの意味をゼロ・イチの意味に読み替えるかどうか
		例：　[0, 3, 5, 7, 0] -> [0, 1, 1, 1, 0]
		ac_demand に適用されることを想定している

Returns:

	スケジュール, [365*96]
*/
func _get_schedule(
	schedule_name_i string,
	noo NumberOfOccupants,
	n_p float64,
	schedule_type string,
	is_proportionable bool,
	is_zero_one bool,
) []float64 {
	// カレンダーの読み込み（日にちの種類（"平日", "休日外", "休日在"））, [365]
	calendar := _load_calendar()

	// スケジュールを記述した辞書の読み込み
	d, err := _load_schedule(schedule_name_i)
	if err != nil {
		panic(err)
	}

	to_float64_array := func(v interface{}) []float64 {
		vv := v.([]interface{})
		ret := make([]float64, len(vv))
		for i, vvv := range vv {
			ret[i] = vvv.(float64)
		}
		return ret
	}

	convert_schedule := func(day_type string) map[string]interface{} {
		if d["schedule_type"] == "number" {
			return map[string]interface{}{
				"schedule_type": "number",
				"1":             to_float64_array(d["schedule"].(map[string]interface{})["1"].(map[string]interface{})[day_type].(map[string]interface{})[schedule_type].([]interface{})),
				"2":             to_float64_array(d["schedule"].(map[string]interface{})["2"].(map[string]interface{})[day_type].(map[string]interface{})[schedule_type].([]interface{})),
				"3":             to_float64_array(d["schedule"].(map[string]interface{})["3"].(map[string]interface{})[day_type].(map[string]interface{})[schedule_type].([]interface{})),
				"4":             to_float64_array(d["schedule"].(map[string]interface{})["4"].(map[string]interface{})[day_type].(map[string]interface{})[schedule_type].([]interface{})),
			}
		} else if d["schedule_type"] == "const" {
			return map[string]interface{}{
				"schedule_type": "const",
				"const":         to_float64_array(d["schedule"].(map[string]interface{})["const"].(map[string]interface{})[day_type].(map[string]interface{})[schedule_type].([]interface{})),
			}
		} else {
			panic(d["schedule_type"])
		}
	}

	d_weekday := convert_schedule("Weekday")
	d_holiday_out := convert_schedule("Holiday_Out")
	d_holiday_in := convert_schedule("Holiday_In")

	interpolated_schedule_map := map[string][]float64{
		"W":  _get_interpolated_schedule(d_weekday, noo, n_p, is_proportionable, is_zero_one),
		"HO": _get_interpolated_schedule(d_holiday_out, noo, n_p, is_proportionable, is_zero_one),
		"HI": _get_interpolated_schedule(d_holiday_in, noo, n_p, is_proportionable, is_zero_one),
	}

	d_365_96 := make([]float64, 365*96)
	off := 0
	for i := 0; i < 365; i++ {
		c := calendar[i]
		interpolated_schedule := interpolated_schedule_map[c]
		for j := 0; j < 96; j++ {
			d_365_96[off] = interpolated_schedule[j]
			off++
		}
	}

	return d_365_96
}

/*
世帯人数で線形補間してリストを返す
Args:

	daily_schedule: スケジュール
		Keyは必ず"1", "2", "3", "4"
		Valueは96個のリスト形式の値（15分インターバル）
		{
			"schedule_type": "number",
			"1": d["schedule"]["1"][day_type][schedule_type],
			"2": d["schedule"]["2"][day_type][schedule_type],
			"3": d["schedule"]["3"][day_type][schedule_type],
			"4": d["schedule"]["4"][day_type][schedule_type],
		}
		または
		{
			"schedule_type": "const",
			"const": d["schedule"]["const"][day_type][schedule_type]
		}
	noo: 居住人数の指定方法
	n_p: 居住人数
	is_proportionable: 按分可能かどうか
		按分可能な場合は居住人数により按分が行われる
		按分可能でない場合は2つの数字のうち大きい方の値が採用される
		按分作業が発生しない場合（schedule_type が const の場合または schedule_type が number でかつ居住人数が auto ではない場合）、本パラメータは無視される。
		これが適用されないのは唯一、ac_setting を想定している。
	is_zero_one: 数字データの意味をゼロ・イチの意味に読み替えるかどうか
		例：　[0, 3, 5, 7, 0] -> [0, 1, 1, 1, 0]
		ac_demand に適用されることを想定している

Returns:

	線形補間したリスト, [96]
*/
func _get_interpolated_schedule(
	daily_schedule map[string]interface{},
	noo NumberOfOccupants,
	n_p float64,
	is_proportionable bool,
	is_zero_one bool,
) []float64 {
	var scd []float64

	if daily_schedule["schedule_type"].(string) == "const" {

		scd = daily_schedule["const"].([]float64)

		if is_zero_one {
			return convert_to_zero_one(scd)
		} else {
			return scd
		}
	} else if daily_schedule["schedule_type"].(string) == "number" {

		switch noo {
		case One, Two, Three, Four:
			scd = daily_schedule[noo.String()].([]float64)

			if is_zero_one {
				return convert_to_zero_one(scd)
			} else {
				return scd
			}
		case Auto:
			ceil_np, floor_np := _get_ceil_floor_np(n_p)

			var ceil_schedule, floor_schedule, interpolate_np_schedule []float64
			if is_zero_one {
				ceil_schedule = convert_to_zero_one(daily_schedule[fmt.Sprint(ceil_np)].([]float64))
				floor_schedule = convert_to_zero_one(daily_schedule[fmt.Sprint(floor_np)].([]float64))
			} else {
				ceil_schedule = daily_schedule[fmt.Sprint(ceil_np)].([]float64)
				floor_schedule = daily_schedule[fmt.Sprint(floor_np)].([]float64)
			}

			interpolate_np_schedule = make([]float64, len(ceil_schedule))
			if is_proportionable {
				for i := 0; i < len(ceil_schedule); i++ {
					interpolate_np_schedule[i] = ceil_schedule[i]*(n_p-float64(floor_np)) + floor_schedule[i]*(float64(ceil_np)-n_p)
				}
			} else {
				for i := 0; i < len(ceil_schedule); i++ {
					interpolate_np_schedule[i] = math.Max(ceil_schedule[i], floor_schedule[i])
				}
			}

			return interpolate_np_schedule
		default:
			panic(noo)
		}
	} else {
		panic(daily_schedule["schedule_type"])
	}
}

func convert_to_zero_one(scd []float64) []float64 {
	data := make([]float64, len(scd))
	for i := 0; i < len(scd); i++ {
		if scd[i] > 0.0 {
			data[i] = 1.0
		} else {
			data[i] = 0.0
		}
	}
	return data
}

/*
世帯人数から切り上げ・切り下げた人数を整数値で返す

Args:

	n_p: 世帯人数

Returns:

	タプル：
		切り上げた世帯人数
		切り下げた世帯人数

Notes:

	1人未満、4人より大の人数を指定した場合はエラーを返す。
*/
func _get_ceil_floor_np(n_p float64) (int, int) {
	var ceil_np, floor_np int
	switch {
	case 1.0 <= n_p && n_p < 2.0:
		ceil_np = 2
		floor_np = 1
	case 2.0 <= n_p && ceil_np < 3.0:
		ceil_np = 3
		floor_np = 2
	case 3.0 <= n_p && ceil_np <= 4.0:
		ceil_np = 4
		floor_np = 3
	default:
		log.Printf("The number of people is out of range.")
		panic(n_p)
	}
	return ceil_np, floor_np
}

/*
床面積の合計から居住人数を計算する。
Args:

	noo: 居住人数の指定方法
	a_floor_is: 室 i の床面積, m2, [i]

Returns:

	居住人数
*/
func _get_n_p(noo NumberOfOccupants, a_floor_is []float64) float64 {
	switch noo {
	case One:
		return 1.0
	case Two:
		return 2.0
	case Three:
		return 3.0
	case Four:
		return 4.0
	case Auto:
		// 床面積の合計, m2
		var total float64
		for _, aFloor := range a_floor_is {
			total += aFloor
		}

		if total < 30.0 {
			return 1.0
		} else if total < 120.0 {
			return total / 30.0
		} else {
			return 4.0
		}
	default:
		panic("Invalid value for NumberOfOccupants")
	}
}
