package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"runtime/pprof"
	"time"
)

type InputJson struct {
	Rooms                  []RoomJson                  `json:"rooms"`
	Common                 CommonJson                  `json:"common"`
	Building               BuildingJson                `json:"building"`
	Boundaries             []BoudaryJson               `json:"boundaries"`
	MechanicalVentilations []MechanicalVentilationJson `json:"mechanical_ventilations"`
	Equipments             EquipmentsJson              `json:"equipments"`
}

type CommonJson struct {
	AcMethod string         `json:"ac_method"`
	AcConfig []ACConfigJson `json:"ac_config"`
}

type BuildingJson struct {
	Infiltration InfiltrationJson `json:"infiltration"`
}

type InfiltrationJson struct {
	Method         string  `json:"method"`
	Story          float64 `json:"story"`
	CValueEstimate string  `json:"c_value_estimate"`
	CValue         float64 `json:"c_value"`
	UAValue        float64 `json:"ua_value"`
	Structure      string  `json:"struct"`
	InsidePressure string  `json:"inside_pressure"`
}

type BoudaryJson struct {
	Id                            float64 `json:"id"`
	Name                          string  `json:"name"`
	SubName                       string  `json:"sub_name"`
	Area                          float64 `json:"area"`
	Direction                     string  `json:"direction"`
	ConnectedRoomId               float64 `json:"connected_room_id"`
	UValue                        float64 `json:"u_value"`
	EtaValue                      float64 `json:"eta_value"`
	GlassAreaRatio                float64 `json:"glass_area_ratio"`
	H_c                           float64 `json:"h_c"`
	BoundaryType                  string  `json:"boundary_type"`
	TempDifCoef                   float64 `json:"temp_dif_coef"`
	OutsideSolarAbsorption        float64 `json:"outside_solar_absorption"`
	OutsideHeatTransferResistance float64 `json:"outside_heat_transfer_resistance"`
	IncidentAngleCharacteristics  string  `json:"incident_angle_characteristics"`
	InsideHeatTransferResistance  float64 `json:"inside_heat_transfer_resistance"`

	RearSurfaceBoundaryId float64              `json:"rear_surface_boundary_id"`
	OutsideEmissivity     float64              `json:"outside_emissivity"`
	IsSolarAbsorbedInside bool                 `json:"is_solar_absorbed_inside"`
	IsSunStrikedOtside    bool                 `json:"is_sun_striked_outside"`
	IsFloor               bool                 `json:"is_floor"`
	Layers                []LayerJson          `json:"layers"`
	SolarShadingPart      SolarShadingPartJson `json:"solar_shading_part"`
}

type SolarShadingPartJson struct {
	Existence   bool    `json:"existence"`
	InputMethod string  `json:"input_method"`
	Depth       float64 `json:"depth"`
	D_h         float64 `json:"d_h"`
	D_e         float64 `json:"d_e"`
	X1          float64 `json:"x1"`
	X2          float64 `json:"x2"`
	X3          float64 `json:"x3"`
	Y1          float64 `json:"y1"`
	Y2          float64 `json:"y2"`
	Y3          float64 `json:"y3"`
	Z_x_pls     float64 `json:"z_x_pls"`
	Z_x_mns     float64 `json:"z_x_mns"`
	Z_y_pls     float64 `json:"z_y_pls"`
	Z_y_mns     float64 `json:"z_y_mns"`
}

type MechanicalVentilationJson struct {
	Id       float64 `json:"id"`
	RootType string  `json:"root_type"`
	Volume   float64 `json:"volume"`
	Root     []int   `json:"root"`
}

type LayerJson struct {
	Name              string  `json:"name"`
	ThermalCapacity   float64 `json:"thermal_capacity"`
	ThermalResistance float64 `json:"thermal_resistance"`
}

type ACConfigJson struct {
	Mode  int     `json:"mode"`
	Lower float64 `json:"lower"`
	Upper float64 `json:"upper"`
}

type RoomJson struct {
	Id          int             `json:"id"`
	Name        string          `json:"name"`
	SubName     string          `json:"sub_name"`
	FloorArea   float64         `json:"floor_area"`
	Schedule    ScheduleJson    `json:"schedule"`
	Volume      float64         `json:"volume"`
	Furniture   FurnitureJson   `json:"furniture"`
	Ventilation VentilationJson `json:"ventilation"`
}

type FurnitureJson struct {
	InputMethod      string  `json:"input_method"`
	HeatCapacity     float64 `json:"heat_capacity"`
	HeatCond         float64 `json:"heat_cond"`
	MoistureCapacity float64 `json:"moisture_capacity"`
	MoistureCond     float64 `json:"moisture_cond"`
}

type ScheduleJson struct {
	Name string `json:"name"`
}

type VentilationJson struct {
	Natural float64 `json:"natural"`
}

type EquipmentsJson struct {
	HeatingEquipments []EquipmentJson `json:"heating_equipments"`
	CoolingEquipments []EquipmentJson `json:"cooling_equipments"`
}

type EquipmentJson struct {
	Id            int                  `json:"id"`
	Name          string               `json:"name"`
	Property      EquipmentPropertJson `json:"property"`
	EquipmentType string               `json:"equipment_type"`
}

type EquipmentPropertJson struct {
	SpaceId         int     `json:"space_id"`
	BoundaryId      int     `json:"boundary_id"`
	Qmin            float64 `json:"q_min"`
	Q_max           float64 `json:"q_max"`
	V_min           float64 `json:"v_min"`
	V_max           float64 `json:"v_max"`
	Bf              float64 `json:"bf"`
	MaxCapacity     float64 `json:"max_capacity"`
	Area            float64 `json:"area"`
	ConvectionRatio float64 `json:"convection_ratio"`
}

/*
負荷計算処理の実行

	Args:
	    logger
	    house_data_path (str): 住宅計算条件JSONファイルへのパス
	    output_data_dir (str): 出力フォルダへのパス
	    weather_specify_method: 気象データの指定方法
	    weather_file_path: 気象データのファイルパス
	    region: 地域の区分
*/
func run(
	house_data_path string,
	output_data_dir string,
	weather_specify_method string,
	weather_file_path string,
	region int,
) {
	// interval currently fixed at 15 minutes
	//itv := 15
	log.SetFlags(log.Lmicroseconds)

	// ---- 事前準備 ----

	// 出力ディレクトリの作成
	if _, err := os.Stat(output_data_dir); os.IsNotExist(err) {
		os.Mkdir(output_data_dir, 0755)
	}

	_, err := os.Stat(output_data_dir)
	if os.IsNotExist(err) {
		log.Fatalf("`%s` is not a directory", output_data_dir)
	}

	// 住宅計算条件JSONファイルの読み込み
	log.Printf("住宅計算条件JSONファイルの読み込み開始")
	var rd InputJson
	if house_data_path[0:4] == "http" {
		resp, err := http.Get(house_data_path)
		if err != nil {
			log.Fatal(err)
		}
		defer resp.Body.Close()
		body, err := ioutil.ReadAll(resp.Body)
		if err != nil {
			log.Fatal(err)
		}
		json.Unmarshal(body, &rd)
	} else {
		file, err := os.Open(house_data_path)
		if err != nil {
			log.Fatal(err)
		}
		defer file.Close()
		bytes, err := ioutil.ReadAll(file)
		if err != nil {
			log.Fatal(err)
		}
		json.Unmarshal(bytes, &rd)
	}

	// 気象データの生成 => weather_for_method_file.csv
	log.Printf("気象データの生成開始")
	w := make_weather(
		weather_specify_method,
		IntervalM15,
		weather_file_path,
		Region(fmt.Sprint(region)),
	)

	log.Printf("スケジュール生成開始")
	scnames := make([]string, 0)
	floor_areas := make([]float64, 0)
	for _, rm := range rd.Rooms {
		scnames = append(scnames, rm.Schedule.Name)
		floor_areas = append(floor_areas, rm.FloorArea)
	}
	scd := get_schedule(
		Auto,
		scnames,
		floor_areas,
	)

	// ---- 計算 ----

	// 計算
	calc(&rd, w, scd, IntervalM15, 4, 365, 365, 183)

	// // ---- 計算結果ファイルの保存 ----

	// // 計算結果（瞬時値）
	// result_detail_i_path = path.join(output_data_dir, 'result_detail_i.csv')
	// logger.info('Save calculation results data (detailed version) to `{}`'.format(result_detail_i_path))
	// dd_i.to_csv(result_detail_i_path, encoding='cp932')

	// // 計算結果（平均・積算値）
	// result_detail_a_path = path.join(output_data_dir, 'result_detail_a.csv')
	// logger.info('Save calculation results data (simplified version) to `{}`'.format(result_detail_a_path))
	// dd_a.to_csv(result_detail_a_path, encoding='cp932')
}

func main() {
	var house_data string
	flag.StringVar(&house_data, "input", "", "計算を実行するJSONファイル")

	var output_data_dir string
	flag.StringVar(&output_data_dir, "o", ".", "出力フォルダ")

	var weather string
	flag.StringVar(&weather, "weather", "ees", "気象データの作成方法を指定します。")

	var weather_path string
	flag.StringVar(&weather_path, "weather_path", "", "気象データの絶対パスを指定します。weatherオプションでfileが指定された場合は必ず指定します。")

	var region int
	flag.IntVar(&region, "region", 6, "地域の区分を指定します。気象データの作成方法として建築物省エネ法を指定した場合には必ず指定します。")

	var pprpf_enable bool
	flag.BoolVar(&pprpf_enable, "pprof", false, "プロファイリングを実行し、cpu.prof ファイルに保存します。")

	// 引数を受け取る
	flag.Parse()

	if house_data == "" {
		log.Fatal("inputオプションを指定してください。")
	}

	if pprpf_enable {
		f, err := os.Create("cpu.prof")
		if err != nil {
			panic(err)
		}
		defer func() {
			if err := f.Close(); err != nil {
				panic(err)
			}
		}()
		if err := pprof.StartCPUProfile(f); err != nil {
			panic(err)
		}
		defer pprof.StopCPUProfile()
	}

	start := time.Now()

	run(
		house_data,
		output_data_dir,
		weather,
		weather_path,
		region,
	)

	elapsedTime := time.Since(start)
	log.Printf("elapsed_time: %v [sec]", elapsedTime)
}
