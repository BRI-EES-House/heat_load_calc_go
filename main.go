package main

import (
	"flag"
	"log"
	"os"
	"runtime/pprof"
	"time"

	"github.com/BRI-EES-House/heat_load_calc_go/heat_load_calc"
)

func main() {
	var house_data string
	flag.StringVar(&house_data, "i", "", "計算を実行するJSONファイル")

	var output_data_dir string
	flag.StringVar(&output_data_dir, "o", "", "出力フォルダ")

	var is_schedule_saved bool
	flag.BoolVar(&is_schedule_saved, "schedule_saved", false, "スケジュールを出力するか否かを指定します。")

	var weather string
	flag.StringVar(&weather, "weather", "ees", "気象データの作成方法を指定します。")

	var weather_path string
	flag.StringVar(&weather_path, "weather_path", "", "気象データの絶対パスを指定します。weatherオプションでfileが指定された場合は必ず指定します。")

	var region int
	flag.IntVar(&region, "region", 0, "地域の区分を指定します。気象データの作成方法として建築物省エネ法を指定した場合には必ず指定します。")

	var is_weather_saved bool
	flag.BoolVar(&is_weather_saved, "weather_saved", false, "気象データを出力するか否かを指定します。")

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

	heat_load_calc.Run(
		house_data,
		output_data_dir,
		is_schedule_saved,
		weather,
		weather_path,
		region,
		is_weather_saved,
		output_data_dir != "",
	)

	elapsedTime := time.Since(start)
	log.Printf("elapsed_time: %v [sec]", elapsedTime)
}
