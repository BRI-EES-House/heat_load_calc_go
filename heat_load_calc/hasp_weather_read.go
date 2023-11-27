package heat_load_calc

import (
	"bufio"
	"math"
	"os"
	"strconv"
)

/*
   HASP形式のテキストファイル（utf-8）から気温、絶対湿度、法線面直達日射量、水平面天空日射量、夜間放射量、風向、風速を読み込む
   Args:
       file_name: HASP形式のファイル名（3カラムのSI形式）

   Returns:
       気温, C
       絶対湿度, kg/kg(DA)
       法線面直達日射量, W/m2
       水平面天空日射量, W/m2
       夜間放射量, W/m2
       風向, -
       風速, m/s
*/
func hasp_read(file_name string) ([]float64, []float64, []float64, []float64, []float64, []float64, []float64, error) {
	// テキストファイルを読み込み
	f, err := os.Open(file_name)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	file, err := os.Open(file_name)
	if err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var lines []string = make([]string, 0, 365)
	for scanner.Scan() {
		line := scanner.Text()
		lines = append(lines, line)
	}

	if err := scanner.Err(); err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}

	// 外気温度
	taData, err := _get_elements(lines, 0)
	if err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}
	// 換算し端数処理
	for i := range taData {
		taData[i] = math.Round((taData[i]-500.0)*0.1*10) / 10
	}

	// 絶対湿度
	xaData, err := _get_elements(lines, 1)
	if err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}
	for i := range xaData {
		xaData[i] = math.Round(xaData[i]*0.0001*10000) / 10000
	}

	// 法線面直達日射量
	idnData, err := _get_elements(lines, 2)
	if err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}
	for i := range idnData {
		idnData[i] = math.Round(idnData[i] * (0.01 * 1.0e6 / 3600.0))
	}

	// 水平面天空日射量
	iskyData, err := _get_elements(lines, 3)
	if err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}
	for i := range iskyData {
		iskyData[i] = math.Round(iskyData[i] * (0.01 * 1.0e6 / 3600.0))
	}

	// 夜間放射量
	rnData, err := _get_elements(lines, 4)
	if err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}
	for i := range rnData {
		rnData[i] = math.Round(rnData[i] * (0.01 * 1.0e6 / 3600.0))
	}

	// 風向
	wdreData, err := _get_elements(lines, 5)
	if err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}

	// 風速
	wvData, err := _get_elements(lines, 6)
	if err != nil {
		return nil, nil, nil, nil, nil, nil, nil, err
	}
	for i := range wvData {
		wvData[i] = math.Round(wvData[i]*0.1*10) / 10
	}

	return taData, xaData, idnData, iskyData, rnData, wdreData, wvData, nil
}

func _get_elements(lines []string, row_number int) ([]float64, error) {
	element_data := make([]float64, 0, 24*len(lines)/7)
	for i := row_number; i < len(lines); i += 7 {
		row := lines[i]
		for j := 0; j < 24*3; j += 3 {
			num, err := strconv.Atoi(row[j : j+3])
			if err != nil {
				return nil, err
			}
			element_data = append(element_data, float64(num))
		}
	}

	return element_data, nil
}
