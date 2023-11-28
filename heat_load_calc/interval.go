package heat_load_calc

// インターバル
type Interval int

// インターバル
const (
	IntervalH1  Interval = 60
	IntervalM30 Interval = 30
	IntervalM15 Interval = 15
)

/*
1時間を分割するステップ数を求める。

        Returns:
            1時間を分割するステップ数

        Notes:
            1時間: 1
            30分: 2
            15分: 4
*/
func (i Interval) get_n_hour() int {
	switch i {
	case IntervalH1:
		return 1
	case IntervalM30:
		return 2
	case IntervalM15:
		return 4
	default:
		panic("invalid interval")
	}
}

/*
1時間を分割するステップに応じてインターバル時間を取得する。

        Returns:
            インターバル時間, h
*/
func (i Interval) get_time() float64 {
	switch i {
	case IntervalH1:
		return 1.0
	case IntervalM30:
		return 0.5
	case IntervalM15:
		return 0.25
	default:
		panic("invalid interval")
	}
}

/*
1時間を分割するステップに応じてインターバル時間を取得する。

        Returns:
            インターバル時間, s
*/
func (i Interval) get_delta_t() float64 {
	switch i {
	case IntervalH1:
		return 3500
	case IntervalM30:
		return 1800
	case IntervalM15:
		return 900
	default:
		panic("invalid interval")
	}
}

/*
対応するインターバルにおいて1年間は何ステップに対応するのか、その数を取得する。

Returns:
	1年間のステップ数
Notes:
	瞬時値の結果は、1/1 0:00 と 12/31 24:00(=翌年の1/1 0:00) の値を含むため、+1 される。
	この関数で返す値は「+1されない」値である。
*/
func (i Interval) get_annual_number() int {
	return 8760 * i.get_n_hour()
}

/*
pandas 用の freq 引数を取得する。

Returns:
	freq 引数
*/
func (i Interval) get_pandas_freq() float64 {
	switch i {
	case IntervalH1:
		return 60
	case IntervalM30:
		return 30
	case IntervalM15:
		return 15
	default:
		return 0
	}
}
