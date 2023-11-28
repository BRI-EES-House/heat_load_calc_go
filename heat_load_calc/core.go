package heat_load_calc

import (
	"log"
)

/*
coreメインプログラム

	Args:
	    rd: 住宅計算条件
	    w: 外界気象条件
	    scd: スケジュール
	    itv: 時間間隔
	    n_step_hourly: 計算間隔（1時間を何分割するかどうか）（デフォルトは4（15分間隔））
	    n_d_main: 本計算を行う日数（デフォルトは365日（1年間））, d
	    n_d_run_up: 助走計算を行う日数（デフォルトは365日（1年間））, d
	    n_d_run_up_build: 助走計算のうち建物全体を解く日数（デフォルトは183日（およそ半年））, d

	Returns:
	    以下のタプル
	        (1) 計算結果（詳細版）をいれたDataFrame
	        (2) 計算結果（簡易版）をいれたDataFrame

	Notes:
	    「助走計算のうち建物全体を解く日数」は「助走計算を行う日数」で指定した値以下でないといけない。
*/
func calc(
	rd *InputJson,
	w *Weather,
	scd *Schedule,
	itv Interval,
	n_step_hourly int,
	n_d_main int,
	n_d_run_up int,
	n_d_run_up_build int,
) (string, string, *Boundaries) {
	log.Printf("計算開始")

	// 本計算のステップ数
	// 助走計算のステップ数
	// 助走計算のうち建物全体を解くステップ数
	n_step_main, n_step_run_up, n_step_run_up_build := get_n_step(n_step_hourly, n_d_main, n_d_run_up, n_d_run_up_build)

	// 時間間隔, s
	//deltaT := itv.get_delta_t()

	// json, csv ファイルからパラメータをロードする。
	// （ループ計算する必要の無い）事前計算を行い, クラス PreCalcParameters, PreCalcParametersGround に必要な変数を格納する。
	sqc := NewSequence(itv, rd, w, scd)

	pp := sqc.pre_calc_parameters

	gc_n := initialize_ground_conditions(sqc.bs.n_ground)

	log.Println("助走計算（土壌のみ）")
	N := scd.ac_demand_is_ns.Len()

	nn := N - n_step_run_up_build
	for n := -n_step_run_up; n < -n_step_run_up_build; n++ {
		gc_n = sqc.run_tick_ground(gc_n, n, nn)
		nn++
	}

	result := NewRecorder(n_step_main, sqc.rms.id_rm_is, sqc.bs.id_bs_js, itv)

	result.pre_recording(sqc.weather, sqc.scd, sqc.bs, pp.q_sol_frt_is_ns, pp.q_s_sol_js_ns)

	// 建物を計算するにあたって初期値を与える
	c_n := initialize_conditions(sqc.rms.n_rm, sqc.bs.n_b)

	// 地盤計算の結果（項別公比法の指数項mの吸熱応答の項別成分・表面熱流）を建物の計算に引き継ぐ
	update_conditions_by_ground_conditions(sqc.bs.is_ground_js, c_n, gc_n)

	log.Printf("助走計算（建物全体） (%d steps)\n", n_step_run_up_build)

	// NOTE: XX_plusの配列はインデックスが負の場合に適切な位置を参照しているのか良く確認する。
	N_plus := N + 1
	nn = N - n_step_run_up_build
	nn_plus := N_plus - n_step_run_up_build
	for n := -n_step_run_up_build; n < 0; n++ {
		c_n = sqc.run_tick(n, nn, nn_plus, c_n, result)
		nn++
	}

	log.Printf("本計算 (%d steps)\n", n_step_main)

	// TODO: recorder に1/1 0:00の瞬時状態値を書き込む
	m := 1
	for n := 0; n < n_step_main; n++ {
		c_n = sqc.run_tick(n, n, n, c_n, result)
		if n == int(float64(n_step_main)/12*float64(m)) {
			log.Printf("%d / 12 calculated.", m)
			m++
		}
	}
	log.Print("12 / 12 calculated.")

	result.post_recording(sqc.rms, sqc.bs, pp.f_mrt_is_js, sqc.es)

	log.Println("ログ作成")

	// dd: data detail, 15分間隔のすべてのパラメータ pd.DataFrame
	dd_a, dd_i := result.export_pd()

	return dd_a, dd_i, sqc.bs
}
