module github.com/BRI-EES-House/heat_load_calc_go

go 1.18

require (
	github.com/gocarina/gocsv v0.0.0-20230406101422-6445c2b15027
	gonum.org/v1/gonum v0.9.3
	gonum.org/v1/netlib v0.0.0-20220323200511-14de99971b2d
)

require golang.org/x/exp v0.0.0-20210220032938-85be41e4509f // indirect

replace github.com/BRI-EES-House/heat_load_calc_go/heat_load_calc => ./heat_load_calc
