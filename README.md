This repository is a simple port of heat_load_calc (12c6ca59db9799f9f5aea14faeb48b42625a53be) to Go. However, some processes such as test code and saving results have been omitted.

For the sake of simplicity, the naming conventions of the Go language are ignored and are faithful to the original; since lower-case beginnings are private in the Go language, there is no package separation and everything is included in the main package.

A large number of power calculations and trigonometric functions are performed in sun position calculations, etc., but the speed reduction caused by these functions cannot be ignored, so they are multiplied or read back to literals as appropriate.

Gonum is used for numerical calculation, but it does not support broadcast unlike numpy, so there are many places to expand with For statement. Gonum is designed to save memory allocation time by processing while overwriting variable areas, but this has not been utilized at this time.

Gonum seems to be very conscious of speeding up the process by using BLAS, but I cannot confirm the operation because the relevant package causes build errors in the development environment. Therefore, various operations are performed only in Go native code.

In order to speed up the process, we believe the following should be kept in mind
1. Always allocate memory for matrices generated as local variables in run_tick, and those that can be allocated in Condition should be allocated in Condition.
2. if a table lookup using map is performed within a function, but the number of executions is large, it should be a switch statement. (Global variables should be avoided. Also, Go does not have static.)
3. constants should be defined as constants.

# Quick Start

```
$ go run main.go -input example/data_example1.json --region 6
```