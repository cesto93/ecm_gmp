#!/bin/bash

path="./test/"
n1="bench_20_25"
ext_csv=".csv"
ext_log="_time.csv"

declare -A log csv bin

log[c]=${path}${n1}${ext_csv}
csv[c]=${path}${n1}${ext_log}
bin[c]="./test.o"

log[wasm]=${path}${n1}_wasm${ext_csv}
csv[wasm]=${path}${n1}_wasm${ext_log}
bin[wasm]="./test.wasm"

run_bench() {
	label=${1}
	rep=10

	echo "c,b1,b2,iter,time" > ${csv["${label}"]}
	echo -e "n\tf1\tf2\tfase\titer" > ${log["${label}"]}

	${bin["${label}"]} fact ${rep} 20 8000:2000:4 ${csv["${label}"]} ${log["${label}"]}
	${bin["${label}"]} fact ${rep} 21 15000:5000:4 ${csv["${label}"]} ${log["${label}"]}
	${bin["${label}"]} fact ${rep} 22 20000:5000:4 ${csv["${label}"]} ${log["${label}"]}
	${bin["${label}"]} fact ${rep} 23 25000:5000:4 ${csv["${label}"]} ${log["${label}"]}
	${bin["${label}"]} fact ${rep} 24 30000:5000:4 ${csv["${label}"]} ${log["${label}"]}
	${bin["${label}"]} fact ${rep} 25 35000:5000:4 ${csv["${label}"]} ${log["${label}"]}
}

run_bench c
run_bench wasm


#${bin} fact 8 21 15000:5000:8 ${fn1_t} ${fn1_f}
#${bin} fact 8 22 20000:5000:8 ${fn1_t} ${fn1_f}
#${bin} fact 8 23 30000:5000:7 ${fn1_t} ${fn1_f}
#${bin} fact 8 24 40000:5000:5 ${fn1_t} ${fn1_f}
#${bin} fact 8 25 45000:5000:4 ${fn1_t} ${fn1_f}
#${bin} fact 8 26 100000:25000:5 ${fn1_t} ${fn1_f}
#${bin} fact 8 27 125000:25000:5 ${fn1_t} ${fn1_f}
#${bin} fact 8 28 150000:25000:5 ${fn1_t} ${fn1_f}
#${bin} fact 8 29 175000:25000:5 ${fn1_t} ${fn1_f}
#${bin} fact 8 30 200000:25000:5 ${fn1_t} ${fn1_f}
