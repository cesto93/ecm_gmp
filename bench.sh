#!/bin/bash

path="./test/"
ext_csv=".csv"
ext_log="_time.csv"

start_iter=20
end_iter=25
sub_iter=1

declare -A log csv bin
declare -a b_values

b_values[20]=8000:2000:${sub_iter}
b_values[21]=15000:5000:${sub_iter}
b_values[22]=20000:5000:${sub_iter}
b_values[23]=30000:5000:${sub_iter}
b_values[24]=45000:5000:${sub_iter}
b_values[25]=50000:5000:${sub_iter}
b_values[26]=10000:5000:${sub_iter}
b_values[27]=150000:5000:${sub_iter}
b_values[28]=20000:5000:${sub_iter}
b_values[29]=25000:5000:${sub_iter}

n1="bench_${start_iter}_${end_iter}"
log[c]=${path}${n1}${ext_csv}
csv[c]=${path}${n1}${ext_log}
bin[c]="./test.o"

log[wasm]=${path}${n1}_wasm${ext_csv}
csv[wasm]=${path}${n1}_wasm${ext_log}
bin[wasm]="./test.wasm"

log[mm]=${path}${n1}_mm${ext_csv}
csv[mm]=${path}${n1}_mm${ext_log}
bin[mm]="./mmtest.o"

run_bench() {
	label=${1}
	rep=10
	local i

	echo "c,b1,b2,iter,time" > ${csv["${label}"]}
	echo -e "n\tf1\tf2\tfase\titer" > ${log["${label}"]}

	for ((i = start_iter; i <= end_iter; i++)); do
		${bin["${label}"]} fact ${rep} ${i} ${b_values[$i]} ${csv["${label}"]} ${log["${label}"]}
	done
}

run_bench c
#run_bench wasm
run_bench mm
