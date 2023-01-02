#!/bin/bash
path="./test/"
ext_t=".csv"
ext_f="_log.txt"
n1="t_12"
fn1_t=${path}${n1}${ext_t}
fn1_f=${path}${n1}${ext_f}
echo "c,b1,b2,iter,time" > ${fn1_t}
echo -e "n\tf1\tf2\tfase\titer" > ${fn1_f}
./test.o fact 8 12 3000:1000:9 ${fn1_t} ${fn1_f}

cat ${fn1_f}
