#!/bin/bash
path="./test/"
ext_t=".csv"
ext_f="_log.txt"
n1="t_11"
fn1_t=${path}${n1}${ext_t}
fn1_f=${path}${n1}${ext_f}

echo "c,b1,b2,iter,time" > ${fn1_t}
echo -e "n\tf1\tf2\tfase\titer" > ${fn1_f}
./test.o fact 10 11 3000:1000:1 ${fn1_t} ${fn1_f}

cat ${fn1_f}
rm ${fn1_f} ${fn1_t} 
