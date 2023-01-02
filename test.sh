#!/bin/bash
path="./test/"
ext_t=".csv"
ext_f="_log.txt"
n1="t_16:25"
fn1_t=$path$n1$ext_t
fn1_f=$path$n1$ext_f
echo "c,b1,b2,iter,time" > $fn1_t
echo -e "n\tf1\tf2\tfase\titer" > $fn1_f
./test.o fact 8 12 3000:1000:9 $fn1_t $fn1_f
#./test.o fact 8 16 3000:1000:9 $fn1_t $fn1_f
#./test.o fact 8 17 4000:1000:9 $fn1_t $fn1_f
#./test.o fact 8 18 5000:1000:8 $fn1_t $fn1_f
#./test.o fact 8 19 6000:1000:8 $fn1_t $fn1_f
#./test.o fact 8 20 7000:1000:7 $fn1_t $fn1_f
#./test.o fact 8 21 15000:5000:8 $fn1_t $fn1_f
#./test.o fact 8 22 20000:5000:8 $fn1_t $fn1_f
#./test.o fact 8 23 30000:5000:7 $fn1_t $fn1_f
#./test.o fact 8 24 40000:5000:5 $fn1_t $fn1_f
#./test.o fact 8 25 45000:5000:4 $fn1_t $fn1_f
#./test.o fact 8 26 100000:25000:5 $fn1_t $fn1_f
#./test.o fact 8 27 125000:25000:5 $fn1_t $fn1_f
#./test.o fact 8 28 150000:25000:5 $fn1_t $fn1_f
#./test.o fact 8 29 175000:25000:5 $fn1_t $fn1_f
#./test.o fact 8 30 200000:25000:5 $fn1_t $fn1_f
