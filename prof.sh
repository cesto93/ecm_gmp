#!/bin/bash
par=$1
if [ "$par" = "static" ]; then
make prof_s
else
make prof
fi
./prof.o 1000000000000000000000001970000000000000000000000871 50000
gprof ./prof.o >prof.txt

