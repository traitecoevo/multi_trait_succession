#!/bin/bash

f='output/data/[60,1]/2trait/base/'
t=2
src/evolve -f $f -t $t -L 1 -s 1
src/evolve -f $f -t $t -L 1 -s 5
src/evolve -f $f -t $t -L 1 -s 10
src/evolve -f $f -t $t -L 1 -s 200

f='output/data/[60,1]/hsp/[-0.2]/base'
t=1
for ((i = 0; i <= 19; i=i+1)); do
    echo $i
    src/evolve -f $f -t $t -L 1 -s $i
done

f='output/data/[60,1]/hsp/[0.325]/base'
t=1
for ((i = 0; i <= 19; i=i+1)); do
    echo $i
    src/evolve -f $f -t $t -L 1 -s $i
done

f='output/data/[60,1]/hsp/[0.625]/base'
t=1
for ((i = 0; i <= 19; i=i+1)); do
    echo $i
    src/evolve -f $f -t $t -L 1 -s $i
done

f='output/data/[60,1]/lcc/[0.8175]/base'
t=0

for ((i = 0; i <= 19; i=i+1)); do
    echo $i
    src/evolve -f $f -t $t -L 1 -s $i
done

f='output/data/[60,1]/lcc/[1.3725]/base'
t=0

for ((i = 0; i <= 19; i=i+1)); do
    echo $i
    src/evolve -f $f -t $t -L 1 -s $i
done

f='output/data/[60,1]/lcc/[1.5575]/base'
t=0

for ((i = 0; i <= 19; i=i+1)); do
    echo $i
    src/evolve -f $f -t $t -L 1 -s $i
done

f='output/data/[60,1]/hsp/[0.54]/[c_r1,0.4]'
t=1
for ((i = 0; i <= 19; i=i+1)); do
    echo $i
    src/evolve -f $f -t $t -L 1 -s $i
done

f='output/data/[60,1]/hsp/[0.54]/[c_r1,0.8]'
t=1
for ((i = 0; i <= 19; i=i+1)); do
    echo $i
    src/evolve -f $f -t $t -L 1 -s $i
done
