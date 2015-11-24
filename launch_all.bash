#!/bin/bash

## Runs simulations for HSP, holding LCC constant, 40 processes about 4 days each
d=60
p=1
t=1
for v in -2.000 -1.925 -1.850 -1.775 -1.700 -1.625 -1.550 -1.475 -1.400 -1.325 -1.250 -1.175 -1.100 -1.025 -0.950 -0.875 -0.800 -0.725 -0.650 -0.575 -0.500 -0.425 -0.350 -0.275 -0.200 -0.125 -0.050  0.025  0.100  0.175  0.250  0.325  0.400  0.475  0.550  0.625  0.700  0.775  0.850  0.925  1.000;   # 60  calls  (R command, seq(-2, 1, length.out=41))  # test for v in -0.8 0.3 0.6;
	do
		nohup src/evolve -p $p -d $d -t $t -v $v > output/out/trait-$t.[$d,$p]-$v.txt 2>&1 &
	done


## Runs simulations for LCC, holding HSP constant, 40 processes about 4 days each

d=60
p=1
t=0
for v in -0.20000 -0.15375 -0.10750 -0.06125 -0.01500  0.03125  0.07750  0.12375  0.17000  0.21625  0.26250  0.30875  0.35500  0.40125  0.44750  0.49375  0.54000  0.58625  0.63250 0.67875  0.72500  0.77125  0.81750  0.86375  0.91000  0.95625  1.00250  1.04875  1.09500  1.14125  1.18750  1.23375  1.28000  1.32625  1.37250  1.41875  1.46500  1.51125 1.55750  1.60375  1.65000; do   # 40  calls  R command, seq(-0.2, 1.65, length.out=41)  # test for v in -1 0 1.4; do
	nohup src/evolve -p $p -d $d -t $t -v $v > output/out/trait-$t.[$d,$p]-$v.txt 2>&1 &
done


## Begin 2D simulations, 5x11 = 55 processes, ranging from few days to > 6 weeks each
t=2
# test for d in 15 30 60 120 240; do for p in 1; do
for d in 7.5  10.6  15  21.2  30  42.4  60  84.8 120 169.7 240; do
 for p in 0.6 0.8 1 1.2 1.4; do
	nohup src/evolve -V 1 -p $p -d $d -t $t > output/out/trait-$t.[$d,$p]-base.txt 2>&1 &
 done
done


## Test sensitivity to different parameter values - 32 processes, about 1-2 weeks each
d=60
p=1
t=2
for a in 0.9 1.1; do
	for x in {0..15}; do
		nohup src/evolve -V 1 -p $p -d $d -x $x -a $a -t $t > output/out/trait-$t.[$d,$p]-elas-x$x-a$a.txt 2>&1 &
	done
done


## Test sensitivity to assumption about maximal reproductive allocation, 5 processes about 1-2 weeks each
d=60
p=1
t=2
x=16
for a in 0.2 0.4 0.6 0.8; do
	nohup src/evolve -V 1 -p $p -d $d -x $x -a $a -t $t > output/out/trait-$t.[$d,$p]-elas-x$x-a$a.txt 2>&1 &
done

## Test sensitivity to assumption about steepness of reproductive allocation function, 6 process about 1-2 weeks each
d=60
p=1
t=2
x=17
for a in 0.1 0.2 0.4 0.6 0.8 1.0; do
	nohup src/evolve -V 1 -p $p -d $d -x $x -a $a -t $t > output/out/trait-$t.[$d,$p]-elas-x$x-a$a.txt 2>&1 &
done

## Test sensitivity to assumption about maximal reproductive allocation
## Runs simulations for HMAT, holding LMA constant
d=60
p=1
t=1
v= 0.54000
x=16 
for a in 0.2 0.4 0.6 0.8; do 
	nohup src/evolve -V 1 -p $p -d $d -x $x -a $a -t $t -v $v > output/out2/trait-$t.[$d,$p]-$v-elas-x$x-a$a.txt 2>&1 &
done
