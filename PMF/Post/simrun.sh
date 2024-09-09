#!/bin/bash

for j in 6750
do
	mkdir $j
	for i in 2.0 2.5 3.0 3.5 4.0 4.5 
	do 
	 	mkdir $j/$i
		python3 writerun.py $j $((10000-$j)) $i 1.0
		mv run.nbs $j/$i
		cp SS1.in ABP.pot AAP.pot BBP.pot $j/$i
		cd $j/$i
		sbatch run.nbs
		cd ../..
	done
done

