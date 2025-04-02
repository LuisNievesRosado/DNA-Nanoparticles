#!/bin/bash

for i in $(seq 12 1.0 50.0)
do
	mkdir $i
	mv DX$i.data  $i/in.data
done

sbatch runAA.sh
