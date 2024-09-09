#!/bin/bash

for i in $(seq 12 1.0 50.0)
do
	mkdir AA/$i
	mv AA/DX$i.data  AA/$i/in.data
done
cp runAA.sh AA
cp hamilAA.in AA
cd AA
sbatch runAA.sh
cd ../
	

for i in $(seq 12 1.0 50.0)
do
	mkdir AB/$i
	mv AB/DX$i.data  AB/$i/in.data
done
cp runAB.sh AB
cp hamilAB.in AB
cd AB
sbatch runAB.sh
cd ../
	

for i in $(seq 12 1.0 50.0)
do
	mkdir BB/$i
	mv BB/DX$i.data  BB/$i/in.data
done
cp runBB.sh BB
cp hamilBB.in BB
cd BB
sbatch runBB.sh
cd ../
	
