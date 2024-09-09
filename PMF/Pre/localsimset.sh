#!/bin/bash 

mkdir 171113

python3 write_run.py 171113
mv hamilAA.in 171113
mv hamilAB.in 171113
mv hamilBB.in 171113
mv runAA.sh 171113
mv runAB.sh 171113
mv runBB.sh 171113
cp simset.sh 171113

mkdir 171113/AA 
mkdir 171113/AB 
mkdir 171113/BB 

python3 makedata_nb.py AA 28 51 14 3
mv PNP.data 171113/AA 
cp base.in 171113/AA 

python3 makedata_nb.py AB 28 51 14 3
mv PNP.data 171113/AB 
cp base.in 171113/AB 

python3 makedata_nb.py BB 28 51 14 3
mv PNP.data 171113/BB 
cp base.in 171113/BB 

mkdir 171207

python3 write_run.py 171207
mv hamilAA.in 171207
mv hamilAB.in 171207
mv hamilBB.in 171207
mv runAA.sh 171207
mv runAB.sh 171207
mv runBB.sh 171207
cp simset.sh 171207

mkdir 171207/AA 
mkdir 171207/AB 
mkdir 171207/BB 

python3 makedata_nb.py AA 29 50 13 2
mv PNP.data 171207/AA 
cp base.in 171207/AA 

python3 makedata_nb.py AB 29 50 13 2
mv PNP.data 171207/AB 
cp base.in 171207/AB 

python3 makedata_nb.py BB 29 50 13 2
mv PNP.data 171207/BB 
cp base.in 171207/BB 

mkdir 164301

python3 write_run.py 164301
mv hamilAA.in 164301
mv hamilAB.in 164301
mv hamilBB.in 164301
mv runAA.sh 164301
mv runAB.sh 164301
mv runBB.sh 164301
cp simset.sh 164301

mkdir 164301/AA 
mkdir 164301/AB 
mkdir 164301/BB 

python3 makedata_nb.py AA 29 51 12 2
mv PNP.data 164301/AA 
cp base.in 164301/AA 

python3 makedata_nb.py AB 29 51 12 2
mv PNP.data 164301/AB 
cp base.in 164301/AB 

python3 makedata_nb.py BB 29 51 12 2
mv PNP.data 164301/BB 
cp base.in 164301/BB 

