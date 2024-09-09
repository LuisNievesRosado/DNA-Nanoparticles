#!/bin/bash 

cp makehistH.py 34548
cp makeprob.py 34548
cp PMF2pot.py 34548
cp wham 34548
cd 34548
python3 makehistH.py 
python3 makeprob.py 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaAA.wham AA.out 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaAB.wham AB.out 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaBB.wham BB.out 
python3 PMF2pot.py 
cd .. 

cp makehistH.py 34802
cp makeprob.py 34802
cp PMF2pot.py 34802
cp wham 34802
cd 34802
python3 makehistH.py 
python3 makeprob.py 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaAA.wham AA.out 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaAB.wham AB.out 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaBB.wham BB.out 
python3 PMF2pot.py 
cd .. 

cp makehistH.py 34803
cp makeprob.py 34803
cp PMF2pot.py 34803
cp wham 34803
cd 34803
python3 makehistH.py 
python3 makeprob.py 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaAA.wham AA.out 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaAB.wham AB.out 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaBB.wham BB.out 
python3 PMF2pot.py 
cd .. 

cp makehistH.py 34813
cp makeprob.py 34813
cp PMF2pot.py 34813
cp wham 34813
cd 34813
python3 makehistH.py 
python3 makeprob.py 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaAA.wham AA.out 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaAB.wham AB.out 
./wham 12.0 50.0 39 0.0000001 1.0 0.0 MetaBB.wham BB.out 
python3 PMF2pot.py 
cd .. 

