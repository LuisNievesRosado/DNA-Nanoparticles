import os
import matplotlib.pyplot as plt
import numpy as np
import sys as s
import statistics
import math as m
import pickle
import networkx as nx
import netrd as nr

INFILE1 = s.argv[1] + '.pickle'
INFILE2 = 'gyrref.pickle'

with open(INFILE1, 'rb') as f1:
    DG1 = pickle.load(f1)

with open(INFILE2, 'rb') as f2:
    DG2 = pickle.load(f2)

G1L = DG1[1]
G2L = DG2[1]

dist = nr.distance.PortraitDivergence() # Metric works very well, also physically intuitive
D = dist.dist(G1L,G2L)





print(D)