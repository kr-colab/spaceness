import numpy as np, multiprocessing as mp, subprocess as sp
import os, sys, re, pyslim, msprime, sklearn, sys
from SLiMomatic import *
import time

wd="/projects/kernlab/cbattey2/spaceness/"
os.chdir(wd)

#initialize the data generator
dg=SLiMomatic(
    SLiMExecutableFilePath="slim",
    SLiMRecipe="flat_map_k10.slim",
    numReps=100,
    Ne=9e3,
    recombinationRate=1e-9,
    params_to_sample=["sigma_i"],
    priorLows=[0.2],
    priorHighs=[1])

#run simulations in parallel
dg.simulateAndProduceTrees(direc=os.path.join(wd,"sims/trees/k10"),nProc=28)
