#smcpp runs for simulated nwf model output
import pyslim, msprime, numpy as np, os
os.chdir("/Users/cj/spaceness/demography/smcpp")

ts=pyslim.load("data/sigma_0.3095957995572229_.trees1500000.trees")
