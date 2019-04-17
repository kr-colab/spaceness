#run msprime simulations with Ne equal to census sizes from spatial SLiM simulations
import pandas, numpy as np, msprime, os
os.chdir("/Users/cj/spaceness/sumstats/")
popsizes=np.loadtxt("ss_spatial_random_W50.txt.popsizes")

for i in range(len(popsizes)):
    trees=msprime.simulate(sample_size=120,
                           Ne=popsizes[i,1],
                           recombination_rate=1e-9,
                           mutation_rate=0,
                           length=1e8)
    outfile=("/Users/cj/spaceness/sims/msp/sigma_"+
             str(popsizes[i,0])+
             "_.trees_"+
             str(int(np.random.uniform(1e7,1e8-1,1))))
    trees.dump(outfile)
