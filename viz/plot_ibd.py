#!/usr/bin/env python3

import sys
import pyslim, msprime
import numpy as np
import glob

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

treefiles = sys.argv[1:]

if len(sys.argv) == 1:
    raise ValueError("Usage:  plot_ibd.py [treefile [treefile]]")

for treefile in treefiles:
    outfile = ".".join(treefile.rsplit(".")[:-1] + ["ibd", "png"])
    print(treefile + " -> " + outfile)

    ts = pyslim.SlimTreeSequence(pyslim.load(treefile).simplify())
    rts = ts.recapitate(recombination_rate = 1e-8, Ne=ts.num_samples)
    samples = np.random.choice(rts.num_individuals, 300, replace=False)
    sample_genomes = [list(rts.individual(ind).nodes) for ind in samples]

    locs = np.array([rts.individual(i).location[:2] for i in samples])
    dists = np.sqrt((locs[:,0][:,np.newaxis] - locs[:,0]) ** 2
                    + (locs[:,1][:,np.newaxis] - locs[:,1]) ** 2)

    the_genomes = [x for y in sample_genomes for x in y]
    remapped_genomes = [[the_genomes.index(u) for u in x] for x in sample_genomes]
    sub_ts = rts.simplify(the_genomes)
    bc = msprime.BranchLengthStatCalculator(sub_ts)
    divs = bc.divergence_matrix(remapped_genomes, windows=[0.0, sub_ts.sequence_length])
    ut = np.triu_indices(dists.shape[0])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    pts = ax.scatter(dists[ut], divs[0][ut]/1000, 
                     alpha  = 0.5, s = 5, marker = '.', linewidths = 0)
    plt.xlabel("geographic distance")
    plt.ylabel("divergence")
    fig.set_size_inches([4,4])
    plt.savefig(outfile, dpi=288)

print("... done!")

