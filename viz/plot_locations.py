#!/usr/bin/env python3

import sys
import pyslim
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

if len(sys.argv) == 1:
    raise ValueError("Usage:  plot_locations.py [treefile [treefile]]")

treefiles = sys.argv[1:]

for treefile in treefiles:
    print(treefile)
    outfile = ".".join(treefile.rsplit(".")[:-1] + ["locs", "png"])
    ts = pyslim.load(treefile).simplify()
    locs = np.array([i.location[:2] for i in ts.individuals()])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    pts = ax.scatter(locs[:,0], locs[:,1], s=5, c="black")
    fig.set_size_inches([3,3])
    plt.savefig(outfile, dpi=288)
    print(".")

print("... done!")

