import pyslim
import numpy as np
import matplotlib
import glob
matplotlib.use('Agg')
from matplotlib import pyplot as plt

treefiles = glob.glob("*.trees")

for treefile in treefiles:
    print(treefile)
    outbase = treefile.split(".")[0]
    ts = pyslim.load(treefile)
    ages = np.array([ts.node(i.nodes[0]).time for i in ts.individuals()])
    locations = np.array([i.location for i in ts.individuals()])
    fig_dims = (max(locations[:,0]), max(locations[:,1]))
    for k, age in enumerate(sorted(list(set(ages)))):
        outfile = outbase + ".{:03d}.png".format(k)
        locs = locations[ages == age, :2]
        fig = plt.figure(figsize=(3, 3))
        plt.scatter(locs[:,0], locs[:,1], marker='.')
        plt.xlim(0, fig_dims[0])
        plt.ylim(0, fig_dims[1])
        plt.savefig(outfile, dpi=288)
        print(".")

    print("... done!")

