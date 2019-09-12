#!/usr/bin/env python3

# Note on talapas you must do
#   module load racs-eb  FFmpeg/3.4.1-intel-2017b
# first

import sys
import pyslim
import numpy as np
import glob

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import animation as animation

treefiles = sys.argv[1:]

if len(sys.argv) == 1:
    raise ValueError("Usage:  anim_locations.py [treefile [treefile]]")

for treefile in treefiles:
    outfile = ".".join(treefile.rsplit(".")[:-1] + ["mp4"])
    print(treefile + " -> " + outfile)

    ts = pyslim.load(treefile)
    births = np.array([ts.node(i.nodes[0]).time for i in ts.individuals()])
    birth_list = sorted(list(set(births)), reverse=True)
    ages = np.array([x.age for x in pyslim.extract_individual_metadata(ts.tables)])
    if all(ages == -1):
        # WF model
        ages[:] = 0

    locations = np.array([i.location for i in ts.individuals()])
    fig_dims = (max(locations[:,0]), max(locations[:,1]))

    # hack to correct for https://github.com/MesserLab/SLiM/issues/23
    zeros = np.logical_and(locations[:,0] < 1e-16, locations[:,1] < 1e-16)
    locations[zeros,0] = np.random.uniform(0, fig_dims[0], sum(zeros))
    locations[zeros,1] = np.random.uniform(0, fig_dims[1], sum(zeros))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_facecolor('black')

    # return locations of individuals alive at time n-in-the-past
    def get_locs(n):
        alive = np.logical_and(births >= n, births - ages <= n)
        return (locations[alive, :2], ages[alive])

    def update_points(n):
        locs, age_n = get_locs(n)
        pts.set_offsets(locs)
        pts.set_color(plt.get_cmap('autumn_r')(age_n))
        return pts

    locs, age_n = get_locs(0)
    pts = ax.scatter(locs[:,0], locs[:,1],
                     s = 5, marker = '.',
                     c = plt.get_cmap('autumn_r')(age_n),
                     linewidths = 0)

    fig.set_size_inches([3,3])
    plt.tight_layout()

    ani = animation.FuncAnimation(fig, update_points, 
                                  frames=range(int(max(birth_list)), -1, -1), 
                                  interval=150)
    writer = animation.writers['ffmpeg'](fps=10)
    ani.save(outfile, writer=writer, dpi=300)

print("... done!")

