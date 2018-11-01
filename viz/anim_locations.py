#!/usr/bin/env python3

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

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # return locations of individuals alive at time n-in-the-past
    def get_locs(n):
        alive = np.logical_and(births >= n, births - ages <= n)
        return (locations[alive, :2], ages[alive])

    def update_points(n):
        locs, age_n = get_locs(n)
        pts.set_offsets(locs)
        pts.set_color(age_n)
        return pts

    locs, age_n = get_locs(1)
    pts = ax.scatter(locs[:,0], locs[:,1], s=5, c=age_n)

    fig.set_size_inches([3,3])
    plt.tight_layout()

    ani = animation.FuncAnimation(fig, update_points, frames=range(max(birth_list), 2, -1), 
                                  interval=200)
    writer = animation.writers['ffmpeg'](fps=10)
    ani.save(outfile, writer=writer, dpi=100)


    #     outfile = outbase + ".{:03d}.png".format(k)
    #     locs = locations[births == birth, :2]
    #     fig = plt.figure(figsize=(3, 3))
    #     plt.scatter(locs[:,0], locs[:,1], marker='.')
    #     plt.xlim(0, fig_dims[0])
    #     plt.ylim(0, fig_dims[1])
    #     plt.savefig(outfile, dpi=288)
    #     print(".")

    print("... done!")

