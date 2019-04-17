#spaceness lineage-through-time / coalescence time plots
from slimtools import *
from matplotlib import pyplot as plt
import matplotlib as mpl

#get n nodes by time bin across different sigma values
os.chdir("/Users/cj/spaceness/sims/slimout/spatial/W50_run3/")
trees=[x for x in os.listdir() if not x.startswith(".")]
nodes_per_timebin=[]
for t in trees:
    sigma=float(re.sub("sigma_|_.trees*|.trees","",t))
    ts=sample_treeseq(infile=t,
                     outfile="",
                     nSamples=100,
                     sampling="random",
                     recapitate=True,
                     recombination_rate=1e-8,
                     write_to_file=False,
                     seed=12345)
    nodetimes=[x.time for x in ts.nodes()]
    n_nodes=0
    nodes_per_timebin.append([0,0,sigma])
    for m in np.arange(1,1.6e6-1,1e3):
        n_nodes=n_nodes+len([x for x in nodetimes if  x >= m and x < m+1e3])
        nodes_per_timebin.append([m+1e3,n_nodes,sigma])
nodes_per_timebin=np.array(nodes_per_timebin)

#line plot version (but how to get the colorbar??)
fig = plt.figure(figsize=(2.5,2),dpi=600)
sigmas=np.unique(nodes_per_timebin[:,2])
colorscale=np.linspace(0,1,len(sigmas))
for i in range(len(sigmas)):
    plt.plot((np.log10(np.array([x for x in nodes_per_timebin if x[2]==sigmas[i]])[:,0]/4)/1e4),
             np.array([x for x in nodes_per_timebin if x[2]==sigmas[i]])[:,1]/np.max(np.array([x for x in nodes_per_timebin if x[2]==sigmas[i]])[:,1]),
             c=plt.cm.RdYlBu(colorscale[i]),clip_on=False,linewidth=.75)
#plt.colorbar(shrink=0.9).set_label("sigma", labelpad=-13, y=1.08, rotation=0)
plt.xlabel('Generations/N',fontsize=7)
plt.ylabel('Proportion Total Coalescences',fontsize=7)
fig.savefig("/Users/cj/Desktop/test2.pdf",bbox_inches='tight')

#scatter plot version
plt.scatter(x=(nodes_per_timebin[:,0]/4)/1e4,y=nodes_per_timebin[:,1],c=nodes_per_timebin[:,2],
            s=10,cmap=plt.cm.get_cmap("RdYlBu"),
            clip_on=False,alpha=0.5)
plt.colorbar(shrink=0.9).set_label("sigma", labelpad=-13, y=1.08, rotation=0)
plt.xlabel('Generation',fontsize=7)
plt.ylabel('N Nodes',fontsize=7)


#get proportion coalesced trees by generation
os.chdir("/Users/cj/spaceness/sims/slimout/10kouts/")
trees=[x for x in os.listdir() if not x.startswith(".")]
trees=os.listdir()
sigmas=[float(re.sub("sigma_|_.trees*|.trees","",x)) for x in trees]
out=[]
for sigma in np.unique(sigmas):
    #sigma=sigmas[1000]
    sigtrees=[]
    for i in range(len(sigmas)):
        if sigmas[i]==sigma:
            sigtrees.append(trees[i])
    for t in sigtrees:
        gen=float(re.split("trees",t)[1])
        #t=pyslim.load(t)
        t=sample_treeseq(infile=t,
                         outfile="",
                         nSamples=100,
                         sampling="random",
                         recapitate=False,
                         recombination_rate=1e-8,
                         write_to_file=False,
                         seed=12345)
        roots_per_tree=[x.num_roots for x in t.trees()]
        prop_coalesced=len([x for x in roots_per_tree if x==1])/len(roots_per_tree)
        out.append([sigma,gen,prop_coalesced])
out=np.array(out)

#plot parameters
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False
plt.rcParams.update({'font.size': 6})

#plot proportion coalesced trees by generation
fig = plt.figure(figsize=(2.5,2),dpi=600)
plt.scatter(x=out[:,1]/4,y=out[:,2],c=out[:,0],s=15,alpha=0.7,
            edgecolors='black',cmap=plt.cm.get_cmap("RdYlBu"),
            linewidths=0.25,clip_on=False)
plt.axis([0,5e4,0,1])
plt.xticks(rotation=90)
plt.colorbar(shrink=0.9).set_label("sigma", labelpad=-13, y=1.08, rotation=0)
plt.xlabel('Generation',fontsize=7)
plt.ylabel('Proportion Coalesced\nGene Trees',fontsize=7)
fig.savefig("/Users/cj/Desktop/test.pdf",bbox_inches='tight')
