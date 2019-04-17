#estimate generation times from short slim runs
import sys,subprocess as sp,re,os,numpy as np,pyslim
from tqdm import tqdm
import argparse
import spatial_slim_ts as sps

parser=argparse.ArgumentParser()
parser.add_argument("--infile",dest="infile",type=str)
args=parser.parse_args()

args=argparse.Namespace(infile="/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.25_.trees_1661344")

sigma=float(re.split("_",os.path.basename(args.infile))[1])

id=int(np.random.uniform(1e7,1e8-1,1))
sp.run("cd ~/spaceness/; slim -d sigma="+str(sigma)+" -d id="+str(id)+" slim_recipes/gen_timer_spatial.slim",
       shell=True)
#ts=pyslim.load("/Users/cj/spaceness/tmp"+str(id)+".trees")
ts = sps.SpatialSlimTreeSequence(pyslim.load("/Users/cj/spaceness/tmp"+str(id)+".trees"), dim=2)

############ average generation time across all individuals
node_times = ts.tables.nodes.time
parent_nodes = ts.tables.edges.parent
child_nodes = ts.tables.edges.child
max_T = 200
use_these = (node_times[parent_nodes] < max_T)
parent_ages = node_times[parent_nodes] - node_times[child_nodes]
genealogical_gen_time=np.mean(parent_ages[use_these])

############ variance in reproductive output
mid_node_ids=np.array([x.id for x in ts.nodes() if x.time in range(100,101)])


np.var(n_offspring)
np.mean(n_offspring)


############# 'effective' generation time sampling inds proportional to final-generation ancestry
alive_nodes=[x.id for x in ts.nodes() if x.time==0] #final-generation nodes
#alive = ts.individuals_alive(0)
#alive_nodes = ts.individual_nodes(np.where(alive)[0])
node_ancestry = ts.proportion_ancestry_nodes(alive_nodes)[0]
ind_ancestry = np.zeros(ts.num_individuals)
for k, ind in enumerate(ts.individuals()):
    ind_ancestry[k] = sum(node_ancestry[ind.nodes])

# now find average parent age for everyone
parents = ts.individual_parents_dict()
mean_parent_ages = np.zeros(ts.num_individuals) - 1.0
for ind in parents:
    ages = [ts.individual(x).time - ts.individual(ind).time for x in set(parents[ind])]
    if len(ages) == 2:
        mean_parent_ages[ind] = np.mean(ages)
# the max generation time is this long -
# we need to avoid looking at generations where we don't have their parents recorded
max_gen_time = 20
parent_age_by_time = np.zeros(num_gens - max_gen_time)
for t in range(num_gens - max_gen_time):
    alive = np.logical_and(ts.individuals_alive(t), mean_parent_ages > 0)
    probs = ind_ancestry * alive
    probs = probs / sum(probs)
    parent_age_by_time[t] = sum(probs * mean_parent_ages)

print("time_ago eff_generation_time")
for t in range(num_gens - max_gen_time):
    print("{} {}", t, parent_age_by_time[t])



#write to file
#out=open("/home/cbattey2/spaceness/W50rm_gentimes.txt","a")
#out.write(str(row[0])+" "+str(row[1])+"\n")
#out.close()

sp.run("cd /home/cbattey2/spaceness/;rm /home/cbattey2/spaceness/tmp"+str(id)+".trees",
       shell=True)
