import os
os.chdir("/Users/cj/spaceness/scripts/")
from slimtools import *
os.chdir("/Users/cj/spaceness/")
import pandas as pd
#get generation time predictor from all-individual simulations

gentimes=pd.read_csv("/Users/cj/spaceness/W50sp_slimMuts_gentimes.txt",sep=" ",header=None,float_precision='round_trip')

slim_muts=[];msp_muts=[];sigmas=[];msp_muts_unscaled=[]
files=[x for x in os.listdir("/Users/cj/spaceness/sims/slimout/spatial/W50_muts/") if not x.startswith(".")]
for sim in tqdm(files):
  sigma=float(re.split("_",os.path.basename(sim))[1])
  gt=gentimes[np.round(gentimes[0],5)==np.round(sigma,5)][1]
  a=pyslim.load(os.path.join("/Users/cj/spaceness/sims/slimout/spatial/W50_muts/",sim))
  slim_muts.append(a.num_mutations)
  msp_a=msp.mutate(a,1e-8/(gt),keep=False)
  msp_muts.append(msp_a.num_mutations)
  msp_b=msp.mutate(a,1e-8/(np.mean(gentimes[1])),keep=False)
  msp_muts_unscaled.append(msp_b.num_mutations)
  sigmas.append(sigma)

from matplotlib import pyplot as plt
plt.plot([300000,600000],[300000,600000],color="black")
plt.scatter(slim_muts,msp_muts,color="green")
plt.scatter(slim_muts,msp_muts_unscaled,color="red")
plt.show() 
plt.close()

plt.scatter(sigmas,slim_muts,color="red")
plt.scatter(sigmas,msp_muts,color="blue")
plt.scatter(sigmas,msp_muts_unscaled,color="green")
plt.show()
plt.close()


# muts=pyslim.load("/Users/cj/Desktop/sigma_0.793843124378123_.trees_7796460_822000")
# muts.num_mutations
# muts=sample_treeseq(infile=muts,
#                outfile="",
#                nSamples=100,
#                recapitate=False,
#                recombination_rate=1e-8,
#                write_to_file=False,
#                sampling="random",
#                sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
#                plot=False,
#                seed=12345)
# muts.num_mutations
# 
# nomuts=pyslim.load("/Users/cj/spaceness/sims/slimout/random_mating/W50_run3/sigma_0.7857423531825174_.trees_9690340")
# nomuts=msp.mutate(nomuts,1e-8/4.8)
# nomuts.num_mutations
# nomuts=sample_treeseq(infile=nomuts,
#                outfile="",
#                nSamples=60,
#                recapitate=False,
#                recombination_rate=1e-8,
#                write_to_file=False,
#                sampling="random",
#                sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
#                plot=False,
#                seed=12345)
# nomuts.num_mutations
