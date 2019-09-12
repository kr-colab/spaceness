#import os;os.chdir("/Users/cj/spaceness/scripts/")
from slimtools import *
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--infile")
parser.add_argument("--outfile")
parser.add_argument("--width",type=int)
parser.add_argument("--pop_config")
args=parser.parse_args()

ts=msp.load(args.infile)
#ts=msp.load("/Users/cj/spaceness/sims/sigma_0.21313291689796326_sim_911380669.trees")
#ts=msp.mutate(ts,1e-8)

#sample 60 individuals
nodes_to_sample=np.random.choice(np.arange(0,499,2),60,replace=False)
node_partners=nodes_to_sample+1
nodes_to_sample=np.concatenate((nodes_to_sample,node_partners))
nodes_to_sample=np.sort(nodes_to_sample)
nodes_to_sample=np.array(nodes_to_sample,dtype="int32")
ts=ts.simplify(nodes_to_sample)
haps=np.array(ts.genotype_matrix())
positions=np.array([s.position for s in ts.sites()])

#get locations etc for all individuals
if args.pop_config=="grid":
    #label=float(re.sub(".trees|sigma_|","",os.path.basename(args.infile)))
    label=float(re.split("_",os.path.basename(args.infile))[1])
    pops=[ts.node(x).population for x in ts.samples()]
    locs=[[(np.floor(x/args.width)+.5)*(50/args.width),((x%args.width)+.5)*(50/args.width)] for x in pops]
    locs=np.array(locs)[np.arange(0,len(locs),2)]
if args.pop_config=="onepop":
    locs=np.array([[np.random.uniform(1,50,1)[0],np.random.uniform(1,50,1)[0]] for x in range(50)]) #random locations for nonspatial sim(?)
    label=float(re.split("_",args.infile)[1])

getSLiMSumStats(haps=haps,
                positions=positions,
                label=label,
                locs=locs,
                outfile=args.outfile,
                maxlen=1e8,
                ibs_tracts=False,
                verbose=True,
                min_len_to_keep=2)

###########~~~~~~~debugging zone~~~~~~~~~###########
# ts=msp.load("/Users/cj/spaceness/sims/msp_grid/0.0010278777692654174.trees")
# ts=msp.mutate(ts,1e-8)
# label=0.06
# pops=[ts.node(x).population for x in ts.samples()]
# locs=[[(np.floor(x/5)+.5)*10,((x%5)+.5)*10] for x in pops]
# locs=np.array(locs)[np.arange(0,len(locs),2)]
# haps=np.array(ts.genotype_matrix())
# positions=np.array([s.position for s in ts.sites()])

# #load in genotypes etc into scikit-allel
# genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
# allele_counts=genotypes.count_alleles()
# genotype_allele_counts=genotypes.to_allele_counts()
# 
# #nonspatial population-wide summary statistics
# segsites=np.shape(genotypes)[0]
# pi=allel.sequence_diversity(positions,allele_counts,start=1,stop=1e8)
# tajD=allel.tajima_d(ac=allele_counts,start=1,stop=maxlen) #NOTE check for 0 v 1 positions
# thetaW=allel.watterson_theta(pos=positions,ac=allele_counts,start=1,stop=maxlen)
# het_o=np.mean(allel.heterozygosity_observed(genotypes))
# fis=np.mean(allel.inbreeding_coefficient(genotypes))
# sfs=allel.sfs(allele_counts[:,1]) #last entry seems to be the highest *non-zero* SFS entry (wtf?) so adding zeros
# sfs=np.append(sfs,[np.repeat(0,np.shape(haps)[1]-len(sfs)+1)])
# 
# gen_dist=allel.pairwise_dxy(pos=positions,
#                             gac=genotype_allele_counts,
#                             start=0,
#                             stop=maxlen)
# sp_dist=np.array(scipy.spatial.distance.pdist(locs))
# np.shape(locs)
# np.shape(genotypes)
# 
# #sp_dist=[x+0.1 for x in sp_dist if x==0] #add 0.1 to avoid 0 distances for intra-population comparisons when taking logs (? better way to deal with this ?)
# gen_sp_corr=np.corrcoef(gen_dist[sp_dist>0],np.log(sp_dist[sp_dist>0]))[0,1]
# gen_dist_skew=scipy.stats.skew(gen_dist)
# gen_dist_var=np.var(gen_dist)
# gen_dist_mean=np.mean(gen_dist)
# 
# #IBS tract length summaries
# pairs=itertools.combinations(range(np.shape(haps)[1]),2)
# ibs=[];spat_dists=[];ibs_mean_pair=[];ibs_95p_pair=[]
# ibs_var_pair=[];ibs_skew_pair=[];ibs_blocks_over_1e6_pair=[]
# ibs_blocks_pair=[]
# locs2=np.repeat(locs,2,0)
# for j in pairs:
#     spdist=scipy.spatial.distance.euclidean(locs2[j[0]],locs2[j[1]])
#     if spdist>0:
#         ibspair=getPairwiseIbsTractLengths(haps[:,j[0]],
#                                            haps[:,j[1]],
#                                            positions,
#                                            1e8,
#                                            2)
#         if len(ibspair)>0:
#             spat_dists.append(spdist)
#             ibs.append(ibspair)
#             ibs_mean_pair.append(np.mean(ibspair))
#             #ibs_95p_pair.append(np.percentile(ibspair,95))
#             #ibs_var_pair.append(np.var(ibspair))
#             ibs_skew_pair.append(scipy.stats.skew(ibspair))
#             ibs_blocks_over_1e6_pair.append(len([x for x in ibspair if x>1e6]))
#             ibs_blocks_pair.append(len(ibspair))
# 
# #ibs stat ~ spatial distance correlations
# ibs_mean_spat_corr=np.corrcoef(ibs_mean_pair,np.log(spat_dists))[0,1]
# ibs_1e6blocks_spat_corr=np.corrcoef(ibs_blocks_over_1e6_pair,np.log(spat_dists))[0,1]
# ibs_skew_spat_corr=np.corrcoef(ibs_skew_pair,np.log(spat_dists))[0,1]
# ibs_blocks_spat_corr=np.corrcoef(ibs_blocks_pair,np.log(spat_dists))[0,1]
# 
# #summaries of the full distribution
# ibs_blocks_over_1e6=np.sum(ibs_blocks_over_1e6_pair)
# ibs_blocks_over_1e6_per_pair=np.mean(ibs_blocks_over_1e6_pair)
# ibs_blocks_per_pair=np.mean(ibs_blocks_pair)
# ibs_flat=np.array([x for y in ibs for x in y])
# ibs_mean=np.mean(ibs_flat)
# ibs_var=np.var(ibs_flat)
# ibs_skew=scipy.stats.skew(ibs_flat)
# 
# if os.path.exists(outfile)==False:
#     sfsnames=["sfs_"+str(x)+" " for x in range(np.shape(haps)[1])]
#     sfsnames.append("sfs_"+str(np.shape(haps)[1]))
#     out=open(outfile,"w")
#     out.write("sigma segsites pi thetaW tajD het_o fis gen_dist_mean gen_dist_var gen_dist_skew gen_sp_corr ibs_mean ibs_var ibs_skew ibs_blocks_per_pair ibs_blocks_over_1e6_per_pair ibs_mean_spat_corr ibs_1e6blocks_spat_corr ibs_skew_spat_corr ibs_blocks_spat_corr "+
#     "".join(sfsnames)+"\n")
# row=[label,segsites,pi,thetaW,tajD,het_o,fis,gen_dist_mean,gen_dist_var,gen_dist_skew,gen_sp_corr,
#      ibs_mean,ibs_var,ibs_skew,ibs_blocks_per_pair,ibs_blocks_over_1e6_per_pair,ibs_mean_spat_corr,ibs_1e6blocks_spat_corr,ibs_skew_spat_corr,
#      ibs_blocks_spat_corr]
# row=np.append(row,sfs)
# row=" ".join(map(str, row))
# out=open(outfile,"a")
# out.write(row+"\n")
# out.close()
