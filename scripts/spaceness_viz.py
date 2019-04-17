from slimtools import *
from numpy.polynomial.polynomial import polyfit
maxlen=1e8
ts=sample_treeseq(infile="/Users/cj/spaceness/sims/slimout/spatial/W35/sigma_1.065973462406566_.trees",
               outfile="",
               nSamples=50,
               recapitate=False,
               recombination_rate=1e-9,
               write_to_file=False,
               sampling="random",
               sampling_locs=[[8,8],[27,27],[27,8],[8,27]],
               plot=True,
               seed=0)

ts=msp.mutate(ts,1e-8,random_seed=12345)
haps,pos,locs=get_ms_outs(ts)
positions=discretize_snp_positions(pos)
genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
genotypes=genotypes[:100000,:]

#heterozygosity by spatial position
ind_hets=[]
for i in tqdm(range(np.shape(genotypes)[1])):
    ind_hets.append(1-sum([genotypes[:,i][j][0]==genotypes[:,i][j][1] for j in range(np.shape(genotypes)[0])])/np.shape(genotypes)[0])
ind_het_z=[(x-np.mean(ind_hets))/np.std(ind_hets) for x in ind_hets]

plt.hexbin(x=locs[:,0],y=locs[:,1],C=ind_het_z,edgecolors='k',gridsize=17)
plt.viridis()
plt.colorbar().ax.set_title("SDs from Mean\nHeterozygosity")
plt.axis([0,35,0,35])
plt.savefig("/Users/cj/Desktop/het_z_by_ind.pdf")


#heterozygosity ~ distance to midpoint
dists_to_middle=[scipy.spatial.distance.euclidean(x,[17.5,17.5]) for x in locs]
non_zero_het=[x for x in range(len(ind_hets)) if ind_hets[x]>1e-3]
ind_hets_nonzero=np.array(ind_hets)[non_zero_het]
dists_to_middle_nonzero=np.array(dists_to_middle)[non_zero_het]
b,m = polyfit(dists_to_middle_nonzero, ind_hets_nonzero, 1)
R2=np.corrcoef(x=dists_to_middle_nonzero,y=ind_hets_nonzero)[0,1]**2
plt.hexbin(x=np.log(dists_to_middle_nonzero),y=ind_hets_nonzero,gridsize=20,mincnt=1)
plt.colorbar().ax.set_title("N\nSamples")
plt.plot(dists_to_middle_nonzero, b + m * np.array(dists_to_middle_nonzero),'-',c="k")
plt.title("Heterozygosity by Distance to Midpoint\n")
plt.annotate(text="slope = -1.9e-3\n$R^2=$"+str(np.round(R2,3)),xy=(17,0.16))
plt.savefig("/Users/cj/Desktop/het_by_dist_to_midpoint.pdf")


#other summary stats
allele_counts=genotypes.count_alleles()
genotype_allele_counts=genotypes.to_allele_counts()
segsites=np.shape(genotypes)[0]
pi=allel.sequence_diversity(positions,allele_counts,start=1,stop=1e8)
tajD=allel.stats.diversity.tajima_d(ac=allele_counts,start=1,stop=maxlen) #NOTE check for 0 v 1 positions
thetaW=allel.stats.diversity.watterson_theta(pos=positions,ac=allele_counts,start=1,stop=maxlen)
het_o=np.mean(allel.heterozygosity_observed(genotypes))
fis=np.mean(allel.inbreeding_coefficient(genotypes))
sfs=allel.sfs(allele_counts[:,1]) #last entry seems to be the highest *non-zero* SFS entry (wtf?) so adding zeros
sfs=np.append(sfs,[np.repeat(0,np.shape(haps)[1]-len(sfs)+1)])
gen_dist=allel.pairwise_dxy(pos=positions,
                            gac=genotype_allele_counts,
                            start=0,stop=maxlen)
sp_dist=np.array(scipy.spatial.distance.pdist(locs))
gen_sp_corr=np.corrcoef(gen_dist,np.log(sp_dist))[0,1]
gen_dist_skew=scipy.stats.skew(gen_dist)
gen_dist_var=np.var(gen_dist)
gen_dist_mean=np.mean(gen_dist)

plt.scatter(np.log(sp_dist),gen_dist,s=80,facecolors='none',edgecolors='k')
plt.axis([0,4,0.0001,0.000225])

min_len_to_keep=2
pairs=itertools.combinations(range(np.shape(haps)[1]),2)
ibs=[];spat_dists=[];ibs_mean_pair=[];ibs_95p_pair=[]
ibs_var_pair=[];ibs_skew_pair=[];ibs_blocks_over_1e6_pair=[]
ibs_blocks_pair=[]
locs2=np.repeat(locs,2,0)
for j in tqdm(pairs):
    spdist=scipy.spatial.distance.euclidean(locs2[j[0]],locs2[j[1]])
    if spdist>0:
        ibspair=getPairwiseIbsTractLengths(haps[:,j[0]],
                                              haps[:,j[1]],
                                              positions,
                                              maxlen,
                                              min_len_to_keep)
        if len(ibspair>0):
            spat_dists.append(spdist)
            ibs.append(ibspair)
            ibs_mean_pair.append(np.mean(ibspair))
            #ibs_95p_pair.append(np.percentile(ibspair,95))
            #ibs_var_pair.append(np.var(ibspair))
            ibs_skew_pair.append(scipy.stats.skew(ibspair))
            ibs_blocks_over_1e6_pair.append(len([x for x in ibspair if x>1e6]))
            ibs_blocks_pair.append(len(ibspair))

#ibs stat ~ spatial distance correlations
#ibs_mean_spat_corr=np.corrcoef(ibs_mean_pair,np.log(spat_dists))[0,1] #v noisy - seems to reflect prob of sampling close relatives...
ibs_1e6blocks_spat_corr=np.corrcoef(ibs_blocks_over_1e6_pair,np.log(spat_dists))[0,1] #better
ibs_skew_spat_corr=np.corrcoef(ibs_skew_pair,np.log(spat_dists))[0,1]
ibs_blocks_spat_corr=np.corrcoef(ibs_blocks_pair,np.log(spat_dists))[0,1]

#summaries of the full distribution
ibs_blocks_over_1e6=np.sum(ibs_blocks_over_1e6_pair)
ibs_flat=np.array([x for y in ibs for x in y])
ibs_95p=np.percentile(ibs_flat,95)
ibs_mean=np.mean(ibs_flat)
ibs_var=np.var(ibs_flat)
ibs_skew=scipy.stats.skew(ibs_flat)

plt.scatter(np.log(spat_dists),ibs_blocks_pair,s=80,facecolors='none',edgecolors='k')
