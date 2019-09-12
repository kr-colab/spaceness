from slimtools import *
os.chdir("/Users/cj/spaceness")

#sample and mutate treeseq
def get_ibs_dist(infile,sampling,nSamples,min_len_to_keep):
    np.random.seed(12345)
    #infile="sims/slimout/spatial/W50_run3/sigma_3.05686209988018_.trees_7807882"
    simname=os.path.basename(infile)
    label=float(re.split("_",simname)[1])
    ts=sample_treeseq(infile=infile,
                      outfile="",
                      nSamples=nSamples,
                      sampling=sampling,
                      recapitate=False,
                      recombination_rate=1e-8,
                      write_to_file=False,
                      sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
                      plot=False,
                      seed=12345)

    #get generation times estimated from short all-individual simulations
    gentimes=np.loadtxt("/Users/cj/spaceness/W50sp_gentimes.txt")
    gentime=[x[0] for x in gentimes if np.round(x[1],5)==np.round(label,5)]
    ts=msp.mutate(ts,1e-8/gentime[0],random_seed=12345)
    haps,positions,locs=get_ms_outs(ts)
    positions=discretize_snp_positions(positions)
    maxlen=1e8
    verbose=True

    #load in genotypes etc into scikit-allel
    genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
    allele_counts=genotypes.count_alleles()
    genotype_allele_counts=genotypes.to_allele_counts()

    #nonspatial population-wide summary statistics
    segsites=np.shape(genotypes)[0]
    #mpd=np.mean(allel.mean_pairwise_difference(allele_counts))
    #pi=(mpd*segsites)/n_accessible_sites                               #alternate pi if not all sites are accessible
    pi=allel.sequence_diversity(positions,allele_counts,start=1,stop=1e8)
    tajD=allel.tajima_d(ac=allele_counts,start=1,stop=maxlen) #NOTE check for 0 v 1 positions
    thetaW=allel.watterson_theta(pos=positions,ac=allele_counts,start=1,stop=maxlen)
    het_o=np.mean(allel.heterozygosity_observed(genotypes))
    fis=np.mean(allel.inbreeding_coefficient(genotypes))
    sfs=allel.sfs(allele_counts[:,1])
    sfs=np.append(sfs,[np.repeat(0,np.shape(haps)[1]-len(sfs)+1)])
    #Isolation by distance
    if(verbose):
        print("calculating pairwise statistics")
    gen_dist=allel.pairwise_dxy(pos=positions,
                                gac=genotype_allele_counts,
                                start=0,
                                stop=maxlen)
    sp_dist=np.array(scipy.spatial.distance.pdist(locs))
    gen_sp_corr=np.corrcoef(gen_dist[sp_dist>0],np.log(sp_dist[sp_dist>0]))[0,1]
    #print(gen_sp_corr)
    gen_dist_skew=scipy.stats.skew(gen_dist)
    gen_dist_var=np.var(gen_dist)
    gen_dist_mean=np.mean(gen_dist)

    #IBS tract length summaries
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
            if len(ibspair)>0:
                spat_dists.append(spdist)
                ibs.append(ibspair)
                ibs_mean_pair.append(np.mean(ibspair))
                #ibs_95p_pair.append(np.percentile(ibspair,95))
                #ibs_var_pair.append(np.var(ibspair))
                ibs_skew_pair.append(scipy.stats.skew(ibspair))
                ibs_blocks_over_1e6_pair.append(len([x for x in ibspair if x>1e6]))
                ibs_blocks_pair.append(len(ibspair))

    #ibs stat ~ spatial distance correlations
    ibs_mean_spat_corr=np.corrcoef(ibs_mean_pair,np.log(spat_dists))[0,1] #v noisy - seems to reflect prob of sampling close relatives...
    ibs_1e6blocks_spat_corr=np.corrcoef(ibs_blocks_over_1e6_pair,np.log(spat_dists))[0,1] #better
    ibs_skew_spat_corr=np.corrcoef(ibs_skew_pair,np.log(spat_dists))[0,1] #good
    ibs_blocks_spat_corr=np.corrcoef(ibs_blocks_pair,np.log(spat_dists))[0,1] #best

    #summaries of the full distribution
    ibs_blocks_over_1e6=np.sum(ibs_blocks_over_1e6_pair)
    ibs_blocks_over_1e6_per_pair=np.mean(ibs_blocks_over_1e6_pair)
    ibs_blocks_per_pair=np.mean(ibs_blocks_pair)
    ibs_flat=np.array([x for y in ibs for x in y])
    ibs_mean=np.mean(ibs_flat)
    ibs_var=np.var(ibs_flat)
    ibs_skew=scipy.stats.skew(ibs_flat)

    return ibs,spat_dists

# ibs_dists=[]
infiles=[#"/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_2.037678004729399_.trees_8681918",
         "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.21467883716606018_.trees_3812242",
         #"/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.5948459719801498_.trees_1208057",
         "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.400670106952987_.trees_1541565",
         #"/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.3027443297023636_.trees_3141112",
         "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_1.0034461761403148_.trees_9532665",
         ]

# for i in infiles:
#     ibs_dists.append(get_ibs_dist(i))
plt.rcParams.update({'font.size': 6,
                     "figure.autolayout":True})
fig = plt.figure(figsize=(3,2.5),dpi=200)
fig.subplots_adjust(hspace=0.3, wspace=0.3)
for i in range(len(infiles)):
    ibs,dists=get_ibs_dist(infiles[i],
                           "random",
                           40,
                           2)
    ibs=np.array(ibs)

    #cumulative distributions
    farIBS=ibs[np.array(dists)>48]
    farIBS=np.array([x for y in farIBS for x in y])
    farIBSx=np.sort(farIBS)[::-1]
    farIBSy=np.array(range(len(farIBS)))/np.sum(np.array(dists)>48)#float(len(farIBS))
    closeIBS=ibs[np.array(dists)<2]
    closeIBS=np.array([x for y in closeIBS for x in y])
    closeIBSx=np.sort(closeIBS)[::-1]
    closeIBSy=np.array(range(len(closeIBS)))/np.sum(np.array(dists)<2)#float(len(closeIBS))

    np.savetxt(X=np.transpose([closeIBSx,closeIBSy]),fname="/Users/cj/Desktop/"+os.path.basename(infiles[i])+".close.IBScdist")
    np.savetxt(X=np.transpose([farIBSx,farIBSy]),fname="/Users/cj/Desktop/"+os.path.basename(infiles[i])+".far.IBScdist")

    ax1=fig.add_subplot(1,1,1)
    ax1.plot(closeIBSx,closeIBSy,c="red",linewidth=0.5,label="dist<2",linestyle=["solid","dashed","dotted"][i])
    ax1.plot(farIBSx,farIBSy,c="blue",linewidth=0.5,label="dist>48",linestyle=["solid","dashed","dotted"][i])
    ax1.legend()
    ax1.set_title(r"$\sigma=$"+str(np.round(float(re.split("_",os.path.basename(infiles[i]))[1]),2)))
    ax1.set_xlabel("IBS tract length")
    ax1.set_ylabel("n IBS tracts > X per pair of samples")
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    fig.savefig(fname="/Users/cj/Desktop/ibs_cumdists2.pdf")

    # #plot
    # ax1=fig.add_subplot(1,3,i+1)
    # ax1.plot(closeIBSx,closeIBSy,c="red",linewidth=0.5,label="dist<2")
    # ax1.plot(farIBSx,farIBSy,c="blue",linewidth=0.5,label="dist>48")
    # ax1.legend()
    # ax1.set_title(r"$\sigma=$"+str(np.round(float(re.split("_",os.path.basename(infiles[i]))[1]),2)))
    # ax1.set_xlabel("IBS tract length")
    # ax1.set_ylabel("n IBS tracts > X per pair of samples")
    # ax1.set_yscale("log")
    # ax1.set_xscale("log")
    # fig.savefig(fname="/Users/cj/Desktop/ibsfig_prob.pdf")
