import sys, msprime as msp, numpy as np
import multiprocessing as mp, subprocess as sp
import pyslim, os, io, shutil
import shlex, scipy, time, re
from shutil import copyfile
import allel
import time
from scipy import spatial
from scipy import stats
from tqdm import tqdm
import itertools
import matplotlib.pyplot as plt

def run_one_slim_sim(sampled_param,
                     min,
                     max,
                     slim_path,
                     slim_recipe,
                     outdir):
    '''
    Run one SLiM simulation pulling parameter values from a uniform
    distribution. Output will be named with the name and value of
    the sampled parameter.

    See run_slim.py for a version of this function that accepts command-line
    arguments.

    '''

    #get sampled param names and values
    val=np.random.uniform(min,max)

    #get output file path
    filename = sampled_param+"_"+str(val)+"_.trees"
    filepath = os.path.join(outdir,filename)

    #set strings for defining SLiM variables
    label_str=sampled_param+"="+str(val)
    output_str = "outpath='"+str(filepath)+"'"

    #format params for subprocess.check_output
    command=[slim_path,
             "-d",output_str,
             "-d",label_str,
             slim_recipe]

    #run it
    print("starting slim simulation with "+args.sampled_param+"="+str(val))
    sp.check_output(command)

    return None

def check_treeseq_coalescence(indir):
    '''
    Read in a set of SLiM tree sequences and record (1) the
    proportion of uncoalesced gene trees and (2) the mean
    number of remaining lineages at the root.
    '''
    treeseqs=[f for f in os.listdir(indir) if not f.startswith(".")]
    out=np.zeros(shape= (2,len(treeseqs)))
    j=0
    for t in tqdm(treeseqs):
        ts=pyslim.load(os.path.join(indir,t))
        sts=ts.simplify()
        nroots=[0 for _ in range(sts.num_trees)]
        i=0
        for t in sts.trees():
            nroots[i]=t.num_roots
            i=i+1
        prop_uncoalesced=np.sum([x>1 for x in nroots])/len(nroots)
        mean_lineages=np.mean(nroots)
        out[:,j]=[prop_uncoalesced,mean_lineages]
        j=j+1
    return treeseqs,out

def get_pop_sizes(indir,outpath):
    trees = trees=[f for f in os.listdir(indir) if not f.startswith(".")]
    labels=[float(re.sub("sigma_|_.trees*|.trees","",x)) for x in trees]
    out=np.empty(shape=(len(trees),2))
    for i in tqdm(range(len(trees))):
        ts=pyslim.load(os.path.join(indir,trees[i]))
        popsize=ts.num_samples/2
        out[i]=[labels[i],popsize]
    np.savetxt(outpath,out)
    return out

def sample_treeseq_directory(indir,outdir,nSamples,recapitate,recombination_rate):
    '''
    Sample individuals from a tree sequence then simplify and (optionally)
    recapitate the tree with msprime.
    '''
    trees=[f for f in os.listdir(indir) if not f.startswith(".")]
    #pop_sizes=[]
    for i in tqdm(range(len(trees))):
        ts=pyslim.load(os.path.join(indir,trees[i]))
        #if(output_pop_sizes):
            #pop_sizes.append(ts.num_samples/2)
        sample_inds=np.unique([ts.node(j).individual for j in ts.samples()]) #final-generation individuals
        subsample=np.random.choice(sample_inds,nSamples,replace=False) #get nSamples random sample inds
        subsample_nodes=[ts.individual(x).nodes for x in subsample] #node numbers for sampled individuals
        subsample_nodes=[a for b in subsample_nodes for a in b] #flatten the list
        subsample_nodes=np.sort(np.array(subsample_nodes))
        o=ts.simplify(subsample_nodes)
        if(recapitate):
            ts=ts.recapitate(recombination_rate=recombination_rate)
        o.dump(os.path.join(outdir,trees[i]))
    #if(output_pop_sizes):
    #    np.savetxt(pop_size_outpath,pop_sizes)
    return None

def mutate_treeseqs(indir,outdir,mu):
    '''
    Add mutations at constant rate "mu" to all tree sequences in directory "indir"
    '''
    trees=[f for f in os.listdir(indir) if not f.startswith(".")]
    for i in tqdm(trees):
        ts=pyslim.load(os.path.join(indir,i))
        ts=msp.mutate(ts,mu)
        ts.dump(os.path.join(outdir,i))

    return None

def get_ms_outs(direc,subset,start=0,stop=0):
    '''
    loop through a directory of mutated tree sequences
    and return the haplotype matrices, positions, labels,
    and spatial locations as four numpy arrays.
    If subset is true, stop and start are the
    range of files to process (0-indexed, stop not included,
    as ordered by "os.listdir(direc)").
    '''
    haps = []
    positions = []
    locs = []

    trees=[f for f in os.listdir(direc) if not f.startswith(".")]
    labels=[float(re.sub("sigma_|_.trees*|.trees","",x)) for x in trees]
    if(subset):
        trees=trees[start:stop]
        labels=labels[start:stop]
    for i in tqdm(trees):
        ts = pyslim.load(os.path.join(direc,i))
        haps.append(ts.genotype_matrix())
        positions.append(np.array([s.position for s in ts.sites()]))
        sample_inds=np.unique([ts.node(j).individual for j in ts.samples()])
        locs.append([[ts.individual(x).location[0],
               ts.individual(x).location[1]] for x in sample_inds])

    haps = np.array(haps)
    positions = np.array(positions)
    locs=np.array(locs)

    return haps,positions,labels,locs

def discretize_snp_positions(ogpositions):
    '''
    Takes an array of SNP positions as floats (ie from msprime) and returns
    integer positions. If two SNPs fall in the same integer position, one is
    shifted to the right by one base pair.
    '''
    count=0
    newpositions=[]
    for i in tqdm(range(len(ogpositions))):
        dpos=[int(x) for x in ogpositions[i]]
        for j in range(len(dpos)-1):
            if(dpos[j]==dpos[j+1]):
                dpos[j+1]=dpos[j+1]+1
                count=count+1
        newpositions.append(np.sort(dpos)) #NOTE:shouldn't need to sort here but one alignment threw an error (i=1 for the 200 sim set) that indices were not monotonically increasing for unknown reasons. Investigate further...
    print(str(count)+" SNPs were shifted one base pair")
    return(np.array(newpositions))

def getPairwiseIbsTractLengths(x,y,positions,maxlen):
    '''
    input:
    x: haplotype as 1D array
    y: haplotype as 1D array
    positions: a 1D array of SNP positions (integer, 1-based)

    Returns:
    1d array listing distances between adjacent SNPs in a pair of sequences.

    Example:
    haps,positions,labels,locs=getHapsPosLabelsLocs(direc="mutated_treeSeq_directory")
    ibs_lengths_1_2=getPairwiseIbsTractLengths(haps[0][:,0],haps[0][:,1],positions[0])
    '''
    snps=~np.equal(x,y)
    snp_positions=positions[snps]
    l=len(snp_positions)
    ibs_tracts=[]
    if(l==0):
        ibs_tracts=[maxlen]
    else:
        if(l>1):
            ibs_tracts=snp_positions[np.arange(1,l,2)]-snp_positions[np.arange(0,l-1,2)] #middle blocks
        np.append(ibs_tracts,snp_positions[0])          #first block
        np.append(ibs_tracts,maxlen-snp_positions[l-1]) #last block
    return ibs_tracts

def getHaplotypeSumStats(haps,positions,labels,locs,outfile,verbose=True):
    for i in tqdm(range(len(haps))):
        if(verbose):
            print("processing simulation "+str(i))

        #load in genotypes etc into scikit-allel
        genotypes=allel.HaplotypeArray(haps[i]).to_genotypes(ploidy=2)
        allele_counts=genotypes.count_alleles()
        genotype_allele_counts=genotypes.to_allele_counts()

        #nonspatial population-wide summary statistics
        segsites=np.shape(genotypes)[0]
        mpd=np.mean(allel.diversity.mean_pairwise_difference(allele_counts))
        pi=(mpd*segsites)/1e8
        tajD=allel.diversity.tajima_d(ac=allele_counts,start=1,stop=1e8)
        thetaW=allel.diversity.watterson_theta(pos=positions[i],ac=allele_counts,start=1,stop=1e8)
        het_o=np.mean(allel.heterozygosity_observed(genotypes))
        fis=np.mean(allel.inbreeding_coefficient(genotypes))
        sfs=allel.sfs(allele_counts[:,1]) #last entry seems to be the highest *non-zero* SFS entry (wtf?) so adding zeros
        sfs=np.append(sfs,[np.repeat(0,np.shape(haps[i])[1]-len(sfs)+1)])

        #pairwise summary stats
        if(verbose):
            print("calculating pairwise statistics")
        gen_dist=allel.pairwise_dxy(pos=positions[i],gac=genotype_allele_counts,start=0,stop=1e8) #SLOOOOOW - issue with noninteger positions?
        sp_dist=np.array(scipy.spatial.distance.pdist(locs[i]))

        #IBS tract length summaries
        if(verbose):
            print("calculating IBS tract lengths")
        pairs=itertools.combinations(range(np.shape(haps[i])[1]),2)
        ibs=[]
        for j in pairs:
            ibs.append(getPairwiseIbsTractLengths(haps[i][:,j[0]],
                                                         haps[i][:,j[1]],
                                                         positions[i],
                                                         1e8))
        ibs_flat=[x for y in ibs for x in y]
        ibs_95p=np.percentile(ibs_flat,95)
        ibs_over_1e6=len([x for x in ibs_flat if x>=1e6])
        ibs_num_blocks_per_pair=len(ibs_flat)/scipy.special.comb(np.shape(haps[i])[1],2)
        ibs_binned=scipy.stats.binned_statistic(x=ibs_flat,
                                                values=ibs_flat,
                                                statistic='count',
                                                bins=np.arange(0,1e8,2e4))

        #summaries of pairwise stats as a function of geographic distance
        gen_sp_corr=np.corrcoef(gen_dist,sp_dist)[0,1]

        #row=np.append(sp_dist,gen_dist,sfs)
        row=[labels[i],segsites,pi,thetaW,tajD,het_o,fis,gen_sp_corr,ibs_95p,ibs_num_blocks_per_pair,ibs_over_1e6]
        row=np.append(row,sfs)
        row=np.append(row,ibs_binned[0])

        if(i==0):
            out=np.zeros(shape=(len(haps),len(row)))
            out[i]=row
        else:
            out[i]=row
        if(outfile):
            np.savetxt(outfile,out)

    return(out)
################ example pipeline use ####################

#check treeseq directory for coalescence
# treeseqs,out=check_treeseq_coalescence("/Users/cj/spaceness/sims/slimout/spatial/W16")
# np.mean(out[0]) #mean proportion uncoalesced gene trees
# np.mean(out[1]) #median roots per gene tree
#
# sample_treeseq_directory(indir="/Users/cj/spaceness/sims/slimout/spatial/W16",
#                          outdir="/Users/cj/spaceness/sims/sampled/spatial/W16",
#                          nSamples=50,
#                          recapitate=True,
#                          recombination_rate=1e-9,
#                          output_pop_sizes=True,
#                          pop_size_outpath="/Users/cj/spaceness/sumstats/popsize_W16_rm.txt")
#
# mutate_treeseqs("/Users/cj/spaceness/sims/sampled/spatial/W16/",
#                 "/Users/cj/spaceness/sims/mutated/spatial/W16/",
#                 1e-8)
#
#haps,positions,labels,locs=get_ms_outs("/Users/cj/spaceness/sims/mutated/spatial/W16",
#                                        False)
#
#positions=discretize_snp_positions(positions)
#
# ss=getHaplotypeSumStats(haps,positions,labels,locs,
#                         outfile="/Users/cj/spaceness/sumstats/ss_spatial_W16.txt",
#                         verbose=False)





# #output proportion uncoalesced & dispersal for 10k outs
# np.savetxt("/Users/cj/spaceness/uncoal_prop.txt",np.transpose(out))
# treeseqs=[f for f in os.listdir("/Users/cj/spaceness/sims/slimout/10kouts/") if not f.startswith(".")]
# tmp=open("/Users/cj/spaceness/uncoal_tsnames.txt","w")
# for i in range(len(treeseqs)):
#     tmp.write(treeseqs[i]+"\n")
# tmp.close()
#
# #tree names to get generations for 10k outs
# treeseqs=os.listdir("/Users/cj/spaceness/sims/10kouts")
# o=open("/Users/cj/spaceness/10kout_treenames.txt")
# for i in range(len(treeseqs)):
#     o.write(treeseqs[i]+"\n")
# o.close()

###################################################
## old pipeline using multiple presampled param values
# def sample_sim_params(params_to_sample,
#                       min,
#                       max,
#                       nreps,
#                       outfile=False):
#     '''
#     Sample parameter values from uniform distributions and write to file.
#     params_to_sample, min, and max should be arrays, even if length 1
#     (eg min=[0.2]).
#     '''
#     vals=np.zeros((len(params_to_sample),nreps))
#     for i in range(len(params_to_sample)):
#         for j in range(nreps):
#             val=np.random.uniform(min[i],max[i])
#             vals[i][j]=val
#     if(outfile):
#         header=" ".join(params_to_sample)
#         if(os.path.exists(outfile)):
#             print("warning! Outfile already exists. Delete and try again.")
#         np.savetxt(outfile,vals,header=header)
#     return vals

## old version of run_one_slim_sim that pulls multiple sampled param values from an existing file
# def run_one_slim_sim(sim_num,
#                      sampled_params,
#                      slim_path,
#                      slim_recipe,
#                      outdir):
#     '''
#     Run one SLiM simulation pulling parameter values from the text file
#     at path "sampled_params" (ie the output of sample_sim_params).
#     Output file names are the index (row) of the sampled parameter files.
#     '''
#
#     #get sampled param names and values
#     with open(sampled_params) as f:
#         params_to_sample = f.readline()
#     params_to_sample=re.sub("\n|# ","",params_to_sample)
#     params_to_sample=params_to_sample.split(" ")
#     vals=np.loadtxt(sampled_params)
#
#     #get output file path
#     filename = str(sim_num) + ".trees"
#     filepath = os.path.join(outdir,filename)
#
#     #set strings for defining SLiM variables
#     label_strs=[params_to_sample[i]+"="+str(vals[i][sim_num])
#                 for i in range(len(params_to_sample))]
#     output_str = "outpath='"+str(filepath)+"'"
#
#     #format params for subprocess.check_output
#     command=[slim_path,
#              "-d",output_str]
#     for i in range(len(params_to_sample)):
#         command.append("-d")
#         command.append(label_strs[i])
#     command.append(slim_recipe)
#
#     #run it
#     subprocess.check_output(command)
#
#     return None
