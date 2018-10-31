import numpy as np, multiprocessing as mp, subprocess as sp
import os, sys, re, pyslim, msprime, sklearn, sys#, allel
from sklearn import ensemble
from sklearn import model_selection
#sys.path.insert(0, "/Users/cj/spaceness/") #shouldn't need this unless running via reticulate
from SLiMomatic import *
import time

wd="/proejcts/kernlab/cbattey2/spaceness/"
os.chdir(wd)

#initialize the data generator
dg=SLiMomatic(
    SLiMExecutableFilePath="slim",
    SLiMRecipe="flat_map.slim",
    numReps=200,
    Ne=4e3,
    recombinationRate=1e-9,
    params_to_sample=["sigma_i"],
    priorLows=[0.2],
    priorHighs=[2])

#run simulations in parallel
dg.simulateAndProduceTrees(direc=os.path.join(wd,"sims/trees/"),nProc=28)

#sample individuals and simplify tree sequences
simplifyTreeSequenceDirectory(indir=wd+"sims/trees/k3/",
                              outdir=wd+"sims/sampled/k3/",
                              nSamples=50)

#add mutations
mutateTrees(treesDirec=wd+"/sims/sampled/k3/",
            outputDirec=wd+"/sims/mutated/k3/",
            numMuts=1,muLow=1e-8,muHigh=1e-8)

#get summary stats from haplotype matrices (probably not working yet)
haps,positions,labels,locs=getHapsPosLabelsLocs(direc=wd+"sims/mutated/k10/")

#discretize snp positions
def discretize_snp_positions(ogpositions):
    '''
    Takes an array of SNP positions as floats (ie from msprime) and returns
    integer positions. If two SNPs fall in the same integer position, one is
    shifted to the right by one base pair.
    '''
    count=0
    newpositions=[]
    for i in range(len(ogpositions)):
        dpos=[int(x) for x in ogpositions[i]]
        for j in range(len(dpos)-1):
            if(dpos[j]==dpos[j+1]):
                dpos[j+1]=dpos[j+1]+1
                count=count+1
        newpositions.append(np.sort(dpos)) #NOTE:shouldn't need to sort here but one alignment threw an error (i=1 for the 200 sim set) that indices were not monotonicall increasing for unknown reasons. Investigate further...
    print(str(count)+" SNPs were shifted one base pair")
    return(np.array(newpositions))

dpositions=discretize_snp_positions(positions)

help(allel.mean_pairwise_difference)
help(allel.diversity.tajima_d)
def getSumStatsFromHaplotypes(haps,positions,labels,locs,outfile,verbose=True):
    out=np.zeros(shape=(len(haps),6+101))
    for i in range(len(haps)):
        if(verbose):
            print("processing simulation "+str(i))
        #load in genotypes etc into scikit-allel
        genotypes=allel.HaplotypeArray(haps[i]).to_genotypes(ploidy=2)
        allele_counts=genotypes.count_alleles()

        #nonspatial population-wide summary statistics
        segsites=np.shape(genotypes)[0]
        mpd=np.mean(allel.diversity.mean_pairwise_difference(allele_counts))
        mpd=(mpd*segsites)/1e8
        tajD=allel.diversity.tajima_d(ac=allele_counts,start=1,stop=1e8)
        thetaW=allel.diversity.watterson_theta(pos=positions[i],ac=allele_counts,start=1,stop=1e8)
        het_o=np.mean(allel.heterozygosity_observed(genotypes))
        fis=np.mean(allel.inbreeding_coefficient(genotypes))
        sfs=allel.sfs(allele_counts[:,1]) #last entry seems to be the highest *non-zero* SFS entry (wtf?) so adding zeros
        sfs=np.append(sfs,[np.repeat(0,np.shape(haps[i])[1]-len(sfs)+1)])

        #pairwise summary stats
        #gen_dist=allel.pairwise_dxy(pos=positions[i],gac=genotypes,start=0,stop=1e8) #SLOOOOOW - issue with noninteger positions?
        #sp_dist=np.array(scipy.spatial.distance.pdist(locs[i]))

        #summaries of pairwise stats as a function of geographic distance
        #gen_sp_corr=np.corrcoef(gen_dist,sp_dist)[0,1]

        #row=np.append(sp_dist,gen_dist,sfs)
        row=[segsites,mpd,thetaW,tajD,het_o,fis]
        row=np.append(row,sfs)

        out[i]=row
        if(outfile):
            np.savetxt(outfile,out)

    return(out)

out=getSumStatsFromHaplotypes(haps,dpositions,labels,locs,
                              outfile="/Users/cj/spaceness/sims/sumstats_k10.txt")

#split training and testing
train_stats=out[0:150]
train_labels=labels[0:150]
test_stats=out[150:200]
test_labels=labels[150:200]

#train extra trees regressor
param_grid_forest = {"max_depth": [5, 10, None],
                     "min_samples_split": [2, 3, 10],
                     "min_samples_leaf": [1, 3, 10]}
clf = sklearn.ensemble.ExtraTreesRegressor(n_estimators=1000,bootstrap=False)
clf_grid = sklearn.model_selection.GridSearchCV(clf,
                                            param_grid=param_grid_forest,
                                            cv=3,n_jobs=-1)
clf_grid.fit(train_stats, train_labels)
clf_grid.best_params_
pytheas=clf_grid.best_estimator_
r2=pytheas.score(test_stats,test_labels)
pred=pytheas.predict(test_stats)
pred_out=open("predictions.txt","w")
for row in pred:
    pred_out.write("%s\n" % row)
pred_out.close()
print(r2)
