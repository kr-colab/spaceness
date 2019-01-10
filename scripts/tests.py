#tests for slimtools.py functions
from slimtools import *
from ggplot import *
import matplotlib.pyplot as plt

hap=[]
blocks=[]
for i in range(10):
    blocks=append(blocks,np.round(np.random.uniform(1,1e5))))
    append(hap,np.repeat(0,blocks[i])
    append(hap,1)

check_treeseq_coalescence("/Users/cj/Desktop/tmp/")

help(sample_treeseq)

ts=pyslim.load("/Users/cj/Desktop/sigma_0.3078561850313077_.trees2000000.trees")
nodes=[x for x in ts.nodes()]
inds=[x for x in ts.individuals()]
nodes[0]
inds[0]
extant_nodes=[x.id for x in ts.nodes() if x.time>1999980]



ts.num_samples
np.shape(extant_nodes)
help(ts.simplify)
ts=ts.simplify(samples=np.array(extant_nodes,dtype='int32'))


ts=sample_treeseq(ts,'',50,False,1e-9,False)
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

sts.num_individuals
sts.num_samples
ts.slim_generation


help(np.random.choice)

infile="/Users/cj/spaceness/sims/slimout/spatial/W35/sigma_0.9936382939920805_.trees"
maxlen=1e8
ts=sample_treeseq(infile="/Users/cj/spaceness/sims/slimout/spatial/W35/sigma_0.3003441596629088_.trees",
               outfile="",
               nSamples=100,
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
allele_counts=genotypes.count_alleles()
genotype_allele_counts=genotypes.to_allele_counts()
segsites=np.shape(genotypes)[0]
pi=allel.sequence_diversity(positions,allele_counts,start=1,stop=1e8)
tajD=allel.diversity.tajima_d(ac=allele_counts,start=1,stop=maxlen) #NOTE check for 0 v 1 positions
thetaW=allel.diversity.watterson_theta(pos=positions,ac=allele_counts,start=1,stop=maxlen)
het_o=np.mean(allel.heterozygosity_observed(genotypes))
fis=np.mean(allel.inbreeding_coefficient(genotypes))
sfs=allel.sfs(allele_counts[:,1]) #last entry seems to be the highest *non-zero* SFS entry (wtf?) so adding zeros
sfs=np.append(sfs,[np.repeat(0,np.shape(haps)[1]-len(sfs)+1)])
gen_dist=allel.pairwise_dxy(pos=positions,
                            gac=genotype_allele_counts,
                            start=0,stop=maxlen)
sp_dist=np.array(scipy.spatial.distance.pdist(locs))
gen_sp_corr=np.corrcoef(gen_dist,sp_dist)[0,1]
gen_dist_skew=scipy.stats.skew(gen_dist)
gen_dist_var=np.var(gen_dist)
gen_dist_mean=np.mean(gen_dist)

plt.scatter(np.log(sp_dist),gen_dist,s=80,facecolors='none',edgecolors='k')
plt.axis([0,4,0,0.00175])


def sample_treeseq(infile,
                   nSamples,
                   sampling,
                   recapitate=False,
                   recombination_rate=1e-9,
                   write_to_file=False,
                   outfile='',
                   sampling_locs=[],
                   plot=False,
                   seed=0):
    '''
    Samples a set of individuals from a SLiM tree sequence and returns a
    tree sequence simplified for those individuals.

    params:
    infile - string or object. input treesequence file or object
    nSamples - integer. total samples to return
    sampling - 'random' (nSamples random individuals);
               'midpoint' (nSamples closest individuals to the middle of the landscape);
               'point' (nSamples/len(sampling_locs) individuals closest to each location in sampling_locs).
    recapitate - boolean. Should simulations be run to coalescence with msprime? False unless you know *exactly* what you're doing hereself.
    recombination_rate - float. Recombination rate to use for recapitation, in recombination events per base per generation.
    write_to_file - boolean.
    outfile - string. File path for output.
    sampling_locs - list. locations to draw samples from when using sampling='point', as [[x1,y1],[x2,y2]]
    plot - boolean. plot location of samples?
    seed - random seed for all numpy operations. 0 uses the system default.
    '''
    if type(infile)==str:
        ts=pyslim.load(infile)
    else:
        ts=infile
    if(not seed==0):
        np.random.seed(seed)
    if(sampling=="random"):
        sample_inds=np.unique([ts.node(j).individual for j in ts.samples()]) #final-generation individuals
        subsample=np.random.choice(sample_inds,nSamples,replace=False) #get nSamples random sample inds
        subsample_nodes=[ts.individual(x).nodes for x in subsample] #node numbers for sampled individuals
        subsample_nodes=[a for b in subsample_nodes for a in b] #flatten the list
        subsample_nodes=np.sort(np.array(subsample_nodes)) #note this is needed to preserve individual order relative to locations - should build better individual ID verifications later...
        o=ts.simplify(subsample_nodes)
        locs=np.array([[x.location[0],x.location[1]] for x in ts.individuals()])
        olocs=np.array([[x.location[0],x.location[1]] for x in o.individuals()])
        if plot:
            plt.scatter(locs[:,0], locs[:,1],s=80,facecolors='none',edgecolors='k')
            plt.scatter(olocs[:,0],olocs[:,1],s=80,facecolors='r',edgecolors='r')
            #plt.axis([0, 35, 0, 35])
    if(sampling=="midpoint"):
        #get final generation individuals
        extant_inds=np.unique([ts.node(j).individual for j in ts.samples()]) #final-generation individuals
        extant_nodes=[ts.individual(x).nodes for x in extant_inds] #node numbers for final-generation individuals
        extant_nodes=[a for b in extant_nodes for a in b]
        ts=ts.simplify(extant_nodes)
        #sample nSamples individuals proportional to their distance from the midpoint where p(sample)=(1/dist**4)/sum(dists_to_middle)
        locs=np.array([[x.location[0],x.location[1]] for x in ts.individuals()])
        inds=np.array([x.id for x in ts.individuals()])
        middle=[np.mean(locs[:,0]),np.mean(locs[:,1])]
        dists_to_middle=[scipy.spatial.distance.euclidean(locs[x],middle) for x in range(len(locs))]
        weights=[1/x**4 for x in dists_to_middle]
        weights=weights/np.sum(weights)
        subsample=np.random.choice(inds,nSamples,p=weights,replace=False)
        #closest_inds=sorted(range(len(dists_to_middle)), key=lambda e: dists_to_middle[e],reverse=True)[-nSamples*2:] #this works for sampling half the closest nSamples individuals
        #closest_inds=np.random.choice(closest_inds,int(len(closest_inds)/2),replace=False)
        #subsample=inds[closest_inds]
        subsample_nodes=[ts.individual(x).nodes for x in subsample] #node numbers for sampled individuals
        subsample_nodes=[a for b in subsample_nodes for a in b] #flatten the list
        subsample_nodes=np.sort(np.array(subsample_nodes))
        o=ts.simplify(subsample_nodes)
        olocs=np.array([[x.location[0],x.location[1]] for x in o.individuals()])
        if plot:
            plt.scatter(locs[:,0],locs[:,1],s=80,facecolors='none',edgecolors='k')
            plt.scatter(olocs[:,0],olocs[:,1],s=80,facecolors='r',edgecolors='r')
            #plt.axis([0, 35, 0, 35])
    if(sampling=="point"):
        #get final generation individuals
        extant_inds=np.unique([ts.node(j).individual for j in ts.samples()]) #final-generation individuals
        extant_nodes=[ts.individual(x).nodes for x in extant_inds] #node numbers for final-generation individuals
        extant_nodes=[a for b in extant_nodes for a in b]
        ts=ts.simplify(extant_nodes)
        #sample nSamples individuals proportional to their distance from the midpoint where p(sample)=(1/dist**4)/sum(dists_to_point)
        locs=np.array([[x.location[0],x.location[1]] for x in ts.individuals()])
        inds=np.array([x.id for x in ts.individuals()])
        if nSamples%len(sampling_locs) == 0:
            inds_per_pt= int(nSamples/len(sampling_locs))
        else:
            print("Error: nSamples must be evenly divisible by len(sampling_locs).")
            return
        subsample=[]
        for i in sampling_locs:
            dists=[scipy.spatial.distance.euclidean(locs[x],i) for x in range(len(locs))]
            weights=[1/x**4 for x in dists]
            weights=weights/np.sum(weights)
            subsample.append(np.random.choice(inds,inds_per_pt,p=weights,replace=False))
        subsample=[a for b in subsample for a in b]
        subsample_nodes=[ts.individual(x).nodes for x in subsample] #node numbers for sampled individuals
        subsample_nodes=[a for b in subsample_nodes for a in b] #flatten the list
        subsample_nodes=np.sort(np.array(subsample_nodes))
        o=ts.simplify(subsample_nodes)
        olocs=np.array([[x.location[0],x.location[1]] for x in o.individuals()])
        if plot:
            plt.scatter(locs[:,0], locs[:,1],s=80,facecolors='none',edgecolors='k')
            plt.scatter(olocs[:,0],olocs[:,1],s=80,facecolors='r',edgecolors='r')
            #plt.axis([0, 35, 0, 35])
    if(recapitate):
        o=o.recapitate(recombination_rate=recombination_rate)
    if(write_to_file):
        o.dump(outfile)
    return o
