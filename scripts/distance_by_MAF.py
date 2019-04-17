from slimtools import *
import sklearn.metrics
from tqdm import tqdm
infiles=["/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.21467883716606018_.trees_3812242",
         "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.400670106952987_.trees_1541565",
         "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.8863583600140017_.trees_3937856",
         "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_1.257878563742616_.trees_1186380",
         "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_2.0108173412933894_.trees_7501871",
         "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_4.038506335646773_.trees_9162471"]
for j in infiles:
    simname=os.path.basename(j)
    label=float(re.split("_",simname)[1])
    ts=sample_treeseq(infile=j,
                      outfile="",
                      nSamples=1000,
                      sampling="random",
                      recapitate=False,
                      recombination_rate=1e-9,
                      write_to_file=False,
                      sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
                      plot=False,
                      seed=12345)

    #get generation times estimated from short all-individual simulations
    gentimes=np.loadtxt("/Users/cj/spaceness/W50sp_gentimes.txt")
    gentime=[x[0] for x in gentimes if np.round(x[1],5)==np.round(label,5)]
    ts=msp.mutate(ts,1e-8/gentime[0],random_seed=12345)

    #get haplotypes and locations
    haps=ts.genotype_matrix()
    genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
    sample_inds=np.unique([ts.node(j).individual for j in ts.samples()]) #add check that nodes corresponding to individuals are sequential for future versions?
    locs=[[ts.individual(x).location[0],ts.individual(x).location[1]] for x in sample_inds]

    allele_counts=genotypes.count_alleles()
    dafs=np.array(allele_counts[:,1]/2000)
    #average distance among allele copies, by MAF bin
    derived_distances=np.zeros(shape=len(dafs))
    for i in tqdm(range(len(dafs))):
        site=genotypes[i,:,:]
        isderived=site[:,1]>0
        derivedlocs=np.array(locs)[isderived,:]
        if len(derivedlocs)==0:
            derived_distances[i]=0
        else:
            derived_distances[i]=np.mean(sklearn.metrics.pairwise_distances(derivedlocs))

    out=np.transpose(np.array([dafs,derived_distances]))
    np.savetxt(X=out,fname=os.path.join("/Users/cj/spaceness/dist_by_daf/",str(label)+".txt"))
