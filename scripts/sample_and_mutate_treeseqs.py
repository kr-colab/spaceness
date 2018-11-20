from slimtools import *

#check treeseq directory for coalescence
treeseqs,out=check_treeseq_coalescence("/Users/cj/spaceness/sims/slimout/spatial/W16")
np.mean(out[0]) #mean proportion uncoalesced gene trees
np.mean(out[1]) #median roots per gene tree

sample_treeseq_directory(indir="/Users/cj/spaceness/sims/slimout/random_mating/W16",
                         outdir="/Users/cj/spaceness/sims/sampled/random_mating/W16",
                         nSamples=50,
                         recapitate=False,
                         recombination_rate=1e-9)

mutate_treeseqs("/Users/cj/spaceness/sims/sampled/random_mating/W16/",
                "/Users/cj/spaceness/sims/mutated/random_mating/W16/",
                1e-8)

haps,positions,labels,locs=get_ms_outs("/Users/cj/spaceness/sims/mutated/random_mating/W16",False)

positions=discretize_snp_positions(positions)

popsizes=get_pop_sizes(indir="/Users/cj/spaceness/sims/slimout/random_mating/W16/",
              outpath="/Users/cj/spaceness/sumstats/popsize_randmates_W16.txt")

ss=getHaplotypeSumStats(haps,positions,labels,locs,
                        outfile="/Users/cj/spaceness/sumstats/ss_spatial_W16_sample2.txt",
                        verbose=False)
