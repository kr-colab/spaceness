from slimtools import *

#check treeseq directory for coalescence
treeseqs,out=check_treeseq_coalescence("/Users/cj/spaceness/sims/slimout/wf/")
np.mean(out[0]) #mean proportion uncoalesced trees
np.mean(out[1]) #mean roots per tree

sample_treeseq_directory(indir="/Users/cj/spaceness/sims/slimout/wf/",
                         outdir="/Users/cj/spaceness/sims/sampled/wf/",
                         nSamples=50,
                         recapitate=False,
                         recombination_rate=1e-8)

mutate_treeseqs("/Users/cj/spaceness/sims/sampled/wf",
                "/Users/cj/spaceness/sims/mutated/wf",
                4e-8)

haps,positions,labels,locs=get_ms_outs("/Users/cj/spaceness/sims/mutated/wf",False)

positions=discretize_snp_positions(positions)

popsizes=get_pop_sizes(indir="/Users/cj/spaceness/sims/slimout/wf/",
              outpath="/Users/cj/spaceness/sumstats/popsize_wfN1000.txt")

ss=getHaplotypeSumStats(haps,positions,labels,locs,
                        outfile="/Users/cj/spaceness/sumstats/sumstats_wf_N1000.txt",
                        ibs_tracts=True,
                        verbose=True,maxlen=1e8)
