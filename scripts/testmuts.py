from slimtools import *

muts=pyslim.load("/Users/cj/Desktop/sigma_0.793843124378123_.trees_7796460_822000")
muts.num_mutations
muts=sample_treeseq(infile=muts,
               outfile="",
               nSamples=60,
               recapitate=False,
               recombination_rate=1e-8,
               write_to_file=False,
               sampling="random",
               sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
               plot=False,
               seed=12345)
muts.num_mutations

nomuts=pyslim.load("/Users/cj/spaceness/sims/slimout/random_mating/W50_run3/sigma_0.7857423531825174_.trees_9690340")
nomuts=msp.mutate(nomuts,1e-8/4.8)
nomuts.num_mutations
nomuts=sample_treeseq(infile=nomuts,
               outfile="",
               nSamples=60,
               recapitate=False,
               recombination_rate=1e-8,
               write_to_file=False,
               sampling="random",
               sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
               plot=False,
               seed=12345)
nomuts.num_mutations

#seems fine so wtf is up with pi
