def sample_treeseq_directory(indir,
                             type="random",
                             outdir,
                             nSamples,
                             recapitate,
                             recombination_rate,
                             n_localities,
                             width):
    '''
    loop sample_treeseq over a directory.
    '''
    trees=os.listdir(indir)
    for i in range(len(trees)):
        ts=pyslim.load(os.path.join(indir,trees[i]))
        if(type=="random"):
            sample_inds=np.unique([ts.node(j).individual for j in ts.samples()]) #final-generation individuals
            subsample=np.random.choice(sample_inds,nSamples,replace=False) #get nSamples random sample inds
            subsample_nodes=[ts.individual(x).nodes for x in subsample] #node numbers for sampled individuals
            subsample_nodes=[a for b in subsample_nodes for a in b] #flatten the list
            subsample_nodes=np.sort(np.array(subsample_nodes))
            o=ts.simplify(subsample_nodes)
        else if(type=="grid"):
            sampling_locs=width/n_localities

        if(recapitate):
            ts=ts.recapitate(recombination_rate=recombination_rate)
        o.dump(os.path.join(outdir,trees[i]))

    return None


indir="/Users/cj/spaceness/sims/slimout/k5/"
type="grid"
nSamples=50
recapitate=True
recombination_rate=1e-9
n_localities=4
width=35

sample_treeseq_directory()
