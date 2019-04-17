from slimtools import *
# import argparse
# parser=argparse.ArgumentParser()
# parser.add_argument("--infile")
# parser.add_argument("--outfile")
# args=parser.parse_args()


#######################utility functions####################################
##migration matrix for square landscapes with or without migration along diagonals
def get_grid_migration_matrix(width,rate,diagonal):
    mig_matrix=np.zeros((width*width,width*width))
    lefts=np.arange(0,width**2,width)
    rights=np.arange(width-1,width**2,width)
    for i in np.arange(0,width**2,1):
        if i in lefts:
            mig_matrix[i][i+1]=rate
            if i+width <= width**2-1:
                mig_matrix[i][i+width]=rate
            if i-width >= 0:
                mig_matrix[i][i-width]=rate
            if diagonal:
                if i-(width-1) >= 0:
                    mig_matrix[i][i-(width-1)]=np.sqrt(2)*rate
                if i+(width+1) <= width**2-1:
                    mig_matrix[i][i+(width+1)]=np.sqrt(2)*rate
        if i in rights:
            mig_matrix[i][i-1]=rate
            if i+width <= width**2-1:
                mig_matrix[i][i+width]=rate
            if i-width >= 0:
                mig_matrix[i][i-width]=rate
            if diagonal:
                if i-(width+1) >= 0:
                    mig_matrix[i][i-(width+1)]=np.sqrt(2)*rate
                if i+(width-1) <= width**2-1:
                    mig_matrix[i][i+(width-1)]=np.sqrt(2)*rate
        if i not in lefts and i not in rights:
            if i+1 <= width**2-1:
                mig_matrix[i][i+1]=rate
            if i-1 >= 0:
                mig_matrix[i][i-1]=rate
            if i+width <= width**2-1:
                mig_matrix[i][i+width]=rate
            if i-width >= 0:
                mig_matrix[i][i-width]=rate
            if diagonal:
                if i-(width+1) >= 0:
                    mig_matrix[i][i-(width+1)]=np.sqrt(2)*rate
                if i+(width-1) <= width**2-1:
                    mig_matrix[i][i+(width-1)]=np.sqrt(2)*rate
                if i-(width-1) >= 0:
                    mig_matrix[i][i-(width-1)]=np.sqrt(2)*rate
                if i+(width+1) <= width**2-1:
                    mig_matrix[i][i+(width+1)]=np.sqrt(2)*rate
    return mig_matrix

##migration matrix for toroidal landscapes (no edges)
def get_torus_migration_matrix(npops,rate,diagonal):
    npops=26
    if not np.sqrt(npops)%1 == 0:
        print("error - npops should make a square")
        return
    width=round(np.sqrt(npops))
    mig_matrix=np.zeros((width*width,width*width))
    mig_matrix[i][i+1]=rate
    mig_matrix[i][i-1]=rate
    mig_matrix[i][i+width]=rate
    mig_matrix[i][i-width]=rate
    if diagonal:
        mig_matrix[i][i-(width+1)]=np.sqrt(2)*rate
        mig_matrix[i][i+(width-1)]=np.sqrt(2)*rate
        mig_matrix[i][i-(width-1)]=np.sqrt(2)*rate
        mig_matrix[i][i+(width+1)]=np.sqrt(2)*rate
    return mig_matrix

#returns coordinates per node (haploid chromosome)
def get_coords_per_sample(width,haploids_per_pop):
    samples=np.arange(0,haploids_per_pop*width**2,1)
    locs=[]
    x=0
    y=0
    for x in range(width):
        for y in range(width):
            for z in range(int(haploids_per_pop)):
                locs.append([x,y])
    return np.array(locs)

#########################set params#################################
#simulation params equivalent to W50 slim runs w/250-ind demes
width=50
haploids_per_pop=10
joint_N=2.5e6
#nm=np.random.uniform(0,100)
nm=100
rate=nm/(joint_N/width**2) #numerator is n migrants per generation

mig_matrix=get_grid_migration_matrix(width,rate,diagonal=True)

#generate population configurations and demographic events
pop_configs=[]
for i in range(width**2):
    pop_configs.append(msp.PopulationConfiguration(sample_size=haploids_per_pop,
                                                   initial_size=joint_N/width**2))

################### run it ##################
trees=msp.simulate(population_configurations=pop_configs,
                   migration_matrix=mig_matrix,
                   mutation_rate=0,
                   recombination_rate=1e-8,
                   length=1e8)

#sample 60 random individuals (where individuals are )
#nodes_to_sample=np.random.choice(np.arange(0,width**2*haploids_per_pop,2),60,replace=False)
#nodes_to_sample=[[x,x+1] for x in nodes_to_sample]
#nodes_to_sample=[x for y in nodes_to_sample for x in y]
#nodes_to_sample=np.array(sorted(nodes_to_sample))
#trees=trees.simplify(nodes_to_sample.astype(np.int32))

trees.dump("/projects/kernlab/cbattey2/spaceness/sims/msp_grid_large/m_"+str(rate)+"_sim_"+str(int(np.random.uniform(1e8,1e9-1,1)))+".trees")


##check IBD correlations
# locs=get_coords_per_sample(5,10)
# haps=trees.genotype_matrix()
# genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
# allele_counts=genotypes.count_alleles()
# genotype_allele_counts=genotypes.to_allele_counts()
# positions=np.array([s.position for s in trees.sites()])
#
# #genetic & spatial distance matrices
# gen_dist=allel.pairwise_dxy(pos=positions,
#                             gac=genotype_allele_counts,
#                             start=0,stop=1e6)
# sp_dist=np.array(scipy.spatial.distance.pdist(locs[np.arange(0,len(locs),2)]))
# np.corrcoef(gen_dist,sp_dist)
# plt.scatter(sp_dist,gen_dist,edgecolor="black",facecolor="None")
# plt.axis([np.min(sp_dist)-np.std(sp_dist),np.max(sp_dist)+np.std(sp_dist),
#           np.min(gen_dist)-np.std(gen_dist),np.max(gen_dist)+np.std(gen_dist)])
