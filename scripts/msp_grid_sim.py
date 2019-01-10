#msprime island model grid simulation
import numpy as np, msprime as msp, os, allel, itertools,scipy
from scipy.spatial.distance import pdist, squareform
import seaborn
from matplotlib import pyplot as plt

#utility functions
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

def get_coords_per_sample(width,haploids_per_pop):
    samples=np.arange(0,haploids_per_pop*width**2,1)
    locs=[]
    x=0
    y=0
    for x in range(width):
        for y in range(width):
            for z in range(int(haploids_per_pop/2)):
                locs.append([x,y])
    return locs

#set params
width=10
rate=1e-5
joint_N=1e6
tcoal=1e6
haploids_per_pop=20

#generate migration matrix, pop configs, dem events
mig_matrix=get_grid_migration_matrix(width,rate,True)
locs=get_coords_per_sample(width,haploids_per_pop)
pop_configs=[]
for i in range(width**2):
    pop_configs.append(msp.PopulationConfiguration(sample_size=haploids_per_pop,
                                                   initial_size=joint_N/width**2))
dem_events=[msp.MigrationRateChange(time=tcoal,rate=0)]
for i in range(width**2-1):
    dem_events.append(msp.MassMigration(time=tcoal+1,
                                        source=i,
                                        destination=i+1,
                                        proportion=1))

#simulate
trees=msp.simulate(population_configurations=pop_configs,
                   demographic_events=dem_events,
                   migration_matrix=mig_matrix,
                   mutation_rate=1e-8,
                   recombination_rate=1e-8,
                   length=1e6)

#analyze output
haps=trees.genotype_matrix()
genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
allele_counts=genotypes.count_alleles()
genotype_allele_counts=genotypes.to_allele_counts()
positions=np.array([s.position for s in trees.sites()])
sp_dist=np.array(scipy.spatial.distance.pdist(locs))

#genetic distance matrix
gen_dist=allel.pairwise_dxy(pos=positions,
                            gac=genotype_allele_counts,
                            start=0,stop=1e6)

#plots
seaborn.scatterplot(x=sp_dist,y=gen_dist)
plt.imshow(squareform(gen_dist))
