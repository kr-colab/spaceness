#msprime island model grid simulation
import numpy as np, msprime as msp, os, allel, itertools, scipy, math
from scipy.spatial.distance import pdist, squareform
import seaborn
from matplotlib import pyplot as plt
import matplotlib
from sklearn.linear_model import LinearRegression
from slimtools import *

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

##get parameters in grid format (xloc,yloc,migration_right,migration_left,migration_up,migration_down)
def get_grid_params(width,rate):
    popy=np.repeat(np.arange(0,width,1),width)
    popx=[np.arange(0,width,1) for x in range(width)] #man... fuck python
    popx=[x for y in popx for x in y]
    mr=np.repeat(rate,width**2)
    ml=np.repeat(rate,width**2)
    md=np.repeat(rate,width**2)
    mu=np.repeat(rate,width**2)
    return(np.transpose([popx,popy,ml,mr,md,mu]))
#grid_params=get_grid_params(10,.1/(1e5/100))

##convert grid parameters to migration matrix.
def grid_params_to_migration_matrix(grid_params):
    npops=np.shape(grid_params)[0]
    width=math.sqrt(npops)
    lefts=np.arange(0,width**2,width)
    rights=np.arange(width-1,width**2,width)
    migmat=np.zeros((npops,npops))
    for i in range(len(grid_params)):
        row=grid_params[i]
        if i not in rights: migmat[i,i+1]=row[2]
        if i not in lefts: migmat[i,i-1]=row[3]
        if i+width <=(npops-1): migmat[i,int(i+width)]=row[4]
        if i-width > -1: migmat[i,int(i-width)]=row[5]
    return(migmat)
#mig_matrix=grid_params_to_migration_matrix(grid_params)

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

#########################running a simulation#################################
#simulation params equivalent to W50 slim runs w/250-ind demes
width=zz
haploids_per_pop=20
joint_N=100000
rate=10/(joint_N/width**2) #numerator is n migrants per generation
barrier_m_scaling=0.01 #scaling of n migrants relative to the base rate across "barrier" regions
tcoal=10*joint_N

#for k in range(100):
#generate grid params for migration matrix
grid_params=get_grid_params(width,rate)
h_barrier_loc=np.random.choice(np.arange(1,width-1),1)[0]
h_barrier_rate=0.1/(joint_N/width**2)
#h_barrier_rate=np.max([0,np.random.normal(rate*barrier_m_scaling,rate*barrier_m_scaling,1)])
#v_barrier_loc=np.random.choice(np.arange(1,width-1),1)
#v_barrier_rate=np.max([0,np.random.normal(rate*barrier_m_scaling,rate*barrier_m_scaling,1)])
grid_params[grid_params[:,1]==h_barrier_loc,5]=h_barrier_rate
grid_params[grid_params[:,1]==h_barrier_loc-1,4]=h_barrier_rate
#grid_params[grid_params[:,0]==v_barrier_loc,3]=v_barrier_rate
#grid_params[grid_params[:,0]==v_barrier_loc-1,2]=v_barrier_rate
mig_matrix=grid_params_to_migration_matrix(grid_params)

paramsfile=open("/Users/cj/Dropbox/pytheas/grid_sim_params.txt","a")
paramsfile.write(str(k)+" "+str(h_barrier_loc)+" "+str(h_barrier_rate))
paramsfile.close()

#mig_matrix=np.loadtxt("/Users/cj/Desktop/mig_matrix_var.txt")
#generate population configurations and demographic events
nodelocs=get_coords_per_sample(width,haploids_per_pop)
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
                   recombination_rate=1e-9,
                   length=1e6)
trees.dump("/Users/cj/Dropbox/pytheas/grid_sims/grid_sim"+str(k)+".trees")

#trees=msp.mutate(trees,1e-8)

#######################analyze simulation output###########################
#sample 100 individuals
nodes_to_sample=np.random.choice(np.arange(0,width**2*haploids_per_pop,2),50,replace=False)
nodes_to_sample=[[x,x+1] for x in nodes_to_sample]
nodes_to_sample=[x for y in nodes_to_sample for x in y]
nodes_to_sample=np.array(sorted(nodes_to_sample))
locs=nodelocs[nodes_to_sample]
ts=trees.simplify(nodes_to_sample.astype(np.int32))

#load genotypes for analysis
haps=ts.genotype_matrix()
genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
allele_counts=genotypes.count_alleles()
genotype_allele_counts=genotypes.to_allele_counts()
positions=np.array([s.position for s in ts.sites()])

#genetic & spatial distance matrices
gen_dist=allel.pairwise_dxy(pos=positions,
                            gac=genotype_allele_counts,
                            start=0,stop=1e6)
sp_dist=np.array(scipy.spatial.distance.pdist(locs[np.arange(0,len(locs),2)]))

#IBD regression
lm=LinearRegression()
lm=lm.fit(sp_dist.reshape(-1,1),gen_dist.reshape(-1,1))
preds=lm.predict(sp_dist.reshape(-1,1))
residuals=[gen_dist[x]-preds[x][0] for x in range(len(preds))]

######################### plot ########################
fig = plt.figure(figsize=(6,2.5),dpi=400)
matplotlib.rcParams.update({'font.size': 6})

#pop configuration
ax=fig.add_axes([0,.1,.3,.8])
ax.scatter(x=grid_params[:,0],
        y=grid_params[:,1],
        edgecolor="black",facecolor="white",
        s=50,clip_on=False,linewidth=0.5)
ax.scatter(x=np.array(locs)[:,0],y=np.array(locs)[:,1],
           clip_on=False,s=35,facecolor="grey")
for i in range(len(grid_params)):
    ax.arrow(grid_params[i,0]+.4,grid_params[i,1],0.23,0,head_width=0.12,width=0.015,head_length=0.13,
             color=plt.cm.RdYlBu(1-grid_params[i,2]/np.max(grid_params[:,2])))
    ax.arrow(grid_params[i,0]-.4,grid_params[i,1],-0.23,0,head_width=0.12,width=0.015,head_length=0.13,
             color=plt.cm.RdYlBu(1-grid_params[i,3]/np.max(grid_params[:,3])))
    ax.arrow(grid_params[i,0],grid_params[i,1]+.4,0,0.23,head_width=0.12,width=0.015,head_length=0.13,
             color=plt.cm.RdYlBu(1-grid_params[i,4]/np.max(grid_params[:,4])))
    ax.arrow(grid_params[i,0],grid_params[i,1]-.4,0,-0.23,head_width=0.12,width=0.015,head_length=0.13,
             color=plt.cm.RdYlBu(1-grid_params[i,5]/np.max(grid_params[:,5])))
ax.axis([-0.1,width-1+0.1,-0.1,width-1+0.1])
ax.axis('off')

#colorbar
cbaxes=fig.add_axes([0.32,0.125,.02,.65])
cmap=matplotlib.cm.RdYlBu_r
norm=matplotlib.colors.Normalize(vmin=barrier_rate,vmax=rate)
cb1=matplotlib.colorbar.ColorbarBase(cbaxes,cmap=cmap,norm=norm,orientation='vertical',format="%.0e")
cb1.set_label('Migration\nRate',labelpad=-20,y=1.15,rotation=0,multialignment='left',)

#IBD scatter plot
ax2=fig.add_axes([0.5,0.625,0.4,0.3])
ax2.scatter(sp_dist,gen_dist,edgecolors='k',facecolors='',linewidths=0.5,clip_on=False,s=15)
ax2.plot(sp_dist,preds,color="red")
ax2.axis([np.min(sp_dist)-.1*np.std(sp_dist),np.max(sp_dist)+.1*np.std(sp_dist),
          np.min(gen_dist)-.5*np.std(gen_dist),np.max(gen_dist)+.5*np.std(gen_dist)])
ax2.set_xlabel('Spatial Distance')
ax2.set_ylabel('Genetic Distance')

#pairwise distance matrices
ax3=fig.add_axes([0.42,0.1,0.25,0.375])
mask =  np.tri(np.shape(squareform(gen_dist))[0], k=-1)
gen_dist_plot = np.ma.array(squareform(gen_dist), mask=mask==0)
sp_dist_plot = np.ma.array(squareform(sp_dist), mask=mask)
ax3.imshow(sp_dist_plot,cmap=matplotlib.cm.YlGnBu)
ax3.imshow(gen_dist_plot,cmap=matplotlib.cm.YlOrBr)
ax3.axis('off')

cbax2=fig.add_axes([0.65,0.1,0.01,0.3])
cb2=matplotlib.colorbar.ColorbarBase(cbax2,
                                     cmap=matplotlib.cm.YlGnBu,
                                     norm=matplotlib.colors.Normalize(vmin=np.min(sp_dist),vmax=np.max(sp_dist)),
                                     orientation='vertical',format="%.0e")
cbax3=fig.add_axes([0.425,0.05,0.225,0.02])
cb3=matplotlib.colorbar.ColorbarBase(cbax3,
                                     cmap=matplotlib.cm.YlOrBr,
                                     norm=matplotlib.colors.Normalize(vmin=np.min(gen_dist),vmax=np.max(gen_dist)),
                                     orientation='horizontal',format="%.0e")

ax4=fig.add_axes([0.75,0.1,0.2,0.4])
ax4.imshow(squareform(residuals))

fig.savefig("/Users/cj/Desktop/testpopconfig.pdf",bbox_inches='tight')

h_barrier_rate
v_barrier_rate
