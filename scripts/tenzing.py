#tenzing: a deep learning approach to genome-wide association studies
from slimtools import *
import numpy as np, subprocess as sp, msprime as msp
import pyslim, os
from tqdm import tqdm
import keras
os.chdir("/Users/cj/spaceness/")

#sample and mutates treeseqs
simdir="sims/slimout/random_mating/W35"
n_muts_per_sim=10
nSamples=100
trees=[f for f in os.listdir(simdir) if not f.startswith(".")]
haplotypes=[];phenotypes=[];phen_snps=[]
for tfile in tqdm(trees[0:10]): #drop the indices to run the full set
    for rep in range(n_muts_per_sim):
        #sample and mutate trees
        ts=sample_treeseq(infile=os.path.join(simdir,tfile),
                          outfile="",
                          nSamples=nSamples,
                          recapitate=True,
                          recombination_rate=0.25e-8,
                          write_to_file=False,
                          seed=12345)
        ts=msp.mutate(ts,0.25e-8,random_seed=12345)
        ind_haps=ts.genotype_matrix()
        haplotypes.append(ind_haps)

        #create phenotypes from a single SNP
        genotype_counts=allel.HaplotypeArray(ind_haps).to_genotypes(ploidy=2).to_allele_counts()
        sim_phen=[]
        phen_snp_af=0
        while (phen_snp_af<0.1 or phen_snp_af>0.2):
            phen_snp=np.random.choice([i for i in range(np.shape(ind_haps)[0])],1)
            phen_snp_af=sum(genotype_counts[phen_snp,:,1][0])/(nSamples*2)
        for i in range(nSamples):
            if genotype_counts[phen_snp,i,1]==0:
                sim_phen.append(np.repeat(np.random.normal(100,10,1)[0],2))
            if  genotype_counts[phen_snp,i,1]==1:
                sim_phen.append(np.repeat(np.random.normal(120,10,1)[0],2)) #start w a dominant snp (or use genotype_counts instead of haplotype?)
            if  genotype_counts[phen_snp,i,1]==2:
                sim_phen.append(np.repeat(np.random.normal(120,10,1)[0],2))
        phenotypes.append(sim_phen)
        phen_snps.append(phen_snp)

#haplotype matrix plot
from matplotlib import pyplot as plt
plt.imshow(np.transpose(haplotypes[0][0:1000]))
plt.show()

#pad haplotype matrices with -1 to keep dimensions constant
maxlen=max([np.shape(x)[0] for x in haplotypes])
padded_haps=keras.preprocessing.sequence.pad_sequences(haplotypes,maxlen,int,'post',value=-1)

#one-hot encode position of causal SNPs
labels=[np.zeros(maxlen) for i in range(np.shape(padded_haps)[0])]
for i in range(len(labels)):
    labels[i][phen_snps[i]]=1
labels=np.array(labels)

#split testing and training sets (why is this so slow?)
test_indices=np.random.choice(np.arange(0,len(padded_haps),1),10,replace=False)
train_indices=np.setdiff1d(np.arange(0,len(padded_haps),1),test_indices)
test_haps=padded_haps[test_indices]
test_labels=labels[test_indices]
train_haps=padded_haps[train_indices]
train_labels=labels[train_indices]

#CNN
def CNN1D(inputShape):

    img_1_inputs = Input(shape=(inputShape[0][1],inputShape[0][2]))

    h = layers.Conv1D(1250, kernel_size=2, activation='relu', name='conv1_1')(img_1_inputs)
    h = layers.Conv1D(512, kernel_size=2, dilation_rate=1, activation='relu')(h)
    h = layers.AveragePooling1D(pool_size=2)(h)
    h = layers.Dropout(0.25)(h)
    h = layers.Conv1D(512, kernel_size=2, activation='relu')(h)
    h = layers.AveragePooling1D(pool_size=2)(h)
    h = layers.Dropout(0.25)(h)
    h = layers.Flatten()(h)

    loc_input = Input(shape=(inputShape[1][1],))
    m2 = layers.Dense(64,name="m2_dense1")(loc_input)
    m2 = layers.Dropout(0.1)(m2)

    h =  layers.concatenate([h,m2])
    h = layers.Dense(128,activation='relu')(h)
    h = Dropout(0.2)(h)
    output = layers.Dense(1,kernel_initializer='normal',name="out_dense",activation='linear')(h)

    model = Model(inputs=[img_1_inputs,loc_input], outputs=[output])
    model.compile(loss='mse', optimizer='adam')
    model.summary()

    return model
