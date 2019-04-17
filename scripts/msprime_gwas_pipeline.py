#msprime gwas pipeline
import msprime, numpy as np, os, re

#neutral coalescent simulation
ts=msprime.simulate(Ne=1e5,sample_size=100,length=1e7,recombination_rate=1e-8,mutation_rate=1e-8)
haps=ts.genotype_matrix()
positions=np.array([s.position for s in ts.sites()])

#choose 100 random SNPs with allele freqs bw 0.1 and 0.9
snp_indices=[]
for i in range(100):
  phen_snp_af=0
  while phen_snp_af<0.1 or phen_snp_af>0.9:
    rsnp=np.random.choice(np.arange(0,np.shape(haps)[0],1))
    phen_snp_af=sum(haps[rsnp,:])/len(haps[rsnp,:])
  snp_indices.append(rsnp)
  
#generate phenotypes by drawing from a normal distribution with sd = 1 and mean proportional to n "phenotype" snps
phenotypes=[]
i=1
for i in range(100):
  snp_count=sum(haps[snp_indices,i])
  phenotypes.append(np.random.normal(snp_count*0.02,1,1)[0]) 

#outputs to use for GWAS-style inference
out=[haps,positions,phenotypes]


