from slimtools import *
import msprime as msp,pyslim,numpy as np,os,subprocess as sp
import argparse
import allel

parser = argparse.ArgumentParser(description='Run a GWAS in plink on individuals \
                                              sampled from a SLiM simulation.\
                                              Requires plink, vcftools, msprime, and pyslim.')
parser.add_argument('--treeseq', dest='treeseq',
                    help='path to simulation tree sequence.')
parser.add_argument('--outdir', dest='outdir',
                   help='Output file directory')
parser.add_argument('--plink_path', dest='plink_path',
                   help='call for plink')
parser.add_argument('--vcftools_path', dest='vcftools_path',
                   help='call for vcftools')
parser.add_argument('--nSamples', dest='nSamples', type=int,
                   help='Number of individuals to sample.')
parser.add_argument('--sampling',dest='sampling',
                   help='Type of sampling to use. Options: \
                   random (nSamples random individuals). \
                   point (nSamples individuals closest to sampling_location). \
                   pair (nSamples individuals closest to sampling_location \
                   and its reflection over the x axis)')
parser.add_argument('--mu', dest='mu',type=float,
                   help='mutation rate per base per unit time (msprime generations or SLiM time steps)')
parser.add_argument('--phenotype', dest='phenotype',
                   help='Type of phenotype to generate. \
                         options: "gaussian" (phenotypes drawn from a normal distribution); \
                                  "transform_coord" (phenotype is a transformation of the x coordinate); \
                                  "corner_bimodal" (phenotype drawn from different normal distributions \
                                  in/outside a corner of the landscape).')
parser.add_argument('--phenotype_mean', dest='phenotype_mean',
                   help='mean for gaussian phenotypes')
parser.add_argument('--phenotype_sd', dest='phenotype_sd',
                   help='standard deviation for gaussian phenotypes')
parser.add_argument('--seed', dest='seed',type=int,
                   help='Random seed for sampling individuals and adding mutations to tree sequences.')

args=parser.parse_args()

### debug
# args = argparse.Namespace(treeseq='/Users/cj/spaceness/sims/slimout/spatial/W50/coalesced/sigma_0.20156213958750416_.trees1500000.trees',
#                           outdir="/Users/cj/spaceness/gwas/out/normal_phen/",
#                           plink_path="plink",
#                           vcftools_path="vcftools",
#                           nSamples=1000,
#                           mu=0.25e-8,
#                           phenotype="transform_coord",
#                           phenotype_mean=175,
#                           phenotype_sd=8)

#sample individuals and add mutations
ts=sample_treeseq(infile=args.infile,
                  outfile="",
                  nSamples=100,
                  recapitate=False,
                  recombination_rate=1e-8,
                  write_to_file=False,
                  seed=args.seed)
ts=msp.mutate(ts,args.mu,random_seed=args.seed)

#get haplotypes and locations
haps=ts.genotype_matrix()
sample_inds=np.unique([ts.node(j).individual for j in ts.samples()])
locs=[[ts.individual(x).location[0],ts.individual(x).location[1]] for x in sample_inds]

#run a PCA
genotype_counts=allel.HaplotypeArray(haps).to_genotypes(ploidy=2).to_allele_counts()
pca=allel.pca(genotype_counts[:,:,0])
pcfile=open(os.path.join(args.outdir,simname)+".pca","w")
for i in range(args.nSamples):
    pcfile.write("msp_"+str(i)+" "+"msp_"+str(i)+" ")
    for j in range(10):
        pcfile.write(str(pca[0][i][j])+" ")
    pcfile.write("\n")
pcfile.close()

#write to VCF
with open(os.path.join(args.outdir,simname)+".vcf","w") as vcf_file:
    ts.write_vcf(vcf_file,2)

#convert vcf to .ped (throwing error for opening temp files when run on command line... switch to manual ped file creation?)
sp.check_output([args.vcftools_path,
                 "--vcf",
                 os.path.join(args.outdir,simname)+".vcf",
                 "--out",
                 os.path.join(args.outdir,simname),
                 "--plink"
                ])
sp.check_output(["rm",
                 os.path.join(args.outdir,simname)+".vcf"])

#create phenotypes
if args.phenotype=="gaussian":
    phenotypes=np.random.normal(args.phenotype_mean,args.phenotype_sd,args.nSamples) #approx US height
    phenfile=open(os.path.join(args.outdir,simname)+".phenotypes","w")
    for i in range(args.nSamples):
        phenfile.write("msp_"+str(i)+" "+"msp_"+str(i)+" "+str(phenotypes[i])+"\n")
    phenfile.close()
if args.phenotype=="transform_coord":
    phenfile=open(os.path.join(args.outdir,simname)+".phenotypes","w")
    for i in range(args.nSamples):
        ind_phenotype=np.random.normal(locs[i][0]/2.5+160,7,1)[0] #if W=50, means will range from 160 - 180 to mimic a cline in cm height
        phenfile.write("msp_"+str(i)+" "+"msp_"+str(i)+" "+str(ind_phenotype)+"\n")
    phenfile.close()
if args.phenotype=="corner_bimodal":
    phenfile=open(os.path.join(args.outdir,simname)+".phenotypes","w")
    for i in range(args.nSamples):
        if locs[i][0]<25 and locs[i][1]<25:
            ind_phenotype=np.random.normal(140,10,1)[0]
        else:
            ind_phenotype=np.random.normal(110,10,1)[0]
        phenfile.write("msp_"+str(i)+" "+"msp_"+str(i)+" "+str(ind_phenotype)+"\n")
    phenfile.close()


#run plink association analysis with PC coords as covariates
if args.pca==True:
    sp.check_output([args.plink_path,
                     "--noweb",
                     "--file",
                     os.path.join(args.outdir,simname),
                     "--pheno",
                     os.path.join(args.outdir,simname)+".phenotypes",
                     "--allow-no-sex",
                     "--out",
                     os.path.join(args.outdir,simname),
                     "--linear",
                     "--covar",
                     os.path.join(args.outdir,simname)+".pca",
                     "--hide-covar"])
else:
    sp.check_output([args.plink_path,
                     "--noweb",
                     "--file",
                     os.path.join(args.outdir,simname),
                     "--pheno",
                     os.path.join(args.outdir,simname)+".phenotypes",
                     "--allow-no-sex",
                     "--out",
                     os.path.join(args.outdir,simname),
                     "--linear"])

#clean up big files
sp.check_output(["rm",
                 os.path.join(args.outdir,simname)+".ped",
                 os.path.join(args.outdir,simname)+".map",
                 os.path.join(args.outdir,simname)+".nosex",
                 ])
