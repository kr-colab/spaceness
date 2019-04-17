from slimtools import *
import msprime as msp,pyslim,numpy as np,os,subprocess as sp
import argparse
import allel


files=os.listdir("/projects/kernlab/cbattey2/spaceness/gwas/sims/")

for i in files:
    print("processing file"+str(i))
    args = argparse.Namespace(treeseq=os.path.join("/projects/kernlab/cbattey2/spaceness/gwas/sims/",i),
                              outdir="/projects/kernlab/cbattey2/spaceness/gwas/out/",
                              plink_path="/projects/kernlab/cbattey2/plink-1.07-x86_64/plink",
                              vcftools_path="/projects/kernlab/cbattey2/vcftools_0.1.13/bin/vcftools",
                              nSamples=1000,
                              mu=0.25e-8,
                              phenotype="transform_coord")

    #sample individuals and add mutations
    simname=os.path.basename(args.treeseq)
    ts=pyslim.load(args.treeseq)
    sample_inds=np.unique([ts.node(j).individual for j in ts.samples()]) #final-generation individuals
    subsample=np.random.choice(sample_inds,args.nSamples,replace=False) #get nSamples random sample inds
    subsample_nodes=[ts.individual(x).nodes for x in subsample] #node numbers for sampled individuals
    subsample_nodes=[a for b in subsample_nodes for a in b] #flatten the list
    subsample_nodes=np.sort(np.array(subsample_nodes))
    ts=ts.simplify(subsample_nodes)
    ts=msp.mutate(ts,args.mu)

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

    #convert vcf to .ped (throwing error for opening temp files when run from command line on mac... switch to manual ped file creation?)
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
            ind_phenotype=np.random.normal(locs[i][0]/2.5+args.phenotype_mean,args.phenotype_sd,1)[0] #if W=50, means will range from 160 - 180 to mimic a cline in cm height
            phenfile.write("msp_"+str(i)+" "+"msp_"+str(i)+" "+str(ind_phenotype)+"\n")
        phenfile.close()

    if args.phenotype=="corner_bimodal":
        phenfile=open(os.path.join(args.outdir,simname)+".phenotypes","w")
        for i in range(args.nSamples):
            if locs[i][0]<20 and locs[i][1]<20:
                ind_phenotype=np.random.normal(args.phenotype_mean,args.phenotype_sd,1)[0]
            else:
                ind_phenotype=np.random.normal(args.phenotype_mean+3*args.phenotype_sd,args.phenotype_sd,1)[0]
            phenfile.write("msp_"+str(i)+" "+"msp_"+str(i)+" "+str(ind_phenotype)+"\n")
        phenfile.close()


    #run plink association analysis with PC coords as covariates
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

    #run plink without covariates
    sp.check_output([args.plink_path,
                     "--noweb",
                     "--file",
                     os.path.join(args.outdir,simname),
                     "--pheno",
                     os.path.join(args.outdir,simname)+".phenotypes",
                     "--allow-no-sex",
                     "--out",
                     os.path.join(args.outdir,simname),
                     "--assoc"])

    #clean up big files
    sp.check_output(["rm",
                     os.path.join(args.outdir,simname)+".ped",
                     os.path.join(args.outdir,simname)+".map",
                     os.path.join(args.outdir,simname)+".nosex",
                     ])
