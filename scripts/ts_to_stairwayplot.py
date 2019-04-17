from slimtools import *
import subprocess
import pandas
import argparse
#os.chdir("/Users/cj/spaceness")

parser = argparse.ArgumentParser(
            description='Convert slim tree sequences to stairwayplot input')
parser.add_argument('--infile', dest='infile',
                    help='full path to tree sequence')
parser.add_argument('--outfile', dest='outfile',
                   help='full path for output summary statistics.')
parser.add_argument('--seed', dest='seed',type=int,
                   help='Random seed for sampling individuals and adding mutations to tree sequences.')
parser.add_argument('--sampling',dest='sampling')
parser.add_argument('--mu', dest='mu',type=float,
                   help='mutation rate per base per unit time in SLiM simulations.')
args=parser.parse_args()

#params for debugging
# ts=msp.simulate(sample_size=120,Ne=1e4,mutation_rate=1e-8,recombination_rate=1e-8,length=1e8)
# args=argparse.Namespace(infile='/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.2357863131214634_.trees_3600449',
#                         outfile='/Users/cj/Desktop/test.txt',
#                         sampling='random',
#                         mu=1e-8,
#                         seed=12345)

#folder for input files
subprocess.run("mkdir /Users/cj/spaceness/demography/stairwayplot/input/"+re.split("\\_",os.path.basename(args.infile))[1],
               shell=True)

#load tree sequence
label=float(re.sub("sigma_|_.trees*|.trees","",os.path.basename(args.infile)))
t=pyslim.load(args.infile)

#get generation time estimated from short all-individual simulations
gentimes=np.loadtxt("/Users/cj/spaceness/W50sp_gentimes.txt")
gentime=[x[0] for x in gentimes if np.round(x[1],5)==np.round(label,5)]

print("sigma="+str(label)+"  gentime="+str(gentime[0]))

#sample and mutate treeseq
ts=sample_treeseq(infile=args.infile,
                  outfile="",
                  nSamples=60,
                  recapitate=False,
                  recombination_rate=1e-8,
                  write_to_file=False,
                  sampling=args.sampling,
                  sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
                  plot=False,
                  seed=args.seed)
ts=msp.mutate(ts,args.mu/gentime[0],random_seed=args.seed)

#get ms style outputs
print("loading haplotype matrices")
haps,pos,locs=get_ms_outs(ts)
#haps=np.array(ts.genotype_matrix()) #for testing with msprime

#write 200 bootstrapped input files by bootstrapping over sites
for i in range(100):
    out=open("/Users/cj/spaceness/demography/stairwayplot/input/"+re.split("\\_",os.path.basename(args.infile))[1]+"/"+str(i)+".txt","w")
    out.write(str(label)+"\t"+str(120)+"\t"+str(int(1e8))+"\t"+str(1)+"\t"+str(119)+"\n")
    #bootstrap over sites and return the SFS
    genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
    genotypes=genotypes[np.random.choice(np.arange(0,np.shape(genotypes)[0]-1),np.shape(genotypes)[0],replace=True),:,:]
    allele_counts=genotypes.count_alleles()
    sfs=allel.sfs(allele_counts[:,1])
    sfs=sfs[1:len(sfs)] #first element is n sites with 0 derived alleles
    for x in sfs:
        out.write(str(int(x))+"\t")
    out.close()


#prep for a folder of simulations with
#for f in ~/spaceness/sims/slimout/spatial/W50_run3/*;do python ~/spaceness/scripts/ts_to_stairwayplot.py --infile $f --outfile "test" --seed 12345 --sampling "random" --mu 1e-8;done
#or
#parallel -j 10 python ~/spaceness/scripts/ts_to_stawayplot.py --infile {1} --outfile "test" --seed 12345 --sampling "random" --mu 1e-8 ::: `ls ~/spaceness/sims/slimout/spatial/W50_run3/*`
