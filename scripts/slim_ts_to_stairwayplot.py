import msprime as msp, numpy as np, allel, argparse, subprocess,os,re
from os.path import basename, splitext, join
from slimtools import *

parser = argparse.ArgumentParser(
  description='Run stairwayplot on a tree sequence output from msprime.\
               Requires gnu-parallel and Stairwayplot including \
               the command-line version of Stairway_plot_output_summary,\
               which can be downloaded from: \
               https://sites.google.com/site/jpopgen/stairway-plot.')
parser.add_argument('--infile',dest='infile',
                    help='path to tree sequence')
parser.add_argument('--sampling',dest='sampling',
                    help='sampling scheme (random, midpoint, or point).')
parser.add_argument('--outdir',dest='outdir',
                    help='full path to directory for writing output files')
parser.add_argument('--nboots',dest='nboots',type=int,
                    help='number of bootstrap replicates to run.')
parser.add_argument('--stairwayplot_dir',dest='stairwayplot_dir',
                    help="path to stairwayplot directory")
parser.add_argument('--cores',dest='cores',
                    help="number of cores to use when fitting demographic \
                    models.")
parser.add_argument('--seed',dest='seed',type=int,
                    help="random number seed for bootstrap resampling.")
parser.add_argument('--mutation_rate',dest='mutation_rate',type=float,
                    help="mutation rate (per base per generation)")
parser.add_argument('--gentimes_file',dest='gentimes_file',
                    help='path to file listing sigmas and generation times.')
parser.add_argument('--generation_time',dest='generation_time',type=float,
                    help="generation time for plotting")
parser.add_argument('--plot',dest='plot',type=bool,
                    help='plot Ne ~ t ? T/F')
args=parser.parse_args()

np.random.seed(args.seed)

#read in and sample treeseq
simname=os.path.basename(args.infile)
label=float(re.split("_",simname)[1])
ts=sample_treeseq(infile=args.infile,
                  outfile="",
                  nSamples=20,
                  sampling=args.sampling,
                  recapitate=False,
                  recombination_rate=1e-9,
                  write_to_file=False,
                  sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
                  plot=False,
                  seed=args.seed)

#get generation times estimated from short all-individual simulations
gentimes=np.loadtxt(args.gentimes_file)
gentime=[x[0] for x in gentimes if np.round(x[1],5)==np.round(label,5)]

ts=msp.mutate(ts,args.mutation_rate/gentime[0],random_seed=args.seed)

#count alleles, bootstrap over sites, return the SFS minus the 0% bin
haps=np.array(ts.genotype_matrix())
genotypes=allel.HaplotypeArray(haps).to_genotypes(ploidy=2)
allele_counts=genotypes.count_alleles()
sfs=allel.sfs(allele_counts[:,1])
sfs=sfs[1:len(sfs)]

#tmp directory for input files
command=("cd "+args.outdir+";"+
         "mkdir infiles")
subprocess.run(command,shell=True)

#write stairwayplot input
out=open((join(args.outdir,
               "infiles",
               splitext(basename(args.infile))[0])+
         "_strwyplt.txt"),
         "w")
out.write(("slim"+"\t"+
           str(ts.num_samples)+"\t"+
           str(int(ts.sequence_length))+"\t"+
           str(1)+"\t"+
           str(ts.num_samples-1)+"\n")) #order is name,n_samples,sequence_length,lowest_sfs_bin,highest_sfs_bin
for x in sfs:
    out.write(str(int(x))+"\t")
out.close()

#write bootstrapped inputs
for i in range(args.nboots):
    nsites=np.shape(allele_counts)[0]
    bootset=np.random.choice(np.arange(0,nsites,1),nsites,replace=True)
    bootac=allele_counts[bootset,:]
    bootsfs=allel.sfs(bootac[:,1])
    bootsfs=bootsfs[1:len(bootsfs)]
    out=open((join(args.outdir,
              "infiles",
              splitext(basename(args.infile))[0])+
            "_strwyplt"+str(i)+".txt"),"w")
    out.write(("msp"+"\t"+
               str(ts.num_samples)+"\t"+
               str(int(ts.sequence_length))+"\t"+
               str(1)+"\t"+
               str(ts.num_samples-1)+"\n"))
    for x in bootsfs:
        out.write(str(int(x))+"\t")
    out.close()

#fit models to bootstraps in parallel
command=("cd "+args.stairwayplot_dir+";"+
         "files="+join(args.outdir,"infiles")+"/*;"+
         "parallel -j "+str(args.cores)+" "+
         "java -cp .:swarmops.jar "+
         "Stairway_plot_theta_estimation02 {1} 1 5000 "+
         "::: $files")
subprocess.run(command,shell=True)

#moving output files...
command=("cd "+args.outdir+";"+
         "mkdir thetas;"+
         "mv "+join(args.outdir,"infiles/")+"*.addTheta thetas;"+
         "rm -r infiles")
subprocess.run(command,shell=True)

#get median and 95% CI for Ne~t
command=("cd "+args.stairwayplot_dir+";"+
         "java Stairway_plot_output_summary_commandline "+
         join(args.outdir,"thetas")+" "+
         str(args.mutation_rate)+" "+
         str(args.generation_time)+" "+
         join(args.outdir,splitext(basename(args.infile))[0]+
         "_estimated_Ne.txt")+";"+
         "rm -r "+join(args.outdir,"thetas"))
subprocess.run(command,shell=True)

if args.plot:
    from matplotlib import pyplot as plt
    import pandas
    nt=pandas.read_csv(join(args.outdir,splitext(basename(args.infile))[0]+"_estimated_Ne.txt"),
              sep="\t",skiprows=5)
    nt=nt[nt['year']>10]
    f,ax = plt.subplots(figsize=(7, 7))
    ax.set(xscale="log", yscale="log",
           ylim=(np.min(nt['Ne_2.5%']-np.std(nt['Ne_2.5%'])),
                 np.max(nt['Ne_97.5%'])+np.std(nt['Ne_97.5%'])))
    ax.plot(nt['year'],nt['Ne_median'],c="red")
    ax.plot(nt['year'],nt['Ne_2.5%'],c='grey')
    ax.plot(nt['year'],nt['Ne_97.5%'],c='grey')
    f.savefig(join(args.outdir,splitext(basename(args.infile))[0]+"_Ne_est.png"),
              bbox_inches='tight')

#example run command
#python scripts/slim_ts_to_stairwayplot.py ~/spaceness/sims/slimout/spatial/W50_run3/sigma_0.21467883716606018_.trees_3812242 ~/spaceness/demography/stairwayplot/ 100 ~/stairwayplot/ 50 1234 1e-8 W50sp_gentimes.txt 1 T
files=~/spaceness/sims/slimout/spatial/W50_run3/*
for f in $files; do python scripts/slim_ts_to_stairwayplot.py --infile $f --sampling random --outdir ~/spaceness/demography/stairwayplot/random_sampling/ --nboots 200 --stairwayplot_dir ~/stairwayplot/ --cores 78 --seed 54321 --mutation_rate 1e-8 --gentimes_file W50sp_gentimes.txt --generation_time 1 --plot T; done
