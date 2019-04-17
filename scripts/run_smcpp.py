#run smcpp analysis on an input tree sequence
import msprime,pyslim,os,itertools,re,glob
import numpy as np
import subprocess as sp
from slimtools import *
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("--infile",dest="infile")
parser.add_argument("--gentimes",dest="gentimes")
parser.add_argument("--outdir",dest="outdir")
args=parser.parse_args()

#debug
# args=argparse.Namespace(infile = "/Users/cj/spaceness/sims/slimout/spatial/W50_run3/sigma_0.3332292532035142_.trees_2054475",
#                         outdir = "/Users/cj/Desktop/")

#temporary folder for intermediate output
simname=os.path.basename(args.infile)
label=float(re.split("_",simname)[1])
id=re.split("_",simname)[3]
tmpout=os.path.join(args.outdir,id)
sp.run(["mkdir",tmpout])

#load treeseq, sample, mutate, export VCF, compress, index
ts=pyslim.load(args.infile)
ts=sample_treeseq(infile=ts,
                  outfile="",
                  nSamples=20,
                  sampling="random",
                  recapitate=False,
                  recombination_rate=1e-9,
                  write_to_file=False,
                  sampling_locs="",
                  plot=False,
                  seed=12345)
gentimes=np.loadtxt(args.gentimes)
gentime=[x[0] for x in gentimes if np.round(x[1],5)==np.round(label,5)]

ts=msp.mutate(ts,args.mu/gentime[0],random_seed=98675)

ts.write_vcf(open(os.path.join(tmpout,simname)+".vcf","w"),2)
sp.check_output(["bgzip",os.path.join(tmpout,simname)+".vcf"])
sp.check_output(["tabix",os.path.join(tmpout,simname)+".vcf"+".gz"])

#convert VCF to smcpp allele counts format, looping over randomly selected "distinguished individuals"
inds=",".join(["msp_"+str(x) for x in range(20)])
#for i in range(5):
command=["smc++",
         "vcf2smc",
         #"-d","msp_"+str(int(np.random.uniform(0,19,1))),"msp_"+str(int(np.random.uniform(0,19,1))),
         "--length","100000000",
         "--cores","2",
         os.path.join(tmpout,simname)+".vcf"+".gz",
         os.path.join(tmpout,simname)+"_set"+str(i)+".smcpp",
         "1",
         "pop1:"+inds]
sp.run(command)

#fit model over all DI subsets
command=("smc++ estimate --cores 2 -o "+
         tmpout+" -p 0 -r 1e-9 1e-8 "+
         os.path.join(tmpout,simname)+"_set"+"*")
sp.run(command,shell=True)
sp.run(["mv",
        os.path.join(tmpout,"model.final.json"),
        os.path.join(args.outdir,simname)+".model.json"])
sp.run(["smc++","plot","--csv",
        os.path.join(args.outdir,simname)+".pdf",
        os.path.join(args.outdir,simname)+".model.json"])

#clean up
sp.run(["rm","-r",tmpout])


#parallel -j 30 python scripts/run_smcpp.py --infile {1} --gentimes ~/spaceness/W50sp_gentimes.txt --outdir ~/spaceness/demography/smcpp/out ::: `ls ~/spaceness/sims/slimout/spatial/W50_run3/`
