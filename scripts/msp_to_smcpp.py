#run an msprime simulation and analyze the genotypes in smcpp
import msprime, sys, os, numpy as np, subprocess as sp
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("i")
args=parser.parse_args()

#args=argparse.Namespace(i=1)
os.chdir("/home/cbattey2/spaceness/demography/smcpp/out_msp_growth")
#os.chdir("/users/cj/spaceness/demography/smcpp/out_msp")


sp.run(["mkdir",str(args.i)])
pop_config=msprime.PopulationConfiguration(sample_size=40,initial_size=6.1e5,growth_rate=1e-3)
tree=msprime.simulate(population_configurations=[pop_config],length=1e8,recombination_rate=1e-9,mutation_rate=1e-8)

tree.write_vcf(open(str(args.i)+"/"+str(args.i)+".vcf","w"),2)
sp.check_output(["bgzip",str(args.i)+"/"+str(args.i)+".vcf"])
sp.check_output(["tabix",str(args.i)+"/"+str(args.i)+".vcf"+".gz"])
inds=",".join(["msp_"+str(x) for x in range(20)])
for j in range(1):
    command=["smc++",
             "vcf2smc",
             #"-d","msp_"+str(int(np.random.uniform(0,19,1))),"msp_"+str(int(np.random.uniform(0,19,1))),
             "--length","100000000",
             "--cores","1",
             str(args.i)+"/"+str(args.i)+".vcf"+".gz",
             str(args.i)+"/"+str(args.i)+"_set"+str(j)+".smcpp",
             "1",
             "pop1:"+inds]
    sp.run(command)
command=("smc++ estimate --cores 1 -o "+
         str(args.i)+" -p 0 -r 1e-9 --unfold 1e-8 "+
         str(args.i)+"/"+str(args.i)+"_set"+"*")
sp.run(command,shell=True)

sp.run(["smc++","plot","--csv",
        str(args.i)+".pdf",
        str(args.i)+"/"+"model.final.json"])
