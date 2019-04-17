"""
Utilities for working with scm++
"""
import logging
import subprocess
import tskit
import msprime
import os
from slimtools import *
def write_smcpp_file(path):
    """
    Writes a smcpp input file given a treesequence
    """
    ts = tskit.load(path)
    # write a vcf intermediate input
    with open(path+".vcf", "w") as vcf_file:
        ts.write_vcf(vcf_file, 2)
    # index the vcf
    cmd = f"bgzip {path}.vcf"
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    vz_file = f"{path}.vcf.gz"
    cmd = f"tabix {vz_file}"
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    # run smc++ for vcf conversion
    smc_file = f"{path}.smc.gz"
    cmd = f"smc++ vcf2smc {vz_file} {smc_file} 1 pop1:"
    for n in range(ts.num_samples // 2):
        cmd = cmd + f"msp_{n},"
    cmd = cmd[0:-1]
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)


def run_smcpp_estimate(input_file, mutation_rate, ncores):
    """
    Runs smc++ estimate on the specified file, resulting in the output being written
    to the file input_file.final.jason".
    """
    cmd = (
        f"smc++ estimate "
        f"{mutation_rate} {input_file} --base {input_file} --cores {ncores}")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)


def run_smcpp_plot(input_file, generation_time):
    """
    Runs smc++ plot on the specified file, resulting in the output being written
    to the file input_file.png".
    """
    cmd = (
        f"smc++ plot {input_file}.png {input_file} -g {generation_time} -c")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)

simdir="/home/cbattey2/spaceness/sims/slimout/spatial/W50_run3"
outdir="/home/cbattey2/spaceness/demography/smcpp/random_sampling/"

os.chdir(outdir)
sims=os.listdir(simdir)
#f=sims[1]
for f in sims:
    simname=os.path.basename(f)
    label=float(re.split("_",simname)[1])
    if label < 4:
        ts=sample_treeseq(infile=os.path.join(simdir,f),
                          outfile="",
                          nSamples=20,
                          sampling="random",
                          recapitate=False,
                          recombination_rate=1e-9,
                          write_to_file=False,
                          sampling_locs=[[12.5,12.5],[12.5,37.5],[37.5,37.5],[37.5,12.5]],
                          plot=False,
                          seed=12345)
        #get generation times estimated from short all-individual simulations
        gentimes=np.loadtxt("/home/cbattey2/spaceness/W50sp_gentimes.txt")
        gentime=[x[0] for x in gentimes if np.round(x[1],5)==np.round(label,5)]
        ts=msp.mutate(ts,1e-8/gentime[0],random_seed=12345)
        ts.dump(os.path.join(outdir,str(label)+".trees"))
        write_smcpp_file(str(label)+".trees")
        run_smcpp_estimate(str(label)+".trees"+".smc.gz",1e-8,50)
        run_smcpp_plot(str(label)+".trees"+".smc.gz"+".final.json",1)
