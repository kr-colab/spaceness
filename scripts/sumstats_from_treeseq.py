from slimtools import *
import argparse
#os.chdir("/Users/cj/spaceness")

parser = argparse.ArgumentParser(
            description='Calculate summary statistics from a  \
                         tree sequence produced by msprime or SLiM')
parser.add_argument('--infile', dest='infile',
                    help='full path to tree sequence')
parser.add_argument('--outfile', dest='outfile',
                   help='full path for output summary statistics.')
parser.add_argument('--seed', dest='seed',type=int,
                   help='Random seed for sampling individuals and adding mutations to tree sequences.')
parser.add_argument('--sampling',dest='sampling',
                   help='Type of sampling to use. Options: \
                   random (nSamples random individuals). \
                   point (nSamples individuals closest to sampling_location). \
                   pair (nSamples individuals closest to sampling_location \
                   and its reflection over the x axis)')
parser.add_argument('--mu', dest='mu',
                   help='mutation rate per base per unit time in SLiM simulations.')
parser.add_argument('--sampling_location', dest='sampling_location',
                   help='Sampling site location given as "x,y"')
args=parser.parse_args()

#params for debugging
args=argparse.Namespace(infile='/Users/cj/spaceness/sims/slimout/spatial/W35/sigma_0.2130878654199541_.trees',
                        outfile='/Users/cj/Desktop/test_sumstats.txt',
                        sampling='random',
                        mu=0.25e-8,
                        seed=12345)

print("sampling individuals and adding mutations")
ts=sample_treeseq(infile=args.infile,
                  outfile="",
                  nSamples=100,
                  recapitate=False,
                  recombination_rate=1e-8,
                  write_to_file=False,
                  seed=args.seed)
ts=msp.mutate(ts,args.mu,random_seed=args.seed)

print("loading haplotype matrices")
haps,pos,locs=get_ms_outs(ts)
label=float(re.sub("sigma_|_.trees*|.trees","",os.path.basename(args.infile)))
pos=discretize_snp_positions(pos)

print("calculating summary statistics")
ss=getSLiMSumStats(haps,
                   pos,
                   label,
                   locs,
                   args.outfile,
                   maxlen=1e8,
                   ibs_tracts=True,
                   verbose=True)
