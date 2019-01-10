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
                   midpoint (nSamples individuals proportional to distance from the midpoint). \
                   point (nSamples/len(sampling_locs) individuals proportional to distance from each sampling location. \
                   and its reflection over the x axis)')
parser.add_argument('--mu', dest='mu',
                   help='mutation rate per base per unit time in SLiM simulations.')
parser.add_argument('--sampling_locs', dest='sampling_locs',
                   help='Currently not active - see hard coding in sample_treeseq call. \
                   Sampling site location given as "x1,y2 x2,y2 x3,y3..."')
args=parser.parse_args()

#params for debugging
# args=argparse.Namespace(infile='/Users/cj/spaceness/sims/slimout/spatial/W35/sigma_0.2130878654199541_.trees',
#                         outfile='/Users/cj/Desktop/test_sumstats.txt',
#                         sampling='random',
#                         mu=0.25e-8,
#                         seed=12345)

print("sampling individuals and adding mutations")
#get census population size
label=float(re.sub("sigma_|_.trees*|.trees","",os.path.basename(args.infile)))
t=pyslim.load(args.infile)
nind=t.num_samples/2
out=open(args.outfile+".popsizes",'a')
out.write(str(label)+" "+str(nind)+"\n")
out.close()

#sample and mutate treeseq
ts=sample_treeseq(infile=args.infile,
               outfile="",
               nSamples=100,
               recapitate=False,
               recombination_rate=1e-8,
               write_to_file=False,
               sampling=args.sampling,
               sampling_locs=[[8,8],[27,27],[27,8],[8,27]],
               plot=False,
               seed=args.seed)

ts=msp.mutate(ts,args.mu,random_seed=args.seed)

#get ms style outputs
print("loading haplotype matrices")
haps,pos,locs=get_ms_outs(ts)
pos=discretize_snp_positions(pos)

#get summary stats
print("calculating summary statistics")
ss=getSLiMSumStats(haps,
                   pos,
                   label,
                   locs,
                   args.outfile,
                   maxlen=1e8,
                   ibs_tracts=True,
                   verbose=True,
                   min_len_to_keep=2)
