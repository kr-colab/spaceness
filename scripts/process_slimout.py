from slimtools import *
import argparse
#os.chdir("/Users/cj/spaceness")

parser = argparse.ArgumentParser(
            description='Calculate summary statistics from a set of \
                         tree sequences produced by msprime or SLiM')
parser.add_argument('--indir', dest='indir',
                    help='full path to directory with tree sequences')
parser.add_argument('--outpath', dest='outdir',
                   help='full path for output file.')
parser.add_argument('--subset', dest='subset',type=bool,
                   help='Boolean. If False, process all files in directory.\
                         If True, start and stop give the 0-based \
                         index bounds of files to process as returned \
                         by os.listdir(indir)')
parser.add_argument('--start', dest='start',type=int,
                   help='First file to process.')
parser.add_argument('--stop', dest='stop',type=int,
                   help='Last file to process.')
args=parser.parse_args()

print("loading haplotype matrices")
haps,pos,labels,locs=get_ms_outs("sims/mutated/spatial/W16/",
                                 subset=args.subset,
                                 start=args.start,
                                 stop=args.stop)
pos=discretize_snp_positions(pos)

print("calculating summary statistics")
ss=getHaplotypeSumStats(haps,
                        pos,
                        labels,
                        locs,
                        args.outdir,
                        verbose=False)
