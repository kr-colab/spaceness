import argparse, os, subprocess, numpy

parser = argparse.ArgumentParser(description='Run a SLiM simulation while sampling a \
                                              parameter from a uniform distribution.')
parser.add_argument('--sampled_param', dest='sampled_param',
                    help='Name of parameter to sample \
                    (should match the variable name in the slim recipe).')
parser.add_argument('--min', dest='min',type=float,
                   help='Minimum of sampled range.')
parser.add_argument('--max', dest='max',type=float,
                   help='Maximum of sampled range.')
parser.add_argument('--slim_path', dest='slim_path',
                   help='call for slim.')
parser.add_argument('--slim_recipe', dest='slim_recipe',
                   help='full path to slim recipe to run.')
parser.add_argument('--outdir', dest='outdir',
                   help='full path to output file directory.')
args=parser.parse_args()

def run_one_slim_sim(sampled_param,
                     min,
                     max,
                     slim_path,
                     slim_recipe,
                     outdir):
    '''
    Run one SLiM simulation pulling parameter values from a uniform
    distribution. Output will be named with the name and value of
    the sampled parameter.
    '''

    #get sampled param names and values
    val=numpy.random.uniform(min,max)

    #get output file path
    filename = sampled_param+"_"+str(val)+"_.trees"
    filepath = os.path.join(outdir,filename)

    #set strings for defining SLiM variables
    label_str=sampled_param+"="+str(val)
    output_str = "outpath='"+str(filepath)+"'"

    #format params for subprocess.check_output
    command=[slim_path,
             "-d",output_str,
             "-d",label_str,
             slim_recipe]

    #run it
    print("starting slim simulation with "+args.sampled_param+"="+str(val))
    subprocess.check_output(command)

    return None

run_one_slim_sim(args.sampled_param,
                 args.min,
                 args.max,
                 args.slim_path,
                 args.slim_recipe,
                 args.outdir)
