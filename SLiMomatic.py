import sys
import msprime as msp
import numpy as np
import subprocess
import pyslim
import keras
import os,io,shutil
import multiprocessing as mp
import subprocess as sp
import shlex
import scipy
import shutil
import random
from shutil import copyfile
import time
#import allel

class SLiMomatic(object):

    '''
    This class sets up the infrastructure to produce and analyze
    SLiM simulations in continuous space under random parameter
    values.

    See run_nwf_sim.py for an example use case.

    Initializing the class draws new labels from a uniform distribution
    bound by 'priorLow' and 'priorHigh' for each specified parameter
    (ie the variable name used for the parameter in your slim recipte).
    'simulateAndProduceTrees' runs SLiM simulations in parallel and
    recapitates the resulting tree sequences.

    You can then use 'simplifyTreeSequencesOnSubset' to
    draw a sample of individuals from the population and simplify the
    tree sequence, 'mutateTrees( )' to add mutations on to
    tree sequences, 'getHapsPosLabelsLocs( )' to output the
    haploytypes/positions/labels/coordinates to
    numpy arrays; and getSumStatsFromHaplotypes( ) to estimate a matrix
    of summary statistics.

    Everything here was adapted from the DataGenerator class in
    deeplearningtools/deepseq. Nothing has been tested in a rigorous way.
    Use at your own risk. God speed.

    '''

    def __init__(self,
        SLiMExecutableFilePath=None,
        SLiMRecipe=None,
        numReps=10,
        Ne = 1000, #this is the Ne used for recapitation
        priorLows=[0.5,0.5],
        priorHighs=[1.5,1.5],
        recombinationRate=1e-8,
        params_to_sample=["sigma_d","sigma_i"]):

        self.SLiMExecutableFilePath = SLiMExecutableFilePath
        self.SLiMRecipe = SLiMRecipe
        self.numReps = numReps
        self.Ne = Ne
        self.RecombinationRate = recombinationRate
        self.priorLows = priorLows
        self.priorHighs = priorHighs
        self.params_to_sample=params_to_sample
        self.postfix = ".trees"

        #get parameter values
        self.vals=np.empty((len(self.params_to_sample),numReps))
        for i in range(len(self.params_to_sample)):
            for j in range(numReps):
                val=np.random.uniform(self.priorLows[i],self.priorHighs[i])
                self.vals[i][j]=val

    def runOneSLiMSim(self,simNum,direc):
        '''
        run one SLiM simulation and put the corresponding treeSequence in treesOutputFilePath

        '''

        filename = str(simNum) + self.postfix
        filepath = os.path.join(direc,filename)

        #get strings for defining SLiM variables (only works for numeric variables)
        label_strs=[self.params_to_sample[i]+"="+str(self.vals[i][simNum])
                    for i in range(len(self.params_to_sample))]
        recomb_str="recomb="+str(self.RecombinationRate)
        output_str = "outpath='"+str(filepath)+"'"

        #format params for subprocess.check_output
        command=[self.SLiMExecutableFilePath,
                 "-d",recomb_str,
                 "-d",output_str]
        for i in range(len(self.params_to_sample)):
            command.append("-d")
            command.append(label_strs[i])
        command.append(self.SLiMRecipe)

        subprocess.check_output(command)

        #recapitate slim sims with msprime
        ts = pyslim.load(filepath)
        ts = ts.recapitate(recombination_rate=self.RecombinationRate, Ne=self.Ne)

        #write to file
        filename = str(simNum) + self.postfix
        filepath = os.path.join(direc,filename)
        ts.dump(filepath)

        return None

    def simulateAndProduceTrees(self,direc,nProc=1):
        '''
        run simulations in parallel
        '''

        #create the directory passed if it doesn't exits
        if not os.path.exists(direc):
            print("directory '",direc,"' does not exist, creating it")
            os.makedirs(direc)

        postfix = ".trees"
        labelsfilename = os.path.join(direc,"labels.txt")
        # if(os.path.isfile(labelsfilename)):
        #     print("Labels file already exists. Move or delete and try again.")
        #     return

        np.savetxt(labelsfilename,np.transpose(self.vals))

        # partition data for multiprocessing
        mpID = range(self.numReps)
        task_q = mp.JoinableQueue()
        result_q = mp.Queue()
        params=[postfix, direc]

        # do the work
        print("Simulating...")
        self.create_procs(nProc, task_q, result_q, params)
        self.assign_task(mpID, task_q, nProc)
        try:
            task_q.join()
        except KeyboardInterrupt:
            print("KeyboardInterrupt")
            sys.exit(0)

        return None


    def assign_task(self, mpID, task_q, nProcs):
        c,i,nth_job=0,0,1
        while (i+1)*nProcs <= len(mpID):
            i+=1
        nP1=nProcs-(len(mpID)%nProcs)
        for j in range(nP1):
            task_q.put((mpID[c:c+i], nth_job))
            nth_job += 1
            c=c+i
        for j in range(nProcs-nP1):
            task_q.put((mpID[c:c+i+1], nth_job))
            nth_job += 1
            c=c+i+1


    def create_procs(self, nProcs, task_q, result_q, params):
        for _ in range(nProcs):
            p = mp.Process(target=self.worker, args=(task_q, result_q, params))
            p.daemon = True
            p.start()


    def worker(self, task_q, result_q, params):
        while True:
            try:
                mpID, nth_job = task_q.get()
                #unpack parameters
                postfix, direc = params
                for i in mpID:
                    t0=time.time()
                    simNum = i
                    self.runOneSLiMSim(simNum,direc)
                    t=time.time()-t0
                    time_out=open("/projects/kernlab/cbattey2/spaceness/run_times.txt","a")
                    time_out.write(str(t)+" "+str(i)+"\n")
                    time_out.close()
                    ##result_q.put(outPut) #were to put output from the worker if needed
            finally:
                task_q.task_done()


#END OF CLASS
#---------------------------------------------------------------------------
#.TREE DIRECTORY HELPER FUNNCTIONS
def simplifyTreeSequenceOnSubSampleSet(ts,nSamples):
    '''
    This function should take in a tree sequence, generate
    a subset the size of numSamples, and return the tree sequence simplified on
    that subset of individuals
    '''
    sample_inds=np.unique([ts.node(j).individual for j in ts.samples()]) #final-generation individuals
    subsample=np.random.choice(sample_inds,nSamples,replace=False) #get nSamples random sample inds
    subsample_nodes=[ts.individual(x).nodes for x in subsample] #node numbers for sampled individuals
    subsample_nodes=[a for b in subsample_nodes for a in b] #flatten the list
    subsample_nodes=np.sort(np.array(subsample_nodes))
    ts=ts.simplify(subsample_nodes)

    return ts

def simplifyTreeSequenceDirectory(indir,outdir,nSamples):
    ntrees=len([f for f in os.listdir(indir) if not f.startswith(".")])-1
    trees=[indir+str(i)+".trees" for i in range(ntrees)]
    for i in range(len(trees)):
        t=pyslim.load(trees[i])
        o=simplifyTreeSequenceOnSubSampleSet(ts=t,nSamples=nSamples)
        o.dump(outdir+str(i)+".trees")
    shutil.copyfile(os.path.join(indir,"labels.txt"),os.path.join(outdir,"labels.txt"))

    return None

def mutateTrees(treesDirec,outputDirec,numMuts,muLow,muHigh):
    '''
    read in .trees files from treesDirec, mutate that tree numMuts seperate times
    using a mutation rate pulled from a uniform dirstribution between muLow and muHigh

    also, re-write the labels file to reflect.
    '''
    if not os.path.exists(outputDirec):
        print("directory '",outputDirec,"' does not exist, creating it")
        os.makedirs(outputDirec)

    # labelsFilename = os.path.join(treesDirec,"labels.txt")
    # lablesFile = open(labelsFilename,"r")
    # labels = np.array(lablesFile.readline().split(),dtype='float32')
    #
    # newLabelsFilename = os.path.join(outputDirec,"labels.txt")
    # newLabelsFile = open(newLabelsFilename,"w")
    # newLabels = []

    #how many trees files are in this directory.
    numReps=len([f for f in os.listdir(treesDirec) if not f.startswith(".")])-1

    for i in range(numReps):
        filename = str(i) + ".trees"
        filepath = os.path.join(treesDirec,filename)
        treeSequence = pyslim.load(filepath)
        blankTreeSequence = msp.mutate(treeSequence,0)
        #rho = labels[i]
        for mut in range(numMuts):
            simNum = (i*numMuts) + mut
            simFileName = os.path.join(outputDirec,str(simNum)+".trees")
            mutationRate = np.random.uniform(muLow,muHigh)
            mutatedTreeSequence = msp.mutate(blankTreeSequence,mutationRate)
            mutatedTreeSequence.dump(simFileName)
            #newLabels.append(rho)


    shutil.copyfile(os.path.join(treesDirec,"labels.txt"),os.path.join(outputDirec,"labels.txt"))

    return None

def getBranchLengthSumStats(direc,n): #don't use this yet
    #read through a directory of SLiM tree sequences and return a matrix
    #of summary stats calculated from branch lengths
    out=[]
    for i in range(n):
        filename = str(i) + ".trees"
        filepath = os.path.join(direc,filename)
        ts = pyslim.load(filepath)
        samples=[x for x in ts.samples()]
        locs=[[ts.individual(ts.node(s).individual).location[0],
               ts.individual(ts.node(s).individual).location[1]] for s in samples]
        bs=msp.BranchLengthStatCalculator(ts)
        gd=np.array(bs.divergence([[x] for x in ts.samples()],
                                 windows=[0.0,ts.sequence_length]))
        gd = gd[np.logical_not(np.isnan(gd))]
        sd = np.array(scipy.spatial.distance.pdist(locs))
        simSFS = np.array(bs.site_frequency_spectrum(samples)[0])
        out.append(np.concatenate((sd,gd,simSFS)))
    return out

def getHapsPosLabelsLocs(direc):
    '''
    loops through a trees directory created by the data generator class
    and returns the repsective genotype matrices, positions, and labels
    '''
    haps = []
    positions = []
    labels=np.loadtxt(os.path.join(direc,"labels.txt"))
    locs = []

    ntrees=np.shape(labels)[0]

    for i in range(ntrees):
        filename = str(i) + ".trees"
        filepath = os.path.join(direc,filename)
        ts = pyslim.load(filepath)
        haps.append(ts.genotype_matrix())
        positions.append(np.array([s.position for s in ts.sites()]))
        sample_inds=np.unique([ts.node(j).individual for j in ts.samples()])
        locs.append([[ts.individual(x).location[0],
               ts.individual(x).location[1]] for x in sample_inds])

    haps = np.array(haps)
    positions = np.array(positions)
    locs=np.array(locs)

    return haps,positions,labels,locs

#discretize positions (ie if SNPs are in the same integer position move one to the right)
def discretize_snp_positions(positions):
    '''
    Takes an array of SNP positions as floats (ie from msprime) and returns
    integer positions. If two SNPs fall in the same integer position, one is
    shifted to the right by one base pair.
    '''
    count=0
    for i in range(len(positions)):
        dpos=[int(x) for x in positions[i]]
        for j in range(len(dpos)):
            if(dpos[j]==dpos[j+1]):
                dpos[j+1]=dpos[j+1]+1
                count=count+1
        positions[i]=dpos
    print(str(count)+" SNPs were shifted one base pair")
    return(positions)


def getPairwiseIbsTractLengths(x,y,positions,maxlen):
    '''
    input:
    x: haplotype as 1D array
    y: haplotype as 1D array
    positions: a 1D array of SNP positions (integer, 1-based)

    Returns:
    1d array listing distances between adjacent SNPs in a pair of sequences.
    Note IBS blocks returned are *not* in order(!)

    Example:
    haps,positions,labels,locs=getHapsPosLabelsLocs(direc="mutated_treeSeq_directory")
    ibs_lengths_1_2=getPairwiseIbsTractLengths(haps[0][:,0],haps[0][:,1],positions[0])
    '''
    snps=~np.equal(x,y)
    snp_positions=positions[snps]
    l=len(snp_positions)
    if(l==0):
        ibs_tracts=maxlen
    else:
        if(l>1):
            ibs_tracts=snp_positions[np.arange(1,l,2)]-snp_positions[np.arange(0,l-1,2)] #middle blocks
        np.append(ibs_tracts,snp_positions[0])          #first block
        np.append(ibs_tracts,maxlen-snp_positions[l-1]) #last block
    return ibs_tracts
