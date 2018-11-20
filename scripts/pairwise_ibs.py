import numpy as np
def getPairwiseIbsTractLengths(x,y,positions,maxlen):
    '''
    input:
    x: haplotype as 1D array
    y: haplotype as 1D array
    positions: a 1D array of SNP positions (integer, 1-based)

    Returns:
    1d array listing distances between adjacent SNPs in a pair of sequences.

    Example:
    haps,positions,labels,locs=getHapsPosLabelsLocs(direc="mutated_treeSeq_directory")
    ibs_lengths_1_2=getPairwiseIbsTractLengths(haps[0][:,0],haps[0][:,1],positions[0])
    '''
    snps=np.equal(x,y)
    snp_positions=positions[snps]
    l=len(snp_positions)
    ibs_tracts=[]
    if(l==0):
        ibs_tracts=[maxlen]
    else:
        if(l>1):
            ibs_tracts=snp_positions[np.arange(1,l,2)]-snp_positions[np.arange(0,l-1,2)] #middle blocks
        np.append(ibs_tracts,snp_positions[0])          #first block
        np.append(ibs_tracts,maxlen-snp_positions[l-1]) #last block
    return ibs_tracts
