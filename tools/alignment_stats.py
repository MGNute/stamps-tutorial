import os, sys
import numpy as np

def read_from_fasta(file_path):
    """
    Reads from a fasta file and returns a dictionary where keys are taxon names and values are strings
    representing the sequence.
    :param file_path (string): The full system-readable path to the fasta file
    :return: fasta (dict)
    """
    output={}
    fasta=open(file_path,'r')
    first=True
    seq=''
    for l in fasta:
        if l[0]=='>':
            if first!=True:
                output[name]=seq
            else:
                first=False
            name=l[1:].strip()
            seq=''
        else:
            seq=seq + l.strip()
    output[name]=seq
    fasta.close()
    return output

def fasta_dict_to_nparray(fasta, taxnames=None, order='C'):
    '''
    Converts a dictionary representing an *aligned* fasta file (keys = seq names, values = seqs),
    into a numpy array of type np.uint8. Returns (array, taxa_names) where taxa_names is a list
    of names, in order, represented by the rows of the array. Should match the order of fasta.keys,
    but just in case...
    :param fasta: 
    :param taxnames: 
    :param order: 
    :return: 
    '''
    ntax = len(fasta.keys())
    ncols = max(map(len,fasta.values()))
    nparr = np.zeros((ntax,ncols),dtype=np.uint8,order=order)
    if taxnames is None:
        taxnames=[]
        # for i in range(ntax):
        i=0
        for k in fasta.keys():
            taxnames.append(k)
            seq = fasta[k]
            ls = len(seq)
            nparr[i,:ls]=np.frombuffer(bytes(seq,'utf8'),dtype=np.uint8)
            if ls < ncols:
                nparr[i,ls:]=45
            i+=1
    else:
        for i in range(len(taxnames)):
            seq = fasta[taxnames[i]]
            ls = len(seq)
            nparr[i,:ls]=np.frombuffer(bytes(seq,'utf8'),dtype=np.uint8)
            if ls < ncols:
                nparr[i,ls:]=45
    if order=='F':
        nparrF=np.array(nparr, copy=True, order='F')
        return taxnames, nparrF
    else:
        return taxnames, nparr

def get_avg_pdistance_of_nparray(fnp, getmax=False, weighted=False):
    '''
    Computes average pairwise p-distance from an NP-array of uint8 representing an MSA.
    :param fnp:
    :param getmax: if True, return the max p-distance isntead of the average (default: False)
    :param weighted: if True, return the weighted p-distance instead of the raw average
    :return: (pd, equal_site_ct, common_site_ct, pair_ct)
    '''

    # f = read_from_fasta(fasta_path)
    # taxn, fnp = fasta_dict_to_nparray(f)

    maxpd = 0.
    run_tot = 0.
    run_ct = 0
    comm_sum=0
    same_sum=0
    for i in range(fnp.shape[0]):
        for j in range(i):
            comm = np.sum((fnp[i,:]!=45) & (fnp[j,:]!=45),dtype=np.float64)
            if comm ==0:
                continue
            same = np.sum((fnp[i,:]!=45) & (fnp[j,:]!=45) & (fnp[i,:]==fnp[j,:]),dtype=np.float64)
            run_tot+= 1.- same / comm
            run_ct += 1
            same_sum += same
            comm_sum += comm
            if (1. - same / comm > maxpd):
                maxpd = 1. - same / comm
            # print('i: %s\tj: %s\tcomm: %s\tsame: %s' % (i,j,comm,same))

    # print ('Avg P-Distance: %s' % pd)
    if getmax==False:
        if weighted:
            pd = 1.0 - same_sum / comm_sum
            return pd, int(same_sum), int(comm_sum), run_ct
        else:
            pd = run_tot / float(run_ct)
            return pd, int(same_sum), int(comm_sum), run_ct
    else:
        return maxpd

if __name__ == "__main__":
    fa_pa = sys.argv[1]
    fa = read_from_fasta(fa_pa)
    tn, fnp = fasta_dict_to_nparray(fa)
    seq_len = fnp.shape[1]
    pd, match_site_pairs, aligned_site_pairs, num_seq_pairs = get_avg_pdistance_of_nparray(fnp)
    print("file: %s" % fa_pa)
    print("  Alignment Length:        %d" % seq_len)
    print('  Avg %% similarity:        %0.4f' % (1.0 - pd))
    print("  # Matched Site-pairs:    %d" % match_site_pairs)
    print("  # Aligned Site-pairs:    %d" % aligned_site_pairs)
