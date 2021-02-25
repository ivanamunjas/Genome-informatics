# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:34:51 2020

@author: Ivana Pesic
"""

#import time, psutil, os, sys
import sys
import numpy as np
import global_alignment
import fm_index
import pickle
import os

path = os.path.dirname(os.path.abspath(__file__))

# read input parameters from command line
def parameters():
    
    args = sys.argv
    reference_name = args[1]
    reads_name = args[2]
    seed = int(args[3])
    match = int(args[4])
    mismatch = int(args[5])
    gap = int(args[6])
    margin = int(args[7])
    
    if margin > 3 or margin < 0:
        print('margin has to be 0-3')
        exit()

    return reference_name, reads_name, seed, match, mismatch, gap, margin

# read genom
def read_reference(reference):

    file = open(reference, 'r')
    file.readline()
    reference = file.read().replace('\n', '')
    reference = reference + '#'
    file.close()
    return reference


def process_read(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)
    
def reverse_complement(s):
    return complement(s[::-1])
        
# extract reads from reads.fastq and make reverse complements
def return_reads_complements(reads_name):
        
    n = 4
    
    with open(reads_name, 'r') as reads_file:
        lines = []
        reads = []
        complements = []
        for line in reads_file:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process_read(lines)
                sys.stderr.write("Record: %s\n" % (str(record)))
                reads.append(record["sequence"])
                complements.append(reverse_complement(record["sequence"]))
                lines = []
                
    return reads, complements            

# calculate FM index for all reads and reverse complements
def bw_fm(reference, reads, complements, seed):

    idx = fm_index.FMCheckpointing(reference)
    
    indices = []
    for i in range(len(reads)):
        indices.append(idx.search(reads[i][0:seed]))
        indices.append(idx.search(complements[i][0:seed]))
    
    return indices


def main():
    
    # ask for input parameters
    reference_name, reads_name, seed, match, mismatch, gap, margin = parameters()
    
    # process reference file 
    reference = read_reference(path + '\\' + reference_name)
    
    # extract reads and get their complements
    reads, complements = return_reads_complements(path + '\\' + reads_name)
    
    # get BW transformations and FM indices
    indices = bw_fm(reference, reads, complements, seed)
    
    # list to store fm indices, alignment scores and edit transcripts
    results = []
    
    # for each read and its complement get results
    for i in range(2*len(reads)):
        inds = []
        alignment_scores = []
        transcripts = []
        for ind in indices[i]:
            if np.mod(i,2) == 0:
                seq = reads[i//2]
            else:
                seq = complements[i//2]
            align = global_alignment.GlobalAlignment(match, mismatch, gap)
            alignment_score, edit_transcript = \
            align.return_parameters(reference[ind:ind+len(seq)+margin], seq)
            inds.append(ind)
            alignment_scores.append(alignment_score)
            transcripts.append(edit_transcript)
        # sort and save results
        alignment_scores = np.array(alignment_scores)
        inds = np.array(inds)
        transcripts = np.array(transcripts)
        ind_sorted = alignment_scores.argsort()
        ind_sorted = np.flip(ind_sorted)
        results.append([inds[ind_sorted], alignment_scores[ind_sorted], transcripts[ind_sorted]])

    return reference, reads, results, indices, match, mismatch, gap


reference, reads, results, indices, match, mismatch, gap = main()

# save results variable
with open(path+ '\\results' + str(match) + str(mismatch) + str(gap) + '.pickle', 'wb') as f:
    pickle.dump(results, f)

