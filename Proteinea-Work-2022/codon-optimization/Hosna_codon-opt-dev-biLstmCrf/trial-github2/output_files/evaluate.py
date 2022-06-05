# -*- coding: utf-8 -*-
"""
Created on Tue May 10 14:27:03 2022

@author: hosna
"""

# Calculate mismatch on character-level not codon level which already denotes error in codon level (if 1 happened then skip next 2 and add 1)
#
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import seaborn as sns
import string
import itertools

def len_checker(file_true, file_pred):
    """Check the length of each sequence from both predicted and true files (sequence files) are same or not

    Args:
        file_true (_io.TextIoWrapper): file of the true primary sequences
        file_pred (_io.TextIoWrapper): file of the remapped predicted primary sequences
    """
    for seq_true, seq_pred in zip(file_true, file_pred):
        # Same length checker
        if len(seq_true) == len(seq_pred): # guaranteed after adding gaps for invalid
            print("Both sequences are of same length")
        else:
            print(f"length of the true test line {len(seq_true)}, while length of the predicted test line is {len(seq_pred)}.") 

def evaluate(file_true, file_pred, fixed_file_pred, file_pred_codon, file_true_codon):
    """Evaluate the predictions of the model by: 
       1. Calculate number of mismatched nucleotides per sequence
       2. Calculate number of invalid mappings (between codon boxes and AAs)
       3. Fix the gapped nucleotides of each sequence (by returning the old codons)
       4. Calculate frequency of codon boxes per sequence (and overall)
       5. 

    Args:
        file_true (_io.TextIoWrapper): file of the true primary sequences
        file_pred (_io.TextIoWrapper): file of the remapped predicted primary sequences
        fixed_file_pred (_io.TextIoWrapper): file to save the fixed predicted primary sequences (no gaps)
        file_pred_codon (pandas.DataFrame): dataframe for the predicted file with AAs and predicted codon boxes
    """
    mismatch = 0
    file_true_list = []
    file_pred_list = []
    # for i, j in enumerate(file_true):
    #     pass
    # print(f'Number of True sequences in sequence file is {i+1}')
    # print('end')
    # for i, j in enumerate(file_pred):
    #     pass
    # print(f'Number of Pred sequences in sequence file is {i+1}')
    x = file_true_codon.groupby(['amino_acid'])['codon_box'].agg(pd.Series.mode)
    for i, (seq_true, seq_pred) in enumerate(zip(file_true, file_pred)):
        print(i)
        seq_pred_copy = (seq_pred + '.')[:-1] # to evade a copy by reference, make new string and slice
        file_true_list.append((seq_true))
        file_pred_list.append((seq_pred))
        print(f"for sequence {i}")
        mismatch = 0

        # Invalid nucleotides
        invalid = (seq_pred.count('-'))/3
        print(f"Number of invalid mappings are: {invalid}")

        # Mismatch and fix
        # mismatch[len(seq_pred)] = sum(1 for x, y in zip(seq_true, seq_pred) if x != y) # (true length of sequence : number of mismatch)
        for i, (nucl_true, nucl_pred) in enumerate(zip(seq_true, seq_pred)):
            # per-nucleotide Mismatch
            if nucl_pred != nucl_true:
                mismatch +=1 
            # Q: if it's a gap, we return the same old codon?
            
            # Invalid fallback => try on BSF
            # DataFrame for true mapping => AA, CodonBox, frequency
            # groupby AA: get top frequent codonboxes for each AA 
            # save in dict or any [AA:codonbox]
            if nucl_pred == '-':
                seq_pred_copy = seq_pred_copy[:i] + nucl_true + seq_pred_copy[i+1:]

        print(f"Number of mismatched nucleotides is {mismatch} out of {len(seq_pred)} total nucleotides")
        fixed_file_pred.writelines(seq_pred_copy + '\n')

    my_list = [len(x)-len(y) for i, (x, y) in enumerate(zip(file_true_list, file_pred_list))]    
    # uni-Frequency of predicted codon-box
    file_pred_codon = file_pred_codon[file_pred_codon['amino_acid'] != '#']
    print('Number of True seqs from AA "X" id {0}'.format(file_true_codon['amino_acid'].value_counts()[-1]))
    print('Number of Pred seqs from AA "X" id {0}'.format(file_pred_codon['amino_acid'].value_counts()[-1]))
    n_pred = file_pred_codon['amino_acid'].value_counts()
    # n_true = 
    # To check the frequency in numbers
    # pred_uni_freq = file_pred_codon['codon_box'].value_counts()
    # true_uni_freq = file_true_codon['codon_box'].value_counts()
    plt.figure(0)
    file_true_codon['codon_box'].hist(alpha = 0.5, grid=False)
    file_pred_codon['codon_box'].hist(alpha=0.5, grid=False)
    plt.legend(["True", "Predicted"])
    plt.title("Unigram Distribution")
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    
    fig = plt.figure(1)
     
    """
    Bigram start
    """
    # initializing string
    str_pred_codon = ''.join(file_pred_codon['codon_box'])
    str_true_codon = ''.join(file_true_codon['codon_box'])
    
    # get all possible combinations between true and predicted
    freq_pair_true = {}
    for idx in range(len(str_true_codon)): # may fail if != lengths; check drosophila
        freq_pair_true[str_pred_codon[idx : idx+2]] = 0
        freq_pair_true[str_true_codon[idx : idx+2]] = 0  
        
    freq_pair_pred = freq_pair_true.copy()
    
    # get values/frequency
    for idx in range(len(str_true_codon)):
        freq_pair_pred[str_pred_codon[idx : idx+2]] += 1
        freq_pair_true[str_true_codon[idx : idx+2]] += 1
    
    # unzip
    k_p, v_p = zip(*freq_pair_pred.items())
    k_t, v_t = zip(*freq_pair_true.items())
    # normalize
    v_p_norm = tuple(((v_n-min(v_p))/(max(v_p)-min(v_p))) for v_n in v_p)
    v_t_norm = tuple(((v_n-min(v_t))/(max(v_t)-min(v_t))) for v_n in v_t)
    
    # what we see here are only all the codon-boxes combination between true and predicted not all the combinations ever
    # all combinations ever would be about 400 (codon box: 20, 64 codons). However, here we have 206
    # plot in dataframe
    bi_p = pd.DataFrame({'1st_codonbox': [x[0] for x in k_p],
                          '2nd_codonbox': [x[-1] for x in k_p],
                          'Frequency': v_p_norm
                     })
    bi_t = pd.DataFrame({'1st_codonbox': [x[0] for x in k_t],
                          '2nd_codonbox': [x[-1] for x in k_t],
                          'Frequency': v_t_norm
                     })

    heatmap_p = bi_p.pivot_table(index='1st_codonbox', columns='2nd_codonbox', values='Frequency')
    sns.heatmap(heatmap_p).set_title('Predicted Bigram')
    plt.show()
    
    heatmap_t = bi_t.pivot_table(index='1st_codonbox', columns='2nd_codonbox', values='Frequency')
    sns.heatmap(heatmap_t).set_title('True Bigram')
    plt.show()
    
    """
    Bigram End
    """

    # plt.fill_between(k_t, v_t_norm, step='post')
    # plt.fill_between(k_p, v_p_norm)
    # plt.plot(k_t, v_t_norm)
    # plt.plot(k_p, v_p_norm)
    # plt.legend(["True", "Predicted"])
    # plt.xticks(color='w')
    # plt.title('Bigram Distribution')
    # plt.xlabel("Codon Pair")
    # plt.ylabel("Frequency")

# file_true = open("dros_all/dros_true_dev_seq.txt") # H: original
# file_pred = open("dros_all/dros_pred_dev_seq.txt") # H: predicted
# fixed_file_pred = open("dros_all/pred_seq_dev_fixed.txt", 'w+') # H: file to save fixed predictions at (no gaps)
# file_pred_codon = pd.read_csv("dros_all/dros_pred_dev_codon.txt", sep = ' ', 
#                                 names = ['amino_acid', 'codon_box']) # H: predicted codons
# file_true_codon = pd.read_csv("dros_all/dros_true_dev_codon.txt", sep = ' ', 
#                                 names = ['amino_acid', 'codon_box']) # H: true codons

file_true = open("bsf_all/bsf_true_seq.txt") # H: original
file_pred = open("bsf_all/bsf_pred_seq.txt") # H: predicted
fixed_file_pred = open("bsf_all/pred_seq_fixed.txt", 'w+') # H: file to save fixed predictions at (no gaps)
file_pred_codon = pd.read_csv("bsf_all/bsf_pred_codon.txt", sep = ' ', 
                                names = ['amino_acid', 'codon_box']) # H: predicted codons
file_true_codon = pd.read_csv("bsf_all/bsf_true_codon.txt", sep = ' ', 
                                names = ['amino_acid', 'codon_box']) # H: true codons

evaluate(file_true, file_pred, fixed_file_pred, file_pred_codon, file_true_codon)


# Calculate the frequency of each character of codon boxes (blank/3: invalid mapping for AA (codon, 3 nucleotide))
# extract the 3rd columns values and save them to list then plot

# todo: change open to with open
# todo: make as argparse to input the path only
# todo: make if __main__

