# -*- coding: utf-8 -*-
"""
Created on Tue May 10 14:27:03 2022

@author: hosna
"""
# For drosophila
# python3 count_mismatch.py -ts dros_true_seq.txt -ps dros_pred_seq.txt -o pred_seq_fix.txt 
# -tc dros_true_codon.txt -pc dros_pred_codon.txt

# Calculate mismatch on character-level not codon level which already denotes error in codon level (if 1 happened then skip next 2 and add 1)
#
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

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

def mismatch_fix(file_true, file_pred, fixed_file_pred):
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
    
    for i, (seq_true, seq_pred) in enumerate(zip(file_true, file_pred)):
        seq_pred_copy = (seq_pred + '.')[:-1] # to evade a copy by reference, make new string and slice
        print(f"for sequence {i}")

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
            if nucl_pred == '-':
                seq_pred_copy = seq_pred_copy[:i] + nucl_true + seq_pred_copy[i+1:]

        print(f"Number of mismatched nucleotides is {mismatch} out of {len(seq_pred)} total nucleotides")
        fixed_file_pred.writelines(seq_pred_copy + '\n')

def codon_distribution(file_pred_codon, file_true_codon):

    # uni-Frequency of predicted codon-box
    file_pred_codon = file_pred_codon[file_pred_codon['amino_acid'] != '#']
    # To check the frequency in numbers
    # pred_uni_freq = file_pred_codon['codon_box'].value_counts()
    # true_uni_freq = file_true_codon['codon_box'].value_counts()
    plt.figure(0)
    file_true_codon['codon_box'].hist(alpha = 0.5, grid=False)
    file_pred_codon['codon_box'].hist(alpha=0.5, grid=False)
    plt.legend(["True", "Predicted"])
    plt.title("Codon's Distribution")
    
    # Bigram Pred
    plt.figure(1) 
    str_pred_codon = ''.join(file_pred_codon['codon_box'])
    pred_bigram = Counter(str_pred_codon[idx : idx + 2] for idx in range(len(str_pred_codon) - 1))
    pred_bigram = dict(pred_bigram)
    k, v = zip(*pred_bigram.items())
    plt.plot(k, v)
    plt.xticks(color='w')
    plt.title('Bigram Distribution')
    plt.xlabel("Codon Pair")
    plt.ylabel("Frequency")
    
    # Bigram True
    plt.figure(2)
    str_true_codon = ''.join(file_true_codon['codon_box'])
    true_bigram = Counter(str_true_codon[idx : idx + 2] for idx in range(len(str_true_codon) - 1))
    true_bigram = dict(true_bigram)
    k, v = zip(*true_bigram.items())
    plt.plot(k, v)
    plt.xticks(color='w')
    plt.title('Bigram Distribution')
    plt.xlabel("Codon Pair")
    plt.ylabel("Frequency")
    
    
    plt.show()

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ts", "--true_seq", default='true_txt.txt',type=str, help="path to true primary sequence")
    parser.add_argument("-ps", "--pred_seq", type=str, help="path to predicted primary sequence")
    parser.add_argument("-o", "--pred_seq_fix", type=str, help="path to output the fixed predicted primary sequence")
    parser.add_argument("-tc", "--true_codon", type=str, help="path to true codon box")
    parser.add_argument("-pc", "--pred_codon", type=str, help="path to predicted codon box")
    args = parser.parse_args()

    return args

def main():

    args = arguments()
    file_true = open(args.true_seq) # H: original
    file_pred = open(args.pred_seq) # H: predicted
    fixed_file_pred = open(args.pred_seq_fix, 'w+') # H: file to save fixed predictions at (no gaps)

    file_true_codon = pd.read_csv(args.true_codon, sep = ' ', names = ['amino_acid', 'codon_box']) # H: true codons
    file_pred_codon = pd.read_csv(args.pred_codon, sep = ' ', names = ['amino_acid', 'codon_box']) # H: predicted codons

    len_checker(file_true, file_pred)
    print('done 1')
    mismatch_fix(file_true, file_pred, fixed_file_pred)
    print('done 2')
    codon_distribution(file_pred_codon, file_true_codon)
    print('done 3')

if __name__ == "__main__":
    main()

# todo: change open to with open


