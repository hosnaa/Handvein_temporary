# 
# Write a file with the amino acid and box symbol for the input fasta file of sequences

# To run this file
## python codonToBox_v5_argparse.py -i cbox.fasta

from Bio import SeqIO
from sklearn.model_selection import train_test_split
from typing import TypeVar
import argparse

seqRecord = TypeVar('seqRecord')  # Custom type hinting, it's a seqRecord type from BioIO 
# key: codon, value: [AA, symbol_for_box]
CODON_BOX_LOOKUP = {'TTT': ['F', 'a'],
                    'TTC': ['F', 'b'], 
                    'TTA': ['L', 'c'],
                    'TTG': ['L', 'd'],
                    'TCT': ['S', 'b'],
                    'TCC': ['S', 'f'],
                    'TCA': ['S', 'h'],
                    'TCG': ['S', 'g'],
                    'TAT': ['Y', 'c'], 
                    'TAC': ['Y', 'h'],
                    'TAA': ['X', 'w'],
                    'TAG': ['X', 'w'],
                    'TGT': ['C', 'd'],
                    'TGC': ['C', 'g'],
                    'TGA': ['X', 'w'], 
                    'TGG': ['W', 'k'],
                    'CTT': ['L', 'b'],
                    'CTC': ['L', 'f'],
                    'CTA': ['L', 'h'],
                    'CTG': ['L', 'g'],
                    'CCT': ['P', 'f'],
                    'CCC': ['P', 'l'],
                    'CCA': ['P', 'm'],
                    'CCG': ['P', 'n'],
                    'CAT': ['H', 'h'],
                    'CAC': ['H', 'm'],
                    'CAA': ['Q', 'o'],
                    'CAG': ['Q', 'r'],
                    'CGT': ['R', 'g'],
                    'CGC': ['R', 'n'],
                    'CGA': ['R', 'r'],
                    'CGG': ['R', 's'],
                    'ATT': ['I', 'c'],
                    'ATC': ['I', 'h'],
                    'ATA': ['I', 'i'],
                    'ATG': ['M', 'j'],
                    'ACT': ['T', 'h'],
                    'ACC': ['T', 'm'],
                    'ACA': ['T', 'o'],
                    'ACG': ['T', 'r'],
                    'AAT': ['N', 'i'],
                    'AAC': ['N', 'o'],
                    'AAG': ['K', 'u'],
                    'AAA': ['K', 't'],
                    'AGT': ['S', 'j'],
                    'AGC': ['S', 'r'],
                    'AGA': ['R', 'u'],
                    'AGG': ['R', 'q'],
                    'GTT': ['V', 'd'],
                    'GTC': ['V', 'g'],
                    'GTA': ['V', 'j'],
                    'GTG': ['V', 'k'],
                    'GCT': ['A', 'g'],
                    'GCC': ['A', 'n'],
                    'GCA': ['A', 'r'],
                    'GCG': ['A', 's'],
                    'GAT': ['D', 'j'],
                    'GAC': ['D', 'r'],
                    'GAA': ['E', 'u'],
                    'GAG': ['E', 'q'],
                    'GGT': ['G', 'k'],
                    'GGC': ['G', 's'],
                    'GGA': ['G', 'q'],
                    'GGG': ['G', 'p']                      
}

def read_split_fasta(input_fasta_path: str, test_size=0.1, val_size=0.1) -> seqRecord:
    """Read the sequences from the fasta file then split them into train, val and test sets.

    Args:
        input_fasta_path (str): path of the fasta file to train on
        test_size (float, optional): size of test set used in splitting between train and test. Defaults to 0.1.
        val_size (float, optional): size of val set used in splitting between the resultant train and val. Defaults to 0.1.

    Returns:
        seqRecord: 3 seqRecords, one for each set where each contain the sequence (.seq) and the ID (.id) of sequence
    """
    seq_tot = [x for x in SeqIO.parse(input_fasta_path,"fasta") if len(x) <= 500] 
    _, _, seq_train, seq_test = train_test_split(seq_tot, seq_tot, test_size = test_size, random_state = 1)
    _, _, seq_train, seq_val = train_test_split(seq_train, seq_train, test_size = val_size, random_state = 1) 

    # (Optional) save test sequences as is for checking afterwards
    orig_test_file = open('sample_data/predictions/orig_test_file.txt', 'w+')
    for i in seq_test:
        orig_test_file.writelines(i + "\n")

    return seq_train, seq_val, seq_test

def codon_to_box(all_sequences: seqRecord, output_path: str):
    """Take the seqRecords (sequence, id) of the whole set and convert the codons of sequence into their 
       corresponding amino acids and a corresponding symbol denoting the box

    Args:
        all_sequences (seqRecord): all the sequences in the set (train/val/test), has .seq and .id
        output_path (str): path of the file to save in it the mapping between amino acid and box symbol
    """
    output_file = open(output_path, 'w+')
    for one_seq in all_sequences:

        i = 0  # Q: Is this the same as initializing i=0 before for loop and then reinitialize it after while loop ends

        # loop over sequences. for each sequence: Extract seq only from each, make mapping for each codon
        # save this mapping in a file 
        while i < len(one_seq.seq): 
            codon = one_seq.seq[i] + one_seq.seq[i+1] + one_seq.seq[i+2] # slice to take consecutive codons (3 nucleotides) from the sequence
            
            try:
                amino_acid, box_symbol = CODON_BOX_LOOKUP[codon]
                """
                Dataframe => AA: [codon, frequency] 
                          => groupby AA, get max(freq) and its corresponding codon
                          => save the AA and max(freq)

                write file: dict{AA, most frequent codon}
                """
                amino_box = amino_acid + " " + box_symbol
                output_file.writelines(amino_box + "\n")
                i = i + 3

                if amino_acid == 'X' :  # X -> Stop codon: denoting end of sequence and we're about to being a new sequence
                    output_file.write('\n')

            except ValueError: # todo: check how to modify these exceptions
                print('Can not find this codon in the LookUp dictionary, please check your sequence')

    output_file.close()

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fasta_path", type=str, default = "dros.fasta",
                        help="Path of the input fasta file")
    args = parser.parse_args()

    return args

def main():
    args = arguments()
    
    # path: "D:/Proteinea-Work-2022/codon-optimization/Hosna_codon-opt-dev/bsf_old_data_cds_filtered_05.fasta"
    input_fasta_path = args.input_fasta_path

    output_train_path = input_fasta_path[ :input_fasta_path.rfind('.')] + "_train.txt"
    output_test_path = input_fasta_path[ :input_fasta_path.rfind('.')] + "_test.txt"
    output_val_path = input_fasta_path[ :input_fasta_path.rfind('.')] + "_dev.txt"

    # Read and Split
    seq_train, seq_val, seq_test = read_split_fasta(input_fasta_path)

    # Mapping 
    codon_to_box(seq_train, output_train_path)
    codon_to_box(seq_val, output_val_path)
    codon_to_box(seq_test, output_test_path)

if __name__ == "__main__":
    main()







