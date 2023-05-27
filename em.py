# -*- coding: utf-8 -*-

'''
goal is to predict whether a sequence of DNA will be bound by a specific transcription factor in a given condition
'''

import numpy as np
import random
import argparse

# Create the argument parser
parser = argparse.ArgumentParser(description='File paths')

# Add the file path arguments
parser.add_argument('bound_file', type=str, help='Path to the bound file')
parser.add_argument('test_file', type=str, help='Path to the test file')

# Parse the command-line arguments
args = parser.parse_args()


#function to read in fasta files 
def read_fasta_file(filename):
    with open(filename, 'r') as f:
        read_id = None
        read_seq = ""
        for line in f:
            if line.startswith(">"):
                if read_id is not None:
                    yield (read_id, read_seq)
                read_id = line.strip()[1:]
                read_seq = ""
            else:
                read_seq += line.strip()
        if read_id is not None:
            yield (read_id, read_seq)

#Computes the pseudocount profile matrix of a given set of motifs.
def compute_pseudocount_profile_matrix(motifs):
    motif_length = len(motifs[0])  # Length of each motif
    num_motifs = len(motifs)  # Number of motifs
    motif_scale = 1 / (num_motifs + 4)  # Frequency scale factor, accounting for pseudocounts

    seq_to_numeric = 'ACGT0123'  # Mapping of nucleotides to numeric values
    numeric_dict = {seq_to_numeric[i]: int(seq_to_numeric[i + 4]) for i in range(4)}

    # Initialize the pseudocount profile matrix with scaled pseudocounts
    pseudocount_matrix = [[motif_scale for _ in range(motif_length)] for __ in range(4)]

    # Count the occurrences of each nucleotide at each position in the motifs and apply pseudocounts
    for motif in motifs:
        for i in range(motif_length):
            nucleotide = motif[i]
            numeric_value = numeric_dict[nucleotide]
            pseudocount_matrix[numeric_value][i] += motif_scale

    return pseudocount_matrix

#function to calculate the proabilities of the kmers given a sequence
def estimate_kmer_probabilities(sequence, PWM):
    k = len(PWM[0])
    seq_length = len(sequence)
    seq_to_numeric = 'ACGT0123'
    numeric_dict = {seq_to_numeric[i]: i for i in range(4)}

    probabilities = []
    for i in range(seq_length - k + 1):
        kmer = sequence[i:i+k]
        kmer_probability = 1
        for j in range(k):
            nucleotide = kmer[j]
            numeric_value = numeric_dict[nucleotide]
            kmer_probability *= PWM[numeric_value][j]
        probabilities.append(kmer_probability)
    
    normalization_factor = sum(probabilities)
    probabilities = [p/normalization_factor for p in probabilities]
    
    kmer_probabilities = {sequence[i:i+k]: probabilities[i] for i in range(seq_length - k + 1)}
    
    return kmer_probabilities



# function performs E-M
def update_PWM(sequences, motif_length, PWM):
    k = motif_length
    seq_to_numeric = 'ACGT0123'
    numeric_dict = {seq_to_numeric[i]: i for i in range(4)}
    alpha = 1.0 # pseudocount of 1

    # initialize new PWM with pseudocounts
    new_PWM = [[alpha for j in range(k)] for i in range(4)]

    for seq in sequences:
        # compute probabilities for each kmer in the sequence
        probabilities = estimate_kmer_probabilities(seq, PWM)

        # update PWM based on weighted kmer counts
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            weight = probabilities[kmer]

            # update counts for each position in the kmer
            for j in range(k):
                nucleotide = kmer[j]
                numeric_value = numeric_dict[nucleotide]
                new_PWM[numeric_value][j] += weight

         # normalize counts to obtain probabilities
    # normalize new PWM based on column counts
    for j in range(k):
        total_count = sum([new_PWM[i][j] for i in range(4)])
        for i in range(4):
            new_PWM[i][j] /= total_count    

    return new_PWM

#function to find the reverse complement of a sequence
def ReverseComplement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_seq = [complement_dict[base] for base in reversed(seq)]
    complement = ''.join(complement_seq)
    return complement

#function to calculate proabilities for kmers 
def calculate_probability(kmer, pwm):
    seq1 = 'ACGT0123'
    seq_dict = {seq1[i]: i for i in range(4)}
    p = 1
    for i in range(len(kmer)):
        nucleotide = kmer[i]
        index = seq_dict[nucleotide]
        p *= pwm[index][i]
    return p

#function to keep the highest probability among a kmer in a sequence 
def find_max_probability(sequence, kmer_length, pwm):
    max_prob = 0.0

    for i in range(len(sequence) - kmer_length + 1):
        kmer = sequence[i:i+kmer_length]
        probability = calculate_probability(kmer, pwm)
        if probability > max_prob:
            max_prob = probability

    return max_prob

#function to compute the highest probability for each forward sequence 
def all_forward_probabilites(test_reads, kmer_length, pwm):
    all_forward_probs = []
    
    for read_id, read_seq in test_reads:
        max_prob = find_max_probability(read_seq, kmer_length, pwm)
        read_probs = {'read_id': read_id, 'max_prob': max_prob}
        all_forward_probs.append(read_probs)
    
    return all_forward_probs


#function to compute the highest probability for each reverse sequence 
def all_reverse_probabilites(test_reads, kmer_length, pwm):
    all_rev_probs = []
    
    for read_id, read_seq in test_reads:
        rev_read = ReverseComplement(read_seq)
        max_prob = find_max_probability(rev_read, kmer_length, pwm)
        read_probs = {'read_id': read_id, 'max_prob': max_prob}
        all_rev_probs.append(read_probs)
    
    return all_rev_probs


# Function to merge and write the top 2000 sequences with the highest probabilities to a file
def merge_and_write_top_sequences(forward_probs, reverse_probs, output_file, top_k=2000):
    # Merge forward and reverse probabilities
    merged_probs = forward_probs + reverse_probs

    # Sort the sequences by max_prob in descending order
    merged_probs.sort(key=lambda x: x['max_prob'], reverse=True)

    # Write the top sequences to the output file
    written_seqs = set()  # Track sequences that have been written
    num_written_seqs = 0  # Track the number of written sequences

    with open(output_file, 'w') as file:
        i = 0
        while num_written_seqs < top_k and i < len(merged_probs):
            sequence = merged_probs[i]['read_id']
            if sequence not in written_seqs:
                file.write(f"{sequence}\n")
                written_seqs.add(sequence)
                num_written_seqs += 1
            i += 1


# Access the file paths using argparse
bound_file = args.bound_file
test_file = args.test_file


bound_reads = list(read_fasta_file(bound_file))
# Extract sequences from bound_reads
sequences = [seq for _, seq in bound_reads]


# randomly initialize a PWM
motif_length = 21

motifs = []

# set seed for reproducibility
random.seed(42)

# randomly select a motif from each sequence
motifs = []
for seq in sequences:
  start = random.randint(0, len(seq)-motif_length)
  motif = seq[start:start+motif_length]
  motifs.append(motif)

PWM = compute_pseudocount_profile_matrix(motifs)

#perfrom the EM step multiple times 
num_iterations = 200
for i in range(num_iterations):
    EM_PWM = update_PWM(sequences, motif_length, PWM)
    PWM = EM_PWM

output_file = 'predictions.txt'

test_reads = read_fasta_file(test_file)
for_probs = all_forward_probabilites(test_reads, motif_length, PWM) 

test_reads = read_fasta_file(test_file)
rev_probs = all_reverse_probabilites(test_reads, motif_length, PWM) 

# Output a predictions.txt file of the top 2000 sequences that have the highest probability of binding to a specific transcription factor
merge_and_write_top_sequences(for_probs, rev_probs, output_file, top_k=2000)
