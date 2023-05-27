# -*- coding: utf-8 -*-
"""
Project goal: build a genome assembler

The input to the project are genome sequencing reads and the output is a predicted genome sequence.

Project builds a spectrum from the sequencing reads and arranges them into a de bruijn graph. The graph is then traversed using DFS. This returns fragments of the predicted reference genome. A simplified overlap-layout-consensus method is used to reconstruct a predicted reference genome. Reads are then aligned onto the reference genome.

the output is a file contains the headers of the reads sorted in the order that they appear in the genome
"""


import random
import argparse
from copy import deepcopy


#set up a command line interface with argparse, create argparse objects
parser = argparse.ArgumentParser(description='Process sequencing reads.')
parser.add_argument('read_file', type=str, help='path to reference file')


args = parser.parse_args()

read_file= args.read_file


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

# Generate a spectrum
def generate_spectrum(reads, k):
    spectrum = {}  # Dictionary to store sequence information
    for read_id, read_seq in reads:
        # Step 1: Break the read into k-mers
        kmers = [read_seq[i:i + k] for i in range(len(read_seq) - k + 1)]
        for kmer in kmers:
            if kmer is not None:
                if kmer in spectrum:
                    spectrum[kmer]["count"] += 1
                    spectrum[kmer]["read_ids"].append(read_id)
                else:
                    spectrum[kmer] = {"count": 1, "read_ids": [read_id]}
    return spectrum

def filter_spectrum(spectrum):
    filtered_spectrum = {kmer: data for kmer, data in spectrum.items() if data["count"] >= 6}
    return filtered_spectrum


#Build a de Bruijn graph:
def de_bruijn(k, patterns):
        adjacency_db = {}
        for p in patterns:
            prefix = p[:k - 1]
            suffix = p[1:]
            adjacency_db.setdefault(prefix, []).append(suffix)
            if suffix not in adjacency_db:
                adjacency_db[suffix] = []
        return adjacency_db
        
def reconstruct_from_path(path):
        return path[0] + ''.join(seq[-1] for seq in path[1:])
        

def dfs(adjacency_db, start_node, visited):
    stack = [start_node]  # Use a stack to keep track of nodes to visit

    while stack:
        node = stack.pop()  # Pop the last node from the stack

        if node not in visited:
            visited.append(node)  # Mark the node as visited

            # Traverse the neighbors (suffixes) of the current node
            for neighbor in adjacency_db.get(node, []):
                stack.append(neighbor)  # Add neighbors to the stack



#BurrowsWheeler  Construction

class suffix_array:
    def __init__(self, text):
        # Constructor that builds the suffix array for a given text string
        self.suff_arr = self._make_suffix_array(text)


    def _sort_array(self, S):
        # Method that sorts the characters of the string S
        l = len(S)
        arrangement = [0] * l
        num = {}
        
        # Count the occurrences of each character in S
        for i in range(l):
            num[S[i]] = num.get(S[i], 0) + 1
        
        # Sort the characters in ascending arrangement
        char_list = sorted(num.keys())
        
        # Compute the starting position of each character in the sorted arrangement
        prev_char = char_list[0]
        for char in char_list[1:]:
            num[char] =  num[char] + num[prev_char]
            prev_char = char
        
        # Compute the arrangement of each suffix based on the starting character
        for i in range(l-1, -1, -1):
            c = S[i]
            num[c] = num[c] - 1
            arrangement[num[c]] = i
        
        return arrangement

    def class_character_arrangement(self, S, arrangement):
        # Method that computes the character classes for each position in the suffix array
        l = len(S)
        class_characters = [0] * l
        class_characters[arrangement[0]] = 0
        
        # Assign the class to each suffix based on whether the starting character is the same as the previous suffix
        for i in range(1, l):
            if S[arrangement[i]] != S[arrangement[i-1]]:
                class_characters[arrangement[i]] = class_characters[arrangement[i-1]] + 1
            else:
                class_characters[arrangement[i]] = class_characters[arrangement[i-1]]
        
        return class_characters

    def _double_suffix_sort(self, S, L, arrangement, class_characters):
        # Method that sorts the doubled suffixes
        string_length = len(S)
        num = [0] * string_length
        new_arrangement = [0] * string_length
        
        # Count the occurrences of each class in the first half of the doubled suffixes
        for i in range(string_length):
            num[class_characters[i]] = num[class_characters[i]] + 1
        
        # Compute the starting position of each class in the sorted arrangement
        for j in range(1, string_length):
            num[j] = num[j] + num[j-1]
        
        # Sort the doubled suffixes based on their second half
        for i in range(string_length-1, -1, -1):
            start = (arrangement[i]-L+string_length) % string_length
            cl = class_characters[start]
            num[cl] = num[cl] - 1
            new_arrangement[num[cl]] = start
        
        return new_arrangement
    
    def _update_classes(self, new_arrangement, class_characters, L):
      # This function updates the character classes based on a new arrangement of indices.
      # new_arrangement: the new arrangement of indices for the string
      # class_characters: a list containing the character classes of the string
      # L: the length of substrings used for sorting
    
      n = len(new_arrangement)
     # n is the length of the new arrangement
    
      class_new = [0] * n
      # Create a new list of character classes with n elements initialized to 0
    
      class_new[new_arrangement[0]] = 0
      # The character class of the first element in the new arrangement is always 0
    
      for i in range(1, n):
          prev = new_arrangement[i-1]
          curr = new_arrangement[i]
          mid = curr + L
          mid_prev = (prev + L) % n
          # Define curr, prev, mid, and mid_prev for easier readability
        
          # Compare the character classes of two adjacent elements in the new arrangement
          # and two adjacent substrings starting at those elements, respectively.
          # If they're different, increment the character class.
          if class_characters[curr] != class_characters[prev] or class_characters[mid] != class_characters[mid_prev]:
              class_new[curr] = class_new[prev] + 1
          else:
              class_new[curr] = class_new[prev]
    
      # Return the updated character classes
      return class_new
    
    def _make_suffix_array(self, S):
      # This function builds the suffix array for a given string S
      # S: the input string
    
      string_length = len(S)
      # The length of S
    
      arrangement = self._sort_array(S)
      # Sort the characters of S and store the resulting arrangement
    
      class_characters = self.class_character_arrangement(S, arrangement)
      # Compute the character classes of S and store them
    
      L = 1
      while L < string_length:
          arrangement = self._double_suffix_sort(S, L, arrangement, class_characters)
          class_characters = self._update_classes(arrangement, class_characters, L)
          L = 2 * L
      # Repeat the process with longer substrings until L >= string_length.
    
      # Return the final suffix array
      return arrangement


    def get_suffix_array(self):
      # This function returns the suffix array of the input string.
    
      return self.suff_arr
      # Return the suffix array stored in self.sa


class BurrowsWheeler:
    def __init__(self, reference_genome):
        self.burrows_wheeler = self.burrows_wheelerFromsuffix_array(reference_genome)
    # Initialize the BurrowsWheeler by calling burrows_wheelerFromsuffix_array with a reference genome
    
    def burrows_wheelerTransform(self, text):
        # This function calculates the Burrows-Wheeler Transform
        # text: the input string
        
        # Generate all possible transfrom of the input string
        transfrom = [text[i:]+text[:i] for i in range(len(text))]
        # Sort the transfrom lexicographically and concatenate their last characters
        burrows_wheeler = ''.join([m[-1] for m in sorted(transfrom)])
        
        return burrows_wheeler

    def burrows_wheelerFromsuffix_array(self, text):
       # This function calculates the Burrows-Wheeler Transform using suffix arrays.
        suff_arr = suffix_array(text).get_suffix_array()
        return ''.join([text[(suff_arr[i]+len(text)-1)%len(text)] for i in range(len(text))])


def HammingDistance(seq1, seq2):
    return len([i for i in range(len(seq1)) if seq1[i] != seq2[i]])


#compute the FirstOccurrence array and CountSymbol function
def FirstOccurrence_CountSymbol(burrows_wheeler, alphabet = ['$','A', 'C', 'G', 'T']):
    l = len(burrows_wheeler)
    CountSymbol = {}
    first_occurances = {}
    for char in alphabet:
        CountSymbol[char] = [0] * (l + 1)
    for i in range(l):
        currChar = burrows_wheeler[i]
        for char, count in CountSymbol.items():
            CountSymbol[char][i+1] = CountSymbol[char][i]
        CountSymbol[currChar][i+1] = CountSymbol[currChar][i+1] + 1
    currIndex = 0
    for char in sorted(alphabet):
        first_occurances[char] = currIndex
        currIndex = currIndex + CountSymbol[char][l]
    return first_occurances, CountSymbol


#perform pattern matching in a Burrows-Wheeler transformed string
#using the better Burrows-Wheeler matching algorithm
def BetterBWMatching(suff_arr, pattern, burrows_wheeler, starts, counts):
    occs = set()
    top = 0
    bottom = len(burrows_wheeler) - 1
    currIndex = len(pattern) - 1
    while top <= bottom:
        if currIndex >= 0:
            symbol = pattern[currIndex]
            currIndex = currIndex - 1
            if counts[symbol][bottom+1] - counts[symbol][top] > 0:
                top = starts[symbol] + counts[symbol][top]
                bottom = starts[symbol] + counts[symbol][bottom+1] - 1
            else:
                break
        else:
            for i in range(top, bottom + 1):
                occs.add(suff_arr[i])
            break

    return occs


#function to align a single read to a genome
def align_read_to_genome(read, reference_genome, k, suffix_array, burrows_wheeler, starts, counts):

    # Step 1: Break the read into k-mers
    kmers = [read[i:i+k] for i in range(len(read)-k+1)]
    #first_kmer = kmers[0]


    # Step 2: Search for matches in the BurrowsWheeler  index and extend the alignment
    best_match = None
    best_score = float('inf')
    for i, kmer in enumerate(kmers):
        positions = BetterBWMatching(suffix_array, kmer, burrows_wheeler, starts, counts)

        for pos in positions:
            # extend the alignment
            offset = i
            alignment_start = pos - offset
            alignment_end = alignment_start + len(read)
            if alignment_start < 0 or alignment_end > len(reference_genome):
                continue # alignment out of bounds, skip to next position
            ref_sequence = reference_genome[alignment_start:alignment_end]
            score = HammingDistance(read, ref_sequence)

            # check if this is the best match so far
            if score < best_score:
                best_score = score
                best_match = alignment_start

    return best_match, best_score

#function to align all reads to the genome
def align_all_reads_to_genome(suffix_array, donor_reads, reference_genome, k, burrows_wheeler, starts, counts):
    results = []
    for read_id, read_seq in donor_reads:
        best_match, best_score = align_read_to_genome(read_seq, reference_genome, k, suffix_array, burrows_wheeler, starts, counts)
        results.append({'donor_read_id': read_id,'sequence' : read_seq ,'best_match': best_match, 'best_score': best_score})
    return results


k_spectrum = 20
reads = read_fasta_file(read_file)
spectrum = generate_spectrum(reads, k_spectrum)
filtered_spectrum = filter_spectrum(spectrum)

de_bruijn = de_bruijn(k_spectrum, filtered_spectrum)
#de_bruijn(k, filtered_spectrum)

# Set the seed for the random number generator
random.seed(123)

visited_dict = {}
fragments = []

# Get a random sample of 1000 keys from de_bruijn dictionary
sample_keys = int(len(de_bruijn.keys())*0.05)
#random_keys = random.sample(de_bruijn.keys(), 100)
random_keys = random.sample(list(de_bruijn.keys()), sample_keys)

sequence_lengths = []  # List to store sequence lengths

for key in random_keys:
    visited = []
    dfs(de_bruijn, key, visited)
    visited_dict[key] = visited

    reconstructed_sequence = reconstruct_from_path(visited)
    sequence_length = len(reconstructed_sequence)

    if sequence_length >= 0:
        fragments.append(reconstructed_sequence)
        sequence_lengths.append(sequence_length)


# Filter fragments with less than desired threshold length
#threshold_length = 18250
#threshold_length = 850
#filtered_fragments = [fragment for fragment in fragments if len(fragment) >= threshold_length]

# Calculate the lengths of the filtered fragments
#fragment_lengths = [len(fragment) for fragment in filtered_fragments]

# Sort the fragments based on their size and only use the top longest sequences
sorted_fragments = sorted(fragments, key=len, reverse=True)
top_longest = sorted_fragments[:1]

# Store the longest sequence as a string
longest_sequence = ''.join(top_longest)



reference_genome = longest_sequence
reference_genome += "$"

#calling the functions
burrows_wheeler = BurrowsWheeler(reference_genome).burrows_wheeler
starts, counts = FirstOccurrence_CountSymbol(burrows_wheeler)
suffix_array = suffix_array(reference_genome).get_suffix_array()

k = 16

#filename = "/content/project3b_20000_reads_without_positions.fasta"

donor_reads = read_fasta_file(read_file)

results = align_all_reads_to_genome(suffix_array, donor_reads, reference_genome, k, burrows_wheeler, starts, counts)
#print(results)

# Filter out results with None as the best_match value
filtered_results = [result for result in results if result['best_match'] is not None]
#print(filtered_results)

# Sort the filtered results based on the 'best_match' value
sorted_results = sorted(filtered_results, key=lambda x: x['best_match'])


# Open a text file for writing
with open('predictions.txt', 'w') as file:
    for result in sorted_results:
        file.write(f'>{result["donor_read_id"]}\n')
