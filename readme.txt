Genome Assembler
================

A Python project that assembles a reference genome from given reads. A spectrum of k is generated from the reads. A de bruijn graph can be built from the spectrum and traversed. However, due to sequencing errors in the reads, we cannot traverse the graph using a singular Eulerian path. However, this algorithm bypasses this by finding the long enough path to map most of the genome. Taking a sample of the start nodes in the graph (aroung 5% of the start nodes), the program will traverse the graph for each node using DFS (Depth First Search). The longest path will be used as the path for the predicted reference genome and the reads will be aligned to this reference genome. 

The input would be a spectrum and the output a is a predicted genome sequence. The output contains the headers of the spectrum in the sorted order that they appear in the genome. This is done by aligning the spectrum back to the predicted assmbled sequence.

Deliverables:
-------------

genomeassembler.py -- code for genome assembly 

predictions.txt -- output of the headers of the spectrum in the sorted order that they appear in the genome

predictions.zip -- zipped csv of predictions.txt


Usage
-----
The program takes in one input, the reads fasta without the genome positions 

To run the program, navigate to the project directory and run:

> python3 euleriangenomeassembler.py reads.fasta

The program takes the following argument:

* `--reads.fasta`: A fasta file of the sequence reads without genome positions 

Examples
--------

Here is an example of how to run the program: 

>python3 genomeassembler.py project3b_20000_reads_without_positions.fasta

Performance
-----------

Runtime is based on laptop computer with the following specifications:
* Processor: 1.1 GHz Quad-Core Intel Core i5
* RAM: 8GB
* Operating System: MacOS Catalina 

For assembling a reference genome of 20000 nucleotides and spectrum alignment is for a given fasta of ~23999 reads is:

real	50m41.759s
user	50m30.366s
sys	0m5.352s

Future Improvements
-------------------
It would be ideal to traverse all the start nodes to find the longest path, but due to the high volume number of reads, this would take quite a while and be impractical. 
