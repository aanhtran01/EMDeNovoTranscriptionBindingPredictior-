EM De Novo Transcription Binding Discovery
==========================================

A Python project that predicts whether a sequence of DNA will be bound by a specific transcription factor in a given condition. 

Given a set of sequences which contains a mix of sequences that are bound by the transcription factor or those that are unbound. The program predicts a subset of the sequences that are likely bound.

Given a given a set of actual DNA sequences where a transcription factor is bound based on a Chromatin Immunoprecipitation (ChIP-seq) experiment, the program using the Expectation-Maximization (EM) algorithm to discover de novo motifs and generate a PWM (positional weight matrix).

The E-M algorithm follows 2 steps:

1. Expection (E) - step: where given a PWM can estimate soft-assignment probabilities for motif instances in a sequence. 

2. Maximization (M) - step: where given given soft-assignment probabilities can easily re-estimate the PWM by taking weight averages. 

The code iterates through the E-M steps 200 times but the number of iterations can be adjusted as necessary in order to reach convergence for different data sets. 


The code outputs the top 2000 sequences with the highest probability of being bound, taking account of reverse sequence complements as well. 

The output of the top 2000 most likely bound sequences will be output in a "predictions.txt" that looks like: 

seq14307
seq1929
seq10381
seq20681
seq10030
seq16200
seq8786
seq2876
seq13921
...

Sequences are ranked from the highest proability of binding to the lowest out of the 2000 top sequences.


Deliverables:
-------------

em.py -- code for predicting bound sequences 

predictions.txt -- output of the top 2000 most likely bound sequences 

predictions.zip -- zipped csv of predictions.txt


Usage
-----

To run the program, navigate to the project directory and run:

> python3 em.py bound.fasta test_reads.fasta 

The program takes the following arguments:

* `--test_reads.fasta`: A fasta file of the test reads that contain bound and unbound sequences.
* `--pwm.txt`: a PWM previously discovered on known sequences.

Examples
--------

Here is an example of how to run the program: 

>python3 em.py bound.fasta test.fasta


Performance
-----------

Runtime is based on laptop computer with the following specifications:
* Processor: 1.1 GHz Quad-Core Intel Core i5
* RAM: 8GB
* Operating System: MacOS Catalina 

For ~25,000 test reads and ~5000 bound reads the run time for motif length of 21 and 200 EM iterations is:

real	23m11.439s
user	23m2.513s
sys	0m1.658s
