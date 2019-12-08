# Uses python3
import sys

# Solving the Imperfect Coverage Problem
'''
We are not guaranteed to be given every possible read at every position,
we cannot expect to use our reads directly and have our de-Bruijn graph have an Eulerian Cycle.
However, our way of combating this is to choose some 𝑘 < ReadLength to break our reads down,
the fragments will be small enough for there to exist an Eulerian Cycle in the graph
'''

'''
Task : Given a list of error-free reads, return an integer 𝑘 such that, when a de Bruijn graph is created from
the 𝑘-length fragments of the reads, the de Bruijn graph has a single possible Eulerian Cycle.

Dataset : The input consist of 400 reads of length 100, each on a separate line. The reads contain no
sequencing errors. Note that you are not given the 100-mer composition of the genome 
(i.e., some 100-mers may be missing).

Output. A single integer 𝑘 on one line
'''

'''
effect on K-mer size
https://www.researchgate.net/figure/Overview-how-to-generate-a-pileup-from-a-read-set-depending-on-the-error-correction_fig3_277405835
https://github.com/rrwick/Bandage/wiki/Effect-of-kmer-size
'''


def is_optimal(k, reads):
    kmers = set()
    for read in reads:
        for i in range(0, len(read) - k + 1):
            kmers.add(read[i:i + k])
    prefixes = set()
    suffixes = set()
    for kmer in kmers:
        prefixes.add(kmer[:-1])
        suffixes.add(kmer[1:])
    return prefixes == suffixes


def read_data():
    input = []
    for i in range(1618):
        read = sys.stdin.readline().strip()
        input.append(read)
    return input


reads = read_data()

for k in range(len(reads[2]), 1, -1):
    if is_optimal(k, reads):
        print(k)
        break
