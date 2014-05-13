#! /usr/bin/python
"""

SUMAC: supermatrix constructor

Will Freyman
freyman@berkeley.edu

Freyman, W.A. 2014. Supermatrix Constructor (SUMAC): a Python package for data mining GenBank and building phylogenetic supermatrices.

Python package that:
1: Downloads GenBank database of the specified GB division (PLN, MAM, etc)
2: Perform exhaustive all-by-all BLAST comparisons of each ingroup and outgroup sequence.
3: Alternatively (much faster) uses a FASTA file of guide sequences to define each cluster. 
   Each ingroup/outgroup sequence is BLASTed against the guide sequences.
4: Build clusters of sequences:
    - use single-linkage hierarchical clustering algorithm
    - distance threshold default BLASTn e-value 1.0e-10 
        - uses sequence length percent similarity cutoff default 0.5
    - discards clusters that are not phylogenetically informative (< 4 taxa)
5: Aligns each cluster of sequences using MUSCLE.
6: Concatenates the clusters creating a supermatrix.

Optimized to run on multicore servers with the use of parallel multiple processes.

Requirements:
    Python 2.7
    Biopython
    MUSCLE
    BLAST+

install: python setup.py install

usage: python -m sumac [-h] [--download_gb DOWNLOAD_GB] [--ingroup INGROUP]
                      [--outgroup OUTGROUP] [--max_outgroup MAX_OUTGROUP]
                      [--evalue EVALUE] [--length LENGTH] [--guide GUIDE]

optional arguments:
  -h, --help            show this help message and exit
  --download_gb DOWNLOAD_GB, -d DOWNLOAD_GB
                        Name of the GenBank division to download (e.g. PLN or
                        MAM).
  --ingroup INGROUP, -i INGROUP
                        Ingroup clade to build supermatrix.
  --outgroup OUTGROUP, -o OUTGROUP
                        Outgroup clade to build supermatrix.
  --max_outgroup MAX_OUTGROUP, -m MAX_OUTGROUP
                        Maximum number of taxa to include in outgroup.
                        Defaults to 10.
  --evalue EVALUE, -e EVALUE
                        BLAST E-value threshold to cluster taxa. Defaults to
                        0.1
  --length LENGTH, -l LENGTH
                        Threshold of sequence length percent similarity to 
            cluster taxa. Defaults to 0.5
  --guide GUIDE, -g GUIDE
                        FASTA file containing sequences to guide cluster 
            construction. If this option is selected then 
            all-by-all BLAST comparisons are not performed.


Example: 
python -m sumac -d pln -i Onagraceae -o Lythraceae

If you already downloaded the GB database:
python -m sumac -i Onagraceae -o Lythraceae

Using guide sequences to cluster:
python -m sumac -i Onagraceae -o Lythraceae -g guides.fasta

"""

import os
import sys
import argparse

from util import Color
from genbank import GenBankSetup
from genbank import GenBankSearch
import distancematrix
import clusters
import alignments


def main():
    # parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--download_gb", "-d", help="Name of the GenBank division to download (e.g. PLN or MAM).")
    parser.add_argument("--ingroup", "-i", help="Ingroup clade to build supermatrix.")
    parser.add_argument("--outgroup", "-o", help="Outgroup clade to build supermatrix.")
    parser.add_argument("--max_outgroup", "-m", help="Maximum number of taxa to include in outgroup. Defaults to 10.")
    parser.add_argument("--evalue", "-e", help="BLAST E-value threshold to cluster taxa. Defaults to 0.1")
    parser.add_argument("--length", "-l", help="Threshold of sequence length percent similarity to cluster taxa. Defaults to 0.5")
    parser.add_argument("--guide", "-g", help="""FASTA file containing sequences to guide cluster construction. If this option is 
                                                 selected then all-by-all BLAST comparisons are not performed.""")
    args = parser.parse_args()
 
    color = Color()

    print("")
    print(color.blue + "SUMAC: supermatrix constructor" + color.done)
    print("")

    # download and set up sqllite db
    if args.download_gb:
        GenBankSetup.download(args.download_gb)
        print(color.yellow + "Setting up SQLite database..." + color.done)
        gb = GenBankSetup.sqlite()
    elif not os.path.exists("genbank/gb.idx"):
        print(color.red + "GenBank database not downloaded. Re-run with the -d option. See --help for more details." + color.done)
        sys.exit(0)
    else:
        print(color.purple + "Genbank database already downloaded. Indexing sequences..." + color.done)
    gb_dir = os.path.abspath("genbank/")
    gb = SeqIO.index_db(gb_dir + "/gb.idx")
    print(color.purple + "%i sequences indexed!" % len(gb) + color.done)

    # check for ingroup and outgroup
    if args.ingroup and args.outgroup:
        ingroup = args.ingroup
        outgroup = args.outgroup
    else:
        print(color.red + "Please specify ingroup and outgroup. See --help for details." + color.done)
    sys.exit(0)
    
    # search db for sequences
    print(color.blue + "Outgroup = " + outgroup)
    print("Ingroup = " + ingroup + color.done)
    print(color.blue + "Searching for ingroup and outgroup sequences..." + color.done)
    search_results = GenBankSearch(gb, ingroup, outgroup)
    ingroup_keys = search_results.ingroup_keys
    outgroup_keys = search_results.outgroup_keys
    all_seq_keys = ingroup_keys + outgroup_keys

    # determine sequence length similarity threshold
    length_threshold = 0.5
    if args.length:
        length_threshold = args.length
    print(color.blue + "Using sequence length similarity threshold " + color.red + str(length_threshold) + color.done)

    # determine e-value threshold
    evalue_threshold = (1.0/10**10)
    if args.evalue:
        evalue_threshold = float(args.evalue)
    print(color.blue + "Using BLAST e-value threshold " + color.red + str(evalue_threshold) + color.done)

    # now build clusters, first checking whether we are using FASTA file of guide sequences
    # or doing all-by-all comparisons
    if args.guide:
        # use FASTA file of guide sequences
        print(color.blue + "Building clusters using the guide sequences..." + color.done)
        clusters = clusters.make_guided_clusters(args.guide, all_seq_keys, length_threshold, evalue_threshold)
    else:
        # make distance matrix
        print(color.blue + "Making distance matrix for all sequences..." + color.done)
        distance_matrix = distancematrix.make_distance_matrix(gb, all_seq_keys, length_threshold)

        # cluster sequences
        print(color.purple + "Clustering sequences..." + color.done)
        clusters = clusters.make_clusters(all_seq_keys, distance_matrix, evalue_threshold)
    
    print(color.purple + "Found " + color.red + str(len(clusters)) + color.purple + " clusters." + color.done)
    if len(clusters) == 0:
        print(color.red + "No clusters found." + color.done)
        sys.exit(0)

    # filter clusters, make FASTA files
    print(color.yellow + "Building sequence matrices for each cluster." + color.done)
    cluster_files, clusters = alignments.assemble_fasta_clusters(gb, clusters)
    print(color.purple + "Kept " + color.red + str(len(clusters)) + color.purple + " clusters, discarded those with < 4 taxa." + color.done)
    if len(clusters) == 0:
        print(color.red + "No clusters left to align." + color.done)
    sys.exit(0)

    # now align each cluster with MUSCLE
    print(color.blue + "Aligning clusters with MUSCLE..." + color.done)
    alignment_files = alignments.align_clusters(cluster_files)
    print_region_data(alignment_files)

    # concatenate alignments
    print(color.purple + "Concatenating alignments..." + color.done)
    final_alignment = alignments.concatenate(alignment_files)
    print(color.yellow + "Final alignment: " + color.red + "alignments/final.fasta" + color.done)

    # TODO:
    # reduce the number of outgroup taxa, make graphs, etc



if __name__ == "__main__":
    main()
