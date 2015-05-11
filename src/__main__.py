#! /usr/bin/python
"""
SUMAC: supermatrix constructor

Copyright 2014 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


import imp
import os
import sys
import argparse
from util import Color
from util import Logger
from genbank import GenBankSetup
from genbank import GenBankSearch
from distancematrix import DistanceMatrixBuilder
from clusters import HACClusterBuilder
from clusters import SLINKClusterBuilder
from clusters import GuidedClusterBuilder
from alignments import Alignments
from supermatrix import Supermatrix



def main():
    # parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--download_gb", "-d", help="Name of the GenBank division to download (e.g. PLN or MAM).")
    parser.add_argument("--path", "-p", help="Absolute path to download GenBank files to. Defaults to ./genbank/")
    parser.add_argument("--ingroup", "-i", help="Ingroup clade to build supermatrix.")
    parser.add_argument("--outgroup", "-o", help="Outgroup clade to build supermatrix.")
    parser.add_argument("--evalue", "-e", help="BLAST E-value threshold to cluster taxa. Defaults to 1e-10")
    parser.add_argument("--length", "-l", help="Threshold of sequence length percent similarity to cluster taxa. Defaults to 0.5")
    parser.add_argument("--max_ingroup", "-m", help="Maximum number of taxa to include in ingroup. Default is none (no maximum limit).") 
    parser.add_argument("--guide", "-g", help="""FASTA file containing sequences to guide cluster construction. If this option is 
                                                 selected then all-by-all BLAST comparisons are not performed.""")
    parser.add_argument("--alignments", "-a", nargs='+', help="List of aligned FASTA files to build supermatrix instead of mining GenBank.")
    parser.add_argument("--salignments", "-sa", nargs='+', help="List of SUMAC alignments to build supermatrix instead of mining GenBank.")
    parser.add_argument("--search", "-s", action='store_true', help="Turn on search and cluster mode. Will not make alignments or supermatrix.")
    parser.add_argument("--decisiveness", "-de", action='store_true', help="Calculate partial decisiveness. For larger matrices this may be slow.")
    parser.add_argument("--hac", action='store_true', help="Use HAC single-linkage clustering algorithm instead of the default SLINK algorithm.")
    args = parser.parse_args()
 
    sys.stdout = Logger()
    color = Color()

    print("")
    print(color.blue + "SUMAC: supermatrix constructor" + color.done)
    print("")

    if args.alignments:
        # if the user provides alignments:
        alignment_files = args.alignments
        alignments = Alignments(alignment_files, "aligned")
    elif args.salignments:
        # if the user inputs SUMAC alignments from previous run
        alignment_files = args.salignments

        # using -sa to input cluster files

        # now align each cluster with MAFFT
        print(color.blue + "Aligning clusters with MAFFT..." + color.done)
        alignments = Alignments(alignment_files, "unaligned")
    
    alignments.print_data()
    alignments.make_gene_region_csv()

    # concatenate alignments
    print(color.purple + "Concatenating alignments..." + color.done)
    supermatrix = Supermatrix(alignments)
   
    try:
        imp.find_module('matplotlib')
        imp.find_module('numpy')
        matplot = True
    except ImportError:
        matplot = False
        print(color.red + "Skipping generating graphs since matplotlib is not installed." + color.done)

    if not args.alignments: # and not args.salignments:
        # only make genbank_csv if the sequences were mined direct from genbank
        supermatrix.make_genbank_csv()
    supermatrix.print_data()
    if matplot:
        supermatrix.make_sequence_data_figure()
    if args.decisiveness:
        supermatrix.print_PD()
        if matplot:
            supermatrix.make_sequence_decisiveness_figure()
        supermatrix.make_decisiveness_csv()
    print(color.yellow + "Final supermatrix: " + color.red + "alignments/supermatrix_concatenated.fasta" + color.done)
    


if __name__ == "__main__":
    main()
