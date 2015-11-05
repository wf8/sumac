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
import multiprocessing
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
    parser.add_argument("--download_gb2", "-d2", help="""Name of the optional second GenBank division to download. Use this
                                                         option if the ingroup and outgroup are in different GenBank divisions.""")
    parser.add_argument("--path", "-p", help="Absolute path to download GenBank files to. Defaults to ./genbank/")
    parser.add_argument("--ingroup", "-i", help="Ingroup clade to build supermatrix.")
    parser.add_argument("--outgroup", "-o", help="Outgroup clade to build supermatrix.")
    parser.add_argument("--cores", "-c", help="The number of CPU cores to use for parallel processing. Defaults to the max available.")
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

    num_cores = multiprocessing.cpu_count()
    if args.cores and int(args.cores) <= num_cores:
        num_cores = int(args.cores) 

    if args.alignments:
        # if the user provides alignments:
        alignment_files = args.alignments
        alignments = Alignments(alignment_files, "aligned", num_cores)
    elif args.salignments:
        # if the user inputs SUMAC alignments from previous run
        alignment_files = args.salignments
        alignments = Alignments(alignment_files, "sumac_aligned", num_cores)
    else:
        if args.search:
            print(color.yellow + "Running in search and cluster mode. Clusters will not be aligned and supermatrix will not assembled." + color.done) 

        # first download and set up sqllite db if necessary
        if args.path:
            gb_dir = args.path
        else:
            gb_dir = os.path.abspath("genbank/")
        # if the user requests downloading
        if args.download_gb:
            divisions = [args.download_gb]
            if args.download_gb2:
                divisions.append(args.download_gb2)
            GenBankSetup.download(divisions, gb_dir)
            print(color.yellow + "Setting up SQLite database..." + color.done)
            gb = GenBankSetup.sqlite(gb_dir)
        # the user didn't request downloading, so check for genbank directory
        elif not os.path.exists(gb_dir):
            print(color.red + "GenBank database not downloaded. Re-run with the -d option. See --help for more details." + color.done)
            sys.exit(0)
        # the genbank directory exists so check for sequences and index them
        else:
            gb = GenBankSetup.sqlite(gb_dir)
        print(color.purple + "%i sequences indexed!" % len(gb) + color.done)

        # check for ingroup and outgroup
        if args.ingroup:
            ingroup = args.ingroup
            if args.outgroup:
                outgroup = args.outgroup
            else:
                outgroup = "NONE"
        else:
            print(color.red + "Please specify ingroup. See --help for details." + color.done)
            sys.exit(0)
        
        # search db for ingroup and outgroup sequences
        print(color.blue + "Ingroup = " + ingroup + color.done)
        if args.outgroup:
            print(color.blue + "Outgroup = " + outgroup + color.done)
        print(color.blue + "Searching for ingroup and outgroup sequences..." + color.done)
        if args.max_ingroup:
            search_results = GenBankSearch(gb, ingroup, outgroup, int(args.max_ingroup))
        else:
            search_results = GenBankSearch(gb, ingroup, outgroup)
        ingroup_keys = search_results.ingroup_keys
        outgroup_keys = search_results.outgroup_keys
        all_seq_keys = ingroup_keys + outgroup_keys
        if len(all_seq_keys) == 0:
            print(color.red + "No sequences found for the ingroup and outgroup!" + color.done)
            sys.exit(0)

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
            cluster_builder = GuidedClusterBuilder(args.guide, all_seq_keys, length_threshold, evalue_threshold, gb_dir, num_cores)
        else:
            # make distance matrix
            print(color.blue + "Making distance matrix for all sequences..." + color.done)
            distance_matrix = DistanceMatrixBuilder(gb, all_seq_keys, length_threshold, gb_dir, num_cores).distance_matrix

            # cluster sequences
            if args.hac:
                print(color.purple + "Clustering sequences using the HAC algorithm..." + color.done)
                cluster_builder = HACClusterBuilder(all_seq_keys, distance_matrix, evalue_threshold)
            else:
                print(color.purple + "Clustering sequences using the SLINK algorithm..." + color.done)
                cluster_builder = SLINKClusterBuilder(all_seq_keys, distance_matrix, evalue_threshold)

        print(color.purple + "Found " + color.red + str(len(cluster_builder.clusters)) + color.purple + " clusters." + color.done)
        if len(cluster_builder.clusters) == 0:
            print(color.red + "No clusters found." + color.done)
            sys.exit(0)

        # filter clusters, make FASTA files
        print(color.yellow + "Building sequence matrices for each cluster." + color.done)
        cluster_builder.assemble_fasta(gb)
        print(color.purple + "Kept " + color.red + str(len(cluster_builder.clusters)) + color.purple + " clusters, discarded those with < 4 taxa." + color.done)
        
        # if we are in search and cluster mode we are done
        if args.search:
            sys.exit(0) 
        
        if len(cluster_builder.clusters) == 0:
            print(color.red + "No clusters left to align." + color.done)
            sys.exit(0)
        # now align each cluster with MUSCLE
        print(color.blue + "Aligning clusters with MAFFT..." + color.done)
        alignments = Alignments(cluster_builder.cluster_files, "unaligned", num_cores)
    
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
