
## SUMAC: supermatrix constructor


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

# Requirements:
    Python 2.7
    Biopython
    MUSCLE
    BLAST+

# To install: 

    python setup.py install

# Usage:

    python -m sumac [-h] [--download_gb DOWNLOAD_GB] [--ingroup INGROUP]
                          [--outgroup OUTGROUP] [--max_outgroup MAX_OUTGROUP]
                          [--evalue EVALUE] [--length LENGTH] [--guide GUIDE]

# Optional arguments:

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

# Examples:

A basic example:

    python -m sumac -d pln -i Onagraceae -o Lythraceae

If you already downloaded the GB database:

    python -m sumac -i Onagraceae -o Lythraceae

Using guide sequences to build clusters:

    python -m sumac -i Onagraceae -o Lythraceae -g guides.fasta

# Citation:
Freyman, W.A. 2014. Supermatrix Constructor (SUMAC): a Python package for data mining GenBank and building phylogenetic supermatrices.


Copyright 2014 Will Freyman - freyman@berkeley.edu  
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
