
# SUMAC: supermatrix constructor


SUMAC (Supermatrix Constructor) is a Python package to data-mine GenBank, construct phylogenetic supermatrices, 
and assess the phylogenetic decisiveness of a matrix given the pattern of missing sequence data. 
SUMAC calculates a novel metric, Missing Sequence Decisiveness Scores (MSDS), which measure how much each 
individual missing sequence contributes to the decisiveness of the matrix. 
MSDS can be used to compare supermatrices and prioritize the acquisition of new sequence data.

SUMAC is designed to be run as a command-line program, though the modules can also be imported and used in other Python scripts. 
SUMAC will assemble supermatrices for any taxonomic group recognized in GenBank, and is optimized to run on multicore systems by utilizing multiple parallel processes.
Please see the [SUMAC user manual](https://rawgit.com/wf8/sumac/master/manual/SUMAC_Manual.pdf) for many more details.

### Requirements:

    Python 2.7
    Biopython
    MAFFT v6.9+
    BLAST+

### To install and use: 

The quick way:

    sudo pip install sumac

For more detailed instructions see the [SUMAC user manual](https://rawgit.com/wf8/sumac/master/manual/SUMAC_Manual.pdf)

### Usuage arguments:

    python -m sumac [-h] [--download_gb DOWNLOAD_GB] [--path PATH]
                   [--ingroup INGROUP] [--outgroup OUTGROUP] [--evalue EVALUE]
                   [--length LENGTH] [--max_ingroup MAX_INGROUP]
                   [--guide GUIDE] [--alignments ALIGNMENTS [ALIGNMENTS ...]]
                   [--salignments SALIGNMENTS [SALIGNMENTS ...]] [--search]
                   [--decisiveness] [--hac]

### Argument details:

    -h, --help            show this help message and exit
    --download_gb DOWNLOAD_GB, -d DOWNLOAD_GB
                          Name of the GenBank division to download (e.g. PLN or
                          MAM).
    --path PATH, -p PATH  Absolute path to download GenBank files to. Defaults
                          to ./genbank/
    --ingroup INGROUP, -i INGROUP
                          Ingroup clade to build supermatrix.
    --outgroup OUTGROUP, -o OUTGROUP
                          Outgroup clade to build supermatrix.
    --cores CORES, -c CORES
                          The number of CPU cores to use for parallel
                          processing. Defaults to the max available.
    --evalue EVALUE, -e EVALUE
                          BLAST E-value threshold to cluster taxa. Defaults to
                          1e-10
    --length LENGTH, -l LENGTH
                          Threshold of sequence length percent similarity to
                          cluster taxa. Defaults to 0.5
    --max_ingroup MAX_INGROUP, -m MAX_INGROUP
                          Maximum number of taxa to include in ingroup. Default
                          is none (no maximum limit).
    --guide GUIDE, -g GUIDE
                          FASTA file containing sequences to guide cluster
                          construction. If this option is selected then all-by-
                          all BLAST comparisons are not performed.
    --alignments ALIGNMENTS [ALIGNMENTS ...], -a ALIGNMENTS [ALIGNMENTS ...]
                          List of aligned FASTA files to build supermatrix
                          instead of mining GenBank.
    --salignments SALIGNMENTS [SALIGNMENTS ...], -sa SALIGNMENTS [SALIGNMENTS ...]
                          List of SUMAC alignments to build supermatrix instead
                          of mining GenBank.
    --search, -s          Turn on search and cluster mode. Will not make
                          alignments or supermatrix.
    --decisiveness, -de   Calculate partial decisiveness. For larger matrices
                          this may be slow.
    --hac                 Use HAC single-linkage clustering algorithm instead of
                            the default SLINK algorithm.

### Citation:

Freyman, W.A. 2015. SUMAC: constructing phylogenetic supermatrices and assessing
partially decisive taxon coverage.

### Other stuff:

Copyright 2015 Will Freyman - freyman@berkeley.edu  
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
