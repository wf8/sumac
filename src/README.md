
# SUMAC 2.22: supermatrix constructor 


SUMAC (Supermatrix Constructor) is a Python package to data-mine GenBank, construct phylogenetic supermatrices, 
and assess the phylogenetic decisiveness of a matrix given the pattern of missing sequence data. 
SUMAC calculates a novel metric, Missing Sequence Decisiveness Scores (MSDS), which measure how much each 
individual missing sequence contributes to the decisiveness of the matrix. 
MSDS can be used to compare supermatrices and prioritize the acquisition of new sequence data.

SUMAC is designed to be run as a command-line program, though the modules can also be imported and used in other Python scripts. 
SUMAC will assemble supermatrices for any taxonomic group recognized in GenBank, and is optimized to run on multicore systems by utilizing multiple parallel processes.
Please see the [SUMAC user manual](https://rawgit.com/wf8/sumac/master/manual/SUMAC_Manual.pdf) for many more details.

SUMAC works on OSX and Linux. Windows is not currently supported.

:exclamation: [USEARCH bug 9.15](http://www.drive5.com/usearch/manual/bugs.html) will cause USEARCH to crash during a SUMAC run. 
This bug affects USEARCH versions v9.2.64 and v9.0.2132. Please ensure you use a different version of USEARCH with SUMAC. Until a new version is released use v8.1 which you can download [here](http://www.drive5.com/usearch/download.html).

### Updates:

- SUMAC v2.22 updates SUMAC to handle SeqRecord comparisons deprecated in Biopython 1.67.
- SUMAC v2.21 handles a bug in Biopython's GenBank parsing that could cause SUMAC to crash.
- SUMAC v2.2 allows multiple taxonomic groups to be searched for (not just a single ingroup and outgroup). See manual for details.
- SUMAC v2.1 greatly improves the handling of taxonomic names, adds the `--min_clusters` flag to set the minimum number of taxa to include in a cluster, 
and other additional error checks.
- SUMAC v2.0 is *significantly* faster than previous versions due to a new clustering algorithm.  

### Citation:

Freyman, W.A. 2015. SUMAC: constructing phylogenetic supermatrices and assessing
partially decisive taxon coverage. *Evolutionary Bioinformatics* 2015:11 263-266

### Other stuff:

Copyright 2016 Will Freyman - freyman@berkeley.edu  
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

### Requirements:

    Python 2.7
    Biopython
    MAFFT v6.9+
    USEARCH
    BLAST+ (optional)

### To install and use: 

The quick way:

    sudo pip install sumac

For more detailed instructions see the [SUMAC user manual](https://rawgit.com/wf8/sumac/master/manual/SUMAC_Manual.pdf)

### Usage:

    python -m sumac [-h] [--download_gb DOWNLOAD_GB]
                   [--download_gb2 DOWNLOAD_GB2] [--path PATH]
                   [--ingroup INGROUP] [--outgroup OUTGROUP] [--cores CORES]
                   [--id ID] [--evalue EVALUE] [--length LENGTH]
                   [--maxlength MAXLENGTH] [--minlength MINLENGTH]
                   [--min_clusters MIN_CLUSTERS] [--max_ingroup MAX_INGROUP] [--guide GUIDE]
                   [--alignments ALIGNMENTS [ALIGNMENTS ...]]
                   [--salignments SALIGNMENTS [SALIGNMENTS ...]] [--search]
                   [--decisiveness] [--hac] [--slink]

### Argument details:

    -h, --help            show this help message and exit                                                                                                                      [8/1029]
    --download_gb DOWNLOAD_GB, -d DOWNLOAD_GB
                          Name of the GenBank division to download (e.g. PLN or
                          MAM).
    --download_gb2 DOWNLOAD_GB2, -d2 DOWNLOAD_GB2
                          Name of the optional second GenBank division to
                          download. Use this option if the ingroup and outgroup
                          are in different GenBank divisions.
    --path PATH, -p PATH  Absolute path to download GenBank files to. Defaults
                          to ./genbank/
    --ingroup INGROUP, -i INGROUP
                          Ingroup clade to build supermatrix.
    --outgroup OUTGROUP, -o OUTGROUP
                          Outgroup clade to build supermatrix.
    --cores CORES, -c CORES
                          The number of CPU cores to use for parallel
                          processing. Defaults to the max available.
    --id ID, -id ID       UCLUST id threshold to cluster taxa. Defaults to 0.50
    --evalue EVALUE, -e EVALUE
                          BLAST E-value threshold to cluster taxa. Defaults to
                          1e-10
    --length LENGTH, -l LENGTH
                          Threshold of sequence length percent similarity to
                          cluster taxa. Defaults to 0.25
    --maxlength MAXLENGTH, -maxl MAXLENGTH
                          Maximum length of sequences to include in UCLUST
                          clusters. Defaults to 5000
    --minlength MINLENGTH, -minl MINLENGTH
                          Minimum length of sequences to include in UCLUST
                          clusters. Defaults to 100
    --min_clusters MIN_CLUSTERS, -minc MIN_CLUSTERS
                          Minimum number of taxa needed for clusters. Defaults
                          to 4
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
                          the default UCLUST algorithm.
    --slink               Use the SLINK clustering algorithm instead of the
                          default UCLUST algorithm.
                          
