#! /usr/bin/python
"""
supermatrix.py

Will Freyman
freyman@berkeley.edu

Python module that:
1: Downloads GenBank database of the specified GB division (PLN, MAM, etc)
2: Build clusters of all sequences for ingroup and outgroup clades
	- use hierarchical agglomerative clustering algorithm
	- distance threshold default blastn e-value 1.0e-10 
3: Determines whether sequences must be reverse complemented.
4: Aligns each cluster of sequences using MUSCLE.
5: Concatenates the clusters creating a supermatrix.
6: Prunes out excess members of the outgroup, keeping those with the
   most sequence data.

Requirements:
    Python 2.7
    Biopython
    MUSCLE
    BLAST+

usage: supermatrix.py [-h] [--download_gb DOWNLOAD_GB] [--ingroup INGROUP]
                      [--outgroup OUTGROUP] [--max_outgroup MAX_OUTGROUP]
                      [--evalue EVALUE]

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


Example: 
supermatrix.py -d pln -i Onagraceae -o Lythraceae

Example for testing:
python supermatrix.py -i Ceratobasidiaceae -o ceratobasidiaceae


"""

import os
import sys
import gzip
import argparse
from ftplib import FTP
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

class color:
    purple = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    done = '\033[0m'

    def disable(self):
        self.purple = ''
        self.blue = ''
        self.green = ''
        self.yellow = ''
        self.red = ''
        self.done = ''


def gettext(ftp, filename, outfile=None):
    """ 
    Fetch a text file. 
    """
    if outfile is None:
        outfile = sys.stdout
    # use a lambda to add newlines to the lines read from the server
    ftp.retrlines("RETR " + filename, lambda s, w=outfile.write: w(s+"\n"))

def getbinary(ftp, filename, outfile=None):
    """ 
    Fetch a binary file. 
    """
    if outfile is None:
        outfile = sys.stdout
    ftp.retrbinary("RETR " + filename, outfile.write)

def download_gb_db(division):
    """
    Downloads and uncompresses files for a GenBank division.
    """
    print(color.yellow + "Connecting to ftp.ncbi.nih.gov..." + color.done)
    ftp = FTP("ftp.ncbi.nih.gov")
    ftp.login()
    print(color.yellow + "Opening directory genbank..." + color.done)
    ftp.cwd("genbank")
    file_list = ftp.nlst()
    i = 1 
    file_name = "gb" + division + str(i) + ".seq.gz"
    if not os.path.exists("genbank"):
        os.makedirs("genbank")
    while file_name in file_list:
        print(color.red + "Downloading file " + file_name + color.done)
        file = open("genbank/" + file_name, "wb")
        getbinary(ftp, file_name, file)
	file.close()
	print(color.yellow + "Uncompressing file " + file_name + color.done)
	file = gzip.open("genbank/" + file_name, "rb")
	file_content = file.read()
	file.close()
        file = open("genbank/" + file_name[:-3], "wb")
	file.write(file_content)
	file.close()
	os.remove("genbank/" + file_name)
	i += 1
	file_name = 'gb' + division + str(i) + '.seq.gz'
    ftp.quit()

def setup_sqlite():
    """
    Sets up the SQLite db for the GenBank division.
    Returns a dictionary of SeqRecord objects.
    """
    print(color.yellow + "Setting up SQLite database..." + color.done)
    os.chdir("genbank/")
    files = os.listdir(".")
    print files
    gb = SeqIO.index_db("gb.idx", files, "genbank")
    os.chdir("..")
    return gb

def get_in_out_groups(gb, ingroup, outgroup):
    """
    Takes as input a dictionary of SeqRecords gb and the names of ingroup and outgroup clades.
    Returns lists of keys to SeqRecords for the ingroup and outgroup.
    """
    os.chdir("genbank/")
    keys = gb.keys()
    ingroup_keys = []
    outgroup_keys = []
    i = 0          # FOR TESTING ONLY
    for key in keys:
	if ingroup in gb[key].annotations['taxonomy']:
            ingroup_keys.append(key)
	elif outgroup in gb[key].annotations['taxonomy']:
	    outgroup_keys.append(key)
	if i == 200:                    #
	    return ingroup_keys, outgroup_keys    #
	else:                            # FOR TESTING ONLY
	    i += 1                      #
    return ingroup_keys, outgroup_keys

def make_distance_matrix(gb, seq_keys):
    """
    Takes as input a dictionary of SeqRecords gb and the keys to all sequences.
    Returns a 2 dimensional list of distances. Distances are blastn e-values.
    """
    dist_matrix = [[99] * len(seq_keys) for i in range(len(seq_keys))]
    i = 0
    j = 0
    for key in seq_keys:
        output_handle = open('subject.fasta', 'w')
        SeqIO.write(gb[key], output_handle, 'fasta')
        output_handle.close()

        for key2 in seq_keys:
	    # only calculate e-values for pairs that have not yet been compared
	    if dist_matrix[i][j] == 99:
	        if key == key2:
	            dist_matrix[i][j] = 0.0
	        else:
                    output_handle = open('query.fasta', 'w')
                    SeqIO.write(gb[key2], output_handle, 'fasta')
                    output_handle.close()

                    blastn_cmd = NcbiblastnCommandline(query='query.fasta', subject='subject.fasta', out='blast.xml', outfmt=5)
                    stdout, stderr = blastn_cmd()
                    blastn_xml = open('blast.xml', 'r')
                    blast_records = NCBIXML.parse(blastn_xml)

                    for blast_record in blast_records:
			if blast_record.alignments:
			    if blast_record.alignments[0].hsps:
                                # blast hit found
			        dist_matrix[i][j] = blast_record.alignments[0].hsps[0].expect
                                dist_matrix[j][i] = blast_record.alignments[0].hsps[0].expect
		        else:
			    # no blast hit found, set distance to default 10.0
		            dist_matrix[i][j] = 10.0
			    dist_matrix[j][i] = 10.0
		    blastn_xml.close()
	    j += 1
        j = 0
	i += 1
	sys.stdout.write(color.blue + '.' + color.done)    
	sys.stdout.flush()    
    sys.stdout.write("\n")
    sys.stdout.flush()
    return dist_matrix

def make_clusters(seq_keys, distance_matrix, e=0.1):
    """
    Input: seq_keys a list of all sequences used in the analysis, distance_matrix based on BLAST e-values, and an optional e-value threshold for clustering.
    Output: a list of clusters (each cluster is itself a list of keys to sequences)
    """
    i = 0
    j = 0
    clusters = []
    for key in seq_keys:
        # skip this sequence if it is already in a cluster
        in_cluster = False
        for cluster in clusters:
	    if key in cluster:
	        in_cluster = True
		break
	if not in_cluster:
	    # it is not in a cluster, so make a new cluster
	    new_cluster = [key]
	    # check to see if other sequences should be in new cluster
	    j = i + 1
            while j < len(seq_keys):
	        key2 = seq_keys[j]
		# skip the other sequences that are already in a cluster
	        key2_in_cluster = False
	        for cluster in clusters:
		    if key2 in cluster:
		        key2_in_cluster = True
			break
		# if the other sequence is not in a cluster and it is below the e-value threshold
		# add it to the new cluster
		if not key2_in_cluster and distance_matrix[i][j] <= e:
                    new_cluster.append(key2)
		j += 1
	    clusters.append(new_cluster)
        i += 1
    return clusters	

def make_matrices(gb, clusters):
    """
    Inputs the dictionary of all GenBank sequences and a list of clustered accessions.
    Must determine whether sequences must be reverse complemented by using BLAST.
    Outputs a list of FASTA files, each file containing an unaligned sequence cluster.
    """
    matrices = []
    # do stuff
    return matrices

def align_matrices(matrices):
    """
    Inputs a list of FASTA files, each file containing an unaligned sequence cluster.
    Outputs a list of FASTA files, each file containing an aligned sequence cluster.
    """
    # do stuff
    # return aligned_matrices

if __name__ == "__main__":
    # parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--download_gb", "-d", help="Name of the GenBank division to download (e.g. PLN or MAM).")
    parser.add_argument("--ingroup", "-i", help="Ingroup clade to build supermatrix.")
    parser.add_argument("--outgroup", "-o", help="Outgroup clade to build supermatrix.")
    parser.add_argument("--max_outgroup", "-m", help="Maximum number of taxa to include in outgroup. Defaults to 10.")
    parser.add_argument("--evalue", "-e", help="BLAST E-value threshold to cluster taxa. Defaults to 0.1")
    args = parser.parse_args()
    
    if args.download_gb:
        gb_division = str(args.download_gb).lower()
        download_gb_db(gb_division)
        # set up SQLite database
        gb = setup_sqlite()

    elif not os.path.exists("genbank/gb.idx"):
        print(color.red + "GenBank database not downloaded. Re-run with the -d option. See --help for more details." + color.done)
        sys.exit(0)
    else:
        print(color.purple + "Genbank database already downloaded. Indexing sequences..." + color.done)
	gb = SeqIO.index_db("genbank/gb.idx")

    print(color.purple + "%i sequences indexed!" % len(gb) + color.done)

    if args.ingroup and args.outgroup:
        ingroup = args.ingroup
        outgroup = args.outgroup
        print(color.blue + "Outgroup = " + outgroup)
	print("Ingroup = " + ingroup + color.done)
        print(color.blue + "Searching for ingroup and outgroup sequences..." + color.done)
        ingroup_keys, outgroup_keys = get_in_out_groups(gb, ingroup, outgroup)
        
	print(color.yellow + "Found " + str(len(ingroup_keys)) + " ingroup sequences." + color.done)
        print(color.yellow + "Found " + str(len(outgroup_keys)) + " outgroup sequences." + color.done)
        
	all_seq_keys = ingroup_keys + outgroup_keys

        print(color.blue + "Making distance matrix for all sequences..." + color.done)
	distance_matrix = make_distance_matrix(gb, all_seq_keys)

        print(color.purple + "Clustering sequences..." + color.done)
	if args.evalue:
	    print(color.purple + "Using e-value threshold " + str(args.evalue) + color.done) 
	    clusters = make_clusters(gb, distance_matrix, float(args.evalue))
	else:
	    print(color.purple + "Using default e-value threshold 0.1" + color.done)
            clusters = make_clusters(all_seq_keys, distance_matrix)
            
	print(color.yellow + "Building sequence matrices for each cluster." + color.done)
	print(color.yellow + "Checking for sequences that must be reverse complemented..." + color.done)
	unaligned_matrices = make_matrices(gb, clusters)


        # now align each cluster, then concantenate, then reduce the number of outgroup taxa
