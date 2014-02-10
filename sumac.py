#! /usr/bin/python
"""

SUMAC: supermatrix constructor
sumac.py

Will Freyman
freyman@berkeley.edu

Python module that:
1: Downloads GenBank database of the specified GB division (PLN, MAM, etc)
2: Build clusters of all sequences for ingroup and outgroup clades
	- use single-linkage hierarchical clustering algorithm
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

usage: sumac.py [-h] [--download_gb DOWNLOAD_GB] [--ingroup INGROUP]
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
  --length LENGTH, -; LENGTH
                        Threshold of sequence length percent similarity to 
			cluster taxa. Defaults to 0.5

Example: 
sumac.py -d pln -i Onagraceae -o Lythraceae

If you already downloaded the GB database:
sumac.py -i Onagraceae -o Lythraceae


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
    print(color.purple + "Connecting to ftp.ncbi.nih.gov..." + color.done)
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
    os.chdir("genbank/")
    files = os.listdir(".")
    gb = SeqIO.index_db("gb.idx", files, "genbank")
    os.chdir(os.pardir)
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
    for key in keys:
	if ingroup in gb[key].annotations['taxonomy']:
            ingroup_keys.append(key)
	elif outgroup in gb[key].annotations['taxonomy']:
	    outgroup_keys.append(key)
	sys.stdout.write('\r' + color.yellow + 'Ingroup sequences found: ' + color.red + str(len(ingroup_keys)) + color.yellow \
	    + '  Outgroup sequences found: ' + color.red + str(len(outgroup_keys)) + color.done)    
	sys.stdout.flush()    
	## FOR TESTING ONLY
	if len(ingroup_keys) == 50:  ##
            sys.stdout.write("\n")   ## FOR TESTING ONLY
            sys.stdout.flush()       ## # FOR TESTING ONLY
	    os.chdir(os.pardir)      ##
	    return ingroup_keys, outgroup_keys    # FOR TESTING ONLY
	## FOR TESTING ONLY
	## remove above
    sys.stdout.write("\n")
    sys.stdout.flush()
    os.chdir(os.pardir)
    return ingroup_keys, outgroup_keys

def make_distance_matrix(gb, seq_keys, length_threshold):
    """
    Takes as input a dictionary of SeqRecords gb and the keys to all sequences.
    length_threshold is the threshold of sequence length percent similarity to cluster taxa.
    For example if length_threshold = 0.25, and one sequence has
    length 100, the other sequence must have length 75 to 125. If the lengths are not similar
    enough the distance is set to 50 (which keeps them from being clustered).
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
		# check sequence lengths
		else:
		    length1 = len(gb[key].seq)
		    length2 = len(gb[key2].seq)
		    if (length1 > length2 * (1 + float(length_threshold))) or (length1 < length2 * (1 - float(length_threshold))):
                        dist_matrix[i][j] = 50.0
			dist_matrix[j][i] = 50.0
		    else:
		        # do the blast comparison
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
	percent = str(round(100 * i/float(len(seq_keys)), 2))
	sys.stdout.write('\r' + color.blue + 'Completed: ' + color.red + str(i) + '/' + str(len(seq_keys)) + ' (' + percent + '%)' + color.done)    
	sys.stdout.flush()    
    sys.stdout.write("\n")
    sys.stdout.flush()
    os.remove("blast.xml")
    os.remove("query.fasta")
    os.remove("subject.fasta")
    return dist_matrix

def make_clusters(seq_keys, distance_matrix, threshold=(1.0/10**10)):
    """
    Input: seq_keys a list of all sequences used in the analysis, distance_matrix based on BLAST e-values, and an optional e-value threshold for clustering.
    Output: a list of clusters (each cluster is itself a list of keys to sequences)
    This function is a wrapper around the recursive function merge_closest_clusters.
    """
    # put each sequence in its own cluster
    clusters = []
    for seq in seq_keys:
        clusters.append([seq])
    return merge_closest_clusters(clusters, distance_matrix, threshold)

def merge_closest_clusters(clusters, distance_matrix, threshold):
    """
    Input: a list of clusters, distance_matrix based on BLAST e-values, and the e-value threshold to stop clustering
    Output: a list of clusters (each cluster is itself a list of keys to sequences)
    Single-linkage hierarchical clustering algorithm.
    """
    cluster1 = 0
    cluster2 = 0
    min_distance = 99
    x = 0
    y = 0
    # find the most similar pair of clusters (or the first pair with distance = 0)
    while x < len(distance_matrix):
        y = x + 1
        while y < len(distance_matrix):
	    if x != y:
		if distance_matrix[x][y] < min_distance:
		    min_distance = distance_matrix[x][y]
		    cluster1 = x
		    cluster2 = y
		if min_distance == 0:
		    break
            y += 1
	if min_distance == 0:
     	    break
	x += 1
	
    # check to see if we are done
    if min_distance > threshold:
        return clusters
    
    # merge the two clusters
    for sequence in clusters[cluster2]:
        clusters[cluster1].append(sequence)
    del clusters[cluster2]

    # update distance matrix
    for i in range(len(distance_matrix[cluster1])):
        if distance_matrix[cluster1][i] > distance_matrix[cluster2][i]:
	    distance_matrix[cluster1][i] = distance_matrix[cluster2][i]
	    distance_matrix[i][cluster1] = distance_matrix[cluster2][i]
    del distance_matrix[cluster2]
    for i in range(len(distance_matrix)):
        del distance_matrix[i][cluster2]

    return merge_closest_clusters(clusters, distance_matrix, threshold)

def assemble_fasta_clusters(gb, clusters):
    """
    Inputs the dictionary of all GenBank sequences and a list of clustered accessions.
    Only make fasta files of clusters containing 4 sequences or more.
    Outputs a list of FASTA files, each file containing an unaligned sequence cluster,
    and an updated list of clustered accessions.
    """
    matrices = []
    if not os.path.exists("clusters"):
        os.makedirs("clusters")
    i = 0
    to_delete = []
    for cluster in clusters:
	if len(cluster) > 3:
	    sequences = []
	    for seq_key in cluster:
                sequences.append(gb[seq_key])
            file_name = "clusters/" + str(i) + ".fasta"
            file = open(file_name, "wb")
	    SeqIO.write(sequences, file, 'fasta')
	    file.close()
	    matrices.append(file_name)
	    i += 1
	else:
	    to_delete.append(cluster)
    for cluster in to_delete:
        del clusters[clusters.index(cluster)]
    return matrices, clusters

def align_matrices(matrices):
    """
    Inputs a list of FASTA files, each file containing an unaligned sequence cluster.
    Must determine whether sequences must be reverse complemented by using BLAST.
    Outputs a list of FASTA files, each file containing an aligned sequence cluster.
    """
    # do stuff
    # return aligned_matrices
    # TODO:
    # blast to check forward/reverse sequences...


def main():
    # parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--download_gb", "-d", help="Name of the GenBank division to download (e.g. PLN or MAM).")
    parser.add_argument("--ingroup", "-i", help="Ingroup clade to build supermatrix.")
    parser.add_argument("--outgroup", "-o", help="Outgroup clade to build supermatrix.")
    parser.add_argument("--max_outgroup", "-m", help="Maximum number of taxa to include in outgroup. Defaults to 10.")
    parser.add_argument("--evalue", "-e", help="BLAST E-value threshold to cluster taxa. Defaults to 0.1")
    parser.add_argument("--length", "-l", help="Threshold of sequence length percent similarity to cluster taxa. Defaults to 0.5")
    args = parser.parse_args()
   
    print("")
    print(color.blue + "SUMAC: supermatrix constructor" + color.done)
    print("")

    # download and set up sqllite db
    if args.download_gb:
        gb_division = str(args.download_gb).lower()
        download_gb_db(gb_division)
        print(color.yellow + "Setting up SQLite database..." + color.done)
        gb = setup_sqlite()
    elif not os.path.exists("genbank/gb.idx"):
        print(color.red + "GenBank database not downloaded. Re-run with the -d option. See --help for more details." + color.done)
        sys.exit(0)
    else:
        print(color.purple + "Genbank database already downloaded. Indexing sequences..." + color.done)
	gb = SeqIO.index_db("genbank/gb.idx")
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
    ingroup_keys, outgroup_keys = get_in_out_groups(gb, ingroup, outgroup)        
    all_seq_keys = ingroup_keys + outgroup_keys

    # make distance matrix
    print(color.blue + "Making distance matrix for all sequences..." + color.done)
    length_threshold = 0.5
    if args.length:
	length_threshold = args.length
    print(color.blue + "Using sequence length similarity threshold " + str(length_threshold) + color.done)
    distance_matrix = make_distance_matrix(gb, all_seq_keys, length_threshold)

    # cluster sequences
    print(color.purple + "Clustering sequences..." + color.done)
    if args.evalue:
	print(color.purple + "Using e-value threshold " + str(args.evalue) + color.done) 
	clusters = make_clusters(all_seq_keys, distance_matrix, float(args.evalue))
    else:
	print(color.purple + "Using default e-value threshold 0.1" + color.done)
        clusters = make_clusters(all_seq_keys, distance_matrix)
    print(color.purple + "Found " + str(len(clusters)) + " clusters." + color.done)

    # make FASTA files
    print(color.yellow + "Building sequence matrices for each cluster." + color.done)
    # print(color.yellow + "Checking for sequences that must be reverse complemented..." + color.done)
    unaligned_matrices, clusters = assemble_fasta_clusters(gb, clusters)

    print(color.purple + "Kept " + str(len(clusters)) + " clusters, discarded those with < 4 sequences." + color.done)

    for cluster in clusters:
	print("***cluster:")
	print(cluster)

    print(unaligned_matrices)

    # now align each cluster, then concantenate, then reduce the number of outgroup taxa

	
if __name__ == "__main__":
    main()
