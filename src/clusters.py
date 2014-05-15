#! /usr/bin/python

import os
import sys
import multiprocessing
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from util import Color


class Clusters(object):
    """
    Class responsible for managing clusters of homologous sequences.
    """

    clusters = []
    seq_keys = []

    def __init__(self, seq_keys):
        self.seq_keys = seq_keys

    def write_fasta(self):
        return True


class DistanceMatrixClusters(Clusters):
    """
    Builds clusters from a distance matrix.
    Inherits from Clusters.
    """

    distance_matrix = []

    def __init__(self, seq_keys, distance_matrix, threshold=(1.0/10**10)):
        Clusters.__init__(self, seq_keys)
        self.distance_matrix = distance_matrix
        return True

class GuidedClusters(Clusters):
    """
    Builds clusters from guide sequences.
    Inherits from Clusters.
    """

    def __init__(self, seq_keys):
        Clusters.__init__(self, seq_keys)
        return True



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
            row = distance_matrix[cluster1]
            row[i] = distance_matrix[cluster2][i]
            distance_matrix[cluster1] = row
            row = distance_matrix[i]
            row[cluster1] = distance_matrix[cluster2][i]
            distance_matrix[i] = row
    del distance_matrix[cluster2]
    for i in range(len(distance_matrix)):
        row = distance_matrix[i]
        del row[cluster2]
        distance_matrix[i] = row

    return merge_closest_clusters(clusters, distance_matrix, threshold)



def make_guided_clusters(guide_seq, all_seq_keys, length_threshold, evalue_threshold):
    """
    Input: name of FASTA file containing guide sequences, dictionary of all GenBank sequences,
    a list of ingroup/outgroup sequences, the e-value threshold to cluster, and the
    threshold of sequence length percent similarity to cluster taxa.
    Returns a list of clusters (each cluster is itself a list of keys to sequences).
    """
    lock = multiprocessing.Lock()
    manager = multiprocessing.Manager()
    already_compared = manager.list()
    clusters = manager.list()

    color = Color()
    # check for fasta file of guide sequences
    if not os.path.isfile(guide_seq):
        print(color.red + "FASTA file of guide sequences not found. Please re-try." + color.done)
        sys.exit(0)
    else:
        # initialize an empty list for each cluster
        guide_sequences = SeqIO.parse(open(guide_seq, "rU"), "fasta")
        for guide in guide_sequences:
            clusters.append([])

    num_cores = multiprocessing.cpu_count()
    print(color.blue + "Spawning " + color.red + str(num_cores) + color.blue + " processes to make clusters." + color.done)
    processes = []
    
    for i in range(num_cores):
        p = multiprocessing.Process(target=make_guided_clusters_worker, args=(guide_seq, all_seq_keys, \
            length_threshold, evalue_threshold, clusters, already_compared, lock, i))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
    
    sys.stdout.write("\n")
    sys.stdout.flush()
    return clusters



def make_guided_clusters_worker(guide_seq, all_seq_keys, length_threshold, evalue_threshold, clusters, already_compared, lock, process_num):
    """
    Worker process for make_guided_clusters(). Each process will compare all the ingroup/outgroup sequences
    to a guide sequence, adding that guide sequence to the already_compared list.
    """
    # each process must load its own sqlite gb
    gb_dir = os.path.abspath("genbank/")
    gb = SeqIO.index_db(gb_dir + "/gb.idx")
    process_num = str(process_num)

    color = Color()

    # open guide fasta file
    if os.path.isfile(guide_seq):
        guide_sequences = list(SeqIO.parse(open(guide_seq, "rU"), "fasta"))
    else:        
        print(color.red + "FASTA file of guide sequences not found. Please re-try." + color.done)
        sys.exit(0)

    # remember how many guide sequences there are
    num_guides = len(list(guide_sequences))

    for guide in guide_sequences:
        # check whether another process is already comparing this guide sequence
        compare_guide = False
        with lock:
            if guide.id not in already_compared:
                already_compared.append(guide.id)
                compare_guide = True
        if compare_guide:
            # loop through each ingroup/outgroup sequence and blast to guide seq
            output_handle = open('subject' + process_num + '.fasta', 'w')
            SeqIO.write(guide, output_handle, 'fasta')
            output_handle.close()
            for key in all_seq_keys:
                record = gb[key]
                length1 = len(guide.seq)
                length2 = len(record.seq)
                # check if length similarity threshold met
                if (length1 < length2 * (1 + float(length_threshold))) and (length1 > length2 * (1 - float(length_threshold))):
                    # do the blast comparison
                    output_handle = open('query' + process_num + '.fasta', 'w')
                    SeqIO.write(record, output_handle, 'fasta')
                    output_handle.close()

                    blastn_cmd = NcbiblastnCommandline(query='query' + process_num + '.fasta', subject='subject' + process_num + \
                        '.fasta', out='blast' + process_num + '.xml', outfmt=5)
                    stdout, stderr = blastn_cmd()
                    blastn_xml = open('blast' + process_num + '.xml', 'r')
                    blast_records = NCBIXML.parse(blastn_xml)

                    for blast_record in blast_records:
                        if blast_record.alignments:
                            if blast_record.alignments[0].hsps:
                                # blast hit found, add sequence to cluster
                                with lock:
                                    temp_cluster = clusters[guide_sequences.index(guide)]
                                    temp_cluster.append(key)
                                    clusters[guide_sequences.index(guide)] = temp_cluster
                    blastn_xml.close()
        # update status
        percent = str(round(100 * len(already_compared)/float(num_guides), 2))
        sys.stdout.write('\r' + color.blue + 'Completed: ' + color.red + str(len(already_compared)) + '/' + str(num_guides) + ' (' + percent + '%)' + color.done)    
        sys.stdout.flush()    
    # done looping through all guides, now clean up
    if os.path.isfile("blast" + process_num + ".xml"):
        os.remove("blast" + process_num + ".xml")
    if os.path.isfile("query" + process_num + ".fasta"):
        os.remove("query" + process_num + ".fasta")
    if os.path.isfile("subject" + process_num + ".fasta"):
        os.remove("subject" + process_num + ".fasta")


