"""
SUMAC: supermatrix constructor

Copyright 2014 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


import os
from os import listdir
from os.path import isfile, join
import subprocess
from subprocess import CalledProcessError
import sys
import multiprocessing
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from util import Color


class ClusterBuilder(object):
    """
    Class responsible for assembling clusters of homologous sequences.
    """

    clusters = []
    seq_keys = []
    cluster_files = []


    def __init__(self, seq_keys):
        self.seq_keys = seq_keys


    def write_fasta(self):
        return True


    def assemble_fasta(self, gb, min_clusters=4):
        """
        Inputs the dictionary of all GenBank sequence.
        Only make fasta files of clusters containing min_clusters taxa or more,
        and delete those clusters with less than min_clusters.
        Generates a list of FASTA files, each file containing an unaligned sequence cluster.
        """
        cluster_files = []
        if not os.path.exists("clusters"):
            os.makedirs("clusters")
        i = 1
        to_delete = []
        for cluster in self.clusters:
            # get all OTUs in cluster
            otus = []
            for seq_key in cluster:
                descriptors = gb[seq_key].description.split(" ")
                otu = descriptors[0] + " " + descriptors[1]
                if otu not in otus:
                    otus.append(otu)
            # make fasta file if >= min_clusters OTUs in cluster
            otus_in_cluster = []
            if len(otus) >= min_clusters:
                sequences = []
                for seq_key in cluster:
                    descriptors = gb[seq_key].description.split(" ")
                    otu = descriptors[0] + " " + descriptors[1]
                    # do not allow duplicate OTUs in cluster
                    if otu not in otus_in_cluster:
                        sequences.append(gb[seq_key])
                        otus_in_cluster.append(otu)
                file_name = "clusters/" + str(i) + ".fasta"
                file = open(file_name, "wb")
                SeqIO.write(sequences, file, 'fasta')
                file.close()
                cluster_files.append(file_name)
                i += 1
            else:
                to_delete.append(cluster)
        for cluster in to_delete:
            del self.clusters[self.clusters.index(cluster)]
        self.cluster_files = cluster_files


    
    def assemble_fasta_uclust(self, min_clusters=4):
        """
        Only make fasta files of clusters containing min_clusters taxa or more,
        and delete those clusters with less than min_clusters.
        Generates a list of FASTA files, each file containing an unaligned sequence cluster.
        """
        cluster_files = []
        if not os.path.exists("clusters"):
            os.makedirs("clusters")
        i = 1
        to_delete = []
        for cluster in self.clusters:
            # get all OTUs in cluster
            otus = []
            for record in SeqIO.parse(open("uclusters/" + cluster, "rU"), "fasta"):
                descriptors = record.id.split("_")
                otu = descriptors[1] + " " + descriptors[2]
                if otu not in otus:
                    otus.append(otu)
            # make fasta file if >= min_clusters OTUs in cluster
            otus_in_cluster = []
            if len(otus) >= min_clusters:
                sequences = []
                for record in SeqIO.parse(open("uclusters/" + cluster, "rU"), "fasta"):
                    descriptors = record.id.split("_")
                    otu = descriptors[1] + " " + descriptors[2]
                    # do not allow duplicate OTUs in cluster
                    if otu not in otus_in_cluster:
                        sequences.append(record)
                        otus_in_cluster.append(otu)
                file_name = "clusters/" + str(i) + ".fasta"
                temp_file_name = "clusters/temp" + str(i) + ".fasta"
                f = open(temp_file_name, "wb")
                SeqIO.write(sequences, f, 'fasta')
                f.close()
                with open(temp_file_name, "r") as f, open(file_name, "w") as fout:
                    for l in f:
                        if ">" in l:
                            l = l.replace("_", " ")    
                        fout.write(l)
                subprocess.check_call(["rm", temp_file_name])
                cluster_files.append(file_name)
                i += 1
            else:
                to_delete.append(cluster)
        for cluster in to_delete:
            del self.clusters[self.clusters.index(cluster)]
        self.cluster_files = cluster_files



class SLINKClusterBuilder(ClusterBuilder):
    """
    Clusters sequences using the SLINK single-linkage clustering algorithm, which has a O(n^2) time complexity.
    Inherits from ClusterBuilder.
    """


    distance_matrix = []
    threshold = (1.0/10**10)


    def __init__(self, seq_keys, distance_matrix, threshold=(1.0/10**10)):
        """
        Input: seq_keys a list of all sequences used in the analysis, distance_matrix based on BLAST e-values, and an optional e-value threshold for clustering.
        Output: a list of clusters (each cluster is itself a list of keys to sequences)
        This function is a wrapper around the recursive function merge_closest_clusters.
        """
        ClusterBuilder.__init__(self, seq_keys)
        self.distance_matrix = distance_matrix
        self.seq_keys = seq_keys
        self.threshold = threshold
        color = Color()

        # pointer representation of cluster hierarchy:
        # Pi[i] is the first cluster that cluster i joins
        # Lambda[i] is the distance between cluster i and cluster Pi[i]
        n = len(seq_keys)
        Pi = [0] * n
        Lambda = [0] * n
        M = [0] * n

        Pi[0] = 0
        Lambda[0] = float("inf")
        for i in range(1, n):
            # update status
            percent = str(round(100 * i/float(n), 2))
            sys.stdout.write('\r' + color.blue + 'Completed: ' + color.red + str(i) + '/' + str(n) + ' (' + percent + '%)' + color.done)    
            sys.stdout.flush()    
            
            Pi[i] = i
            Lambda[i] = float("inf")

            for j in range(i):
                M[j] = distance_matrix[i][j]

            for j in range(i):
                if Lambda[j] >= M[j]:
                    M[Pi[j]] = min(M[Pi[j]], Lambda[j])
                    Lambda[j] = M[j]
                    Pi[j] = i
                else:
                    M[Pi[j]] = min(M[Pi[j]], M[j])

            for j in range(i):
                if Lambda[j] >= Lambda[Pi[j]]:
                    Pi[j] = i

        # convert from pointer representation to list of sequence clusters
        sys.stdout.write("\n")
        sys.stdout.flush()
        print(color.blue + "Finalizing clusters..." + color.done)
        temp_clusters = [[] for _ in range(n)]

        for i in range(n):
            if Lambda[i] <= self.threshold:
                temp_clusters[Pi[i]].append(i)
                if len(temp_clusters[i]) > 0:
                    for j in temp_clusters[i]:
                        temp_clusters[Pi[i]].append(j)
                    temp_clusters[i] = []
            else:
                temp_clusters[i].append(i)

        for temp_cluster in temp_clusters:
            if len(temp_cluster) > 0:
                temp_cluster_seq = []
                for i in temp_cluster:
                    temp_cluster_seq.append(seq_keys[i])
                self.clusters.append(temp_cluster_seq)



class HACClusterBuilder(ClusterBuilder):
    """
    Clusters sequences using the naive hierarchical agglomerative clustering (HAC) algorithm, 
    which has a O(n^3) time complexity.
    Inherits from ClusterBuilder.
    """


    distance_matrix = []
    threshold = (1.0/10**10)


    def __init__(self, seq_keys, distance_matrix, threshold=(1.0/10**10)):
        """
        Input: seq_keys a list of all sequences used in the analysis, distance_matrix based on BLAST e-values, and an optional e-value threshold for clustering.
        Output: a list of clusters (each cluster is itself a list of keys to sequences)
        This function is a wrapper around the recursive function merge_closest_clusters.
        """
        ClusterBuilder.__init__(self, seq_keys)
        self.distance_matrix = distance_matrix
        self.seq_keys = seq_keys
        self.threshold = threshold

        # put each sequence in its own cluster
        for seq in seq_keys:
            self.clusters.append([seq])
        
        sys.setrecursionlimit(len(self.clusters) + 10)
        self.merge_closest_clusters(self.clusters, distance_matrix)


    def merge_closest_clusters(self, clusters, distance_matrix):
        """
        Input: a list of clusters, distance_matrix based on BLAST e-values
        Output: a list of clusters (each cluster is itself a list of keys to sequences)
        Recursive function for the naive single-linkage hierarchical clustering algorithm.
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
        if min_distance > self.threshold:
            self.clusters = clusters
            return
        
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

        self.merge_closest_clusters(clusters, distance_matrix)



class UCLUSTClusterBuilder(ClusterBuilder):
    """
    Clusters sequences using the UCLUST algorithm. 
    Inherits from ClusterBuilder.
    """


    threshold = 0.75
    error = False

    def __init__(self, gb, seq_keys, gb_dir, num_cores, minlength, maxlength, length_thres=0.5, threshold=0.75, evalue=(1.0/10**10)):
        """
        Input: gb dictionary of SeqRecords, keys to all sequences, and an optional threshold for clustering.
        Output: a list of cluster files from UCLUST
        """
        ClusterBuilder.__init__(self, seq_keys)
        self.seq_keys = seq_keys
        self.threshold = threshold
        color = Color()

        if not os.path.exists("uclusters"):
            os.makedirs("uclusters")
        
        # write sequences to fasta
        sequences = []
        for seq_key in seq_keys:
            record = gb[seq_key]
            record.description = record.annotations["organism"] + " " + record.description
            if "sp." not in record.annotations["organism"]:
                sequences.append(record)
        file_name = "_sumac"
        f = open(file_name, "wb")
        SeqIO.write(sequences, f, 'fasta')
        f.close()
        with open("_sumac", "r") as f, open("_sumac_filtered", "w") as fout:
            for l in f:
                if ">" in l:
                    l = l.replace(" ", "_")    
                fout.write(l)

        # call UCLUST
        sort_sequences = ["usearch", "-sortbylength", "_sumac_filtered", "-fastaout", "_sumac_sorted",
                          "-minseqlength", str(minlength), "-maxseqlength", str(maxlength)]
        uclust = ["usearch", "-cluster_fast", "_sumac_sorted", "-id", str(threshold),
                  "-minsl", str(length_thres), "-strand", "both", "-threads", str(num_cores), 
                  "-clusters", "uclusters/", "-fulldp", "-evalue", str(evalue)]
        try:
            subprocess.check_call(sort_sequences)
            subprocess.check_call(uclust)
        except CalledProcessError as e:
            print(color.red + "UCLUST error: " + str(e) + color.done)
            print(color.red + "Trying SLINK instead..." + color.done)
            self.error = True
            return
        except OSError as e:
            print(color.red + "UCLUST is not installed correctly." + color.done)
            print(color.red + "OS error: " + str(e) + color.done)
            print(color.red + "Trying SLINK instead..." + color.done)
            self.error = True
            return
        finally:
            if os.path.exists("_sumac"):
                subprocess.check_call(["rm", "_sumac"])
            if os.path.exists("_sumac_filtered"):
                subprocess.check_call(["rm", "_sumac_filtered"])
            if os.path.exists("_sumac_sorted"):
                subprocess.check_call(["rm", "_sumac_sorted"])
        cluster_files = [ f for f in listdir("uclusters/") if isfile(join("uclusters/", f)) ]
        for f in cluster_files:
            self.clusters.append(f)




class GuidedClusterBuilder(ClusterBuilder):
    """
    Builds clusters from guide sequences.
    Inherits from ClusterBuilder.
    """


    def __init__(self, guide_seq, all_seq_keys, length_threshold, evalue_threshold, gb_dir, num_cores):
        """
        Input: name of FASTA file containing guide sequences, dictionary of all GenBank sequences,
        a list of ingroup/outgroup sequences, the e-value threshold to cluster, and the
        threshold of sequence length percent similarity to cluster taxa,
        and the GenBank directory.
        Generates a list of clusters (each cluster is itself a list of keys to sequences).
        """
        ClusterBuilder.__init__(self, all_seq_keys)
        
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

        # make blast database
        gb = SeqIO.index_db(gb_dir + "/gb.idx")
        output_handle = open('blast_db.fasta', 'w')
        records = []
        for key in all_seq_keys:
            record = gb[key]
            records.append(record)
        SeqIO.write(records, output_handle, 'fasta')
        output_handle.close()

        # spawn processes
        print(color.blue + "Spawning " + color.red + str(num_cores) + color.blue + " processes to make clusters." + color.done)
        processes = []
        
        for i in range(num_cores):
            p = multiprocessing.Process(target=self.make_guided_clusters_worker, args=(guide_seq, all_seq_keys, \
                length_threshold, evalue_threshold, clusters, already_compared, lock, i, gb_dir))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()
        
        sys.stdout.write("\n")
        sys.stdout.flush()
        self.clusters = clusters
        if os.path.isfile("blast_db.fasta"):
            os.remove("blast_db.fasta")



    def make_guided_clusters_worker(self, guide_seq, all_seq_keys, length_threshold, evalue_threshold, clusters, already_compared, lock, process_num, gb_dir):
        """
        Worker process for make_guided_clusters(). Each process will compare all the ingroup/outgroup sequences
        to a guide sequence, adding that guide sequence to the already_compared list.
        """

        color = Color()
        process_num = str(process_num)

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
                
                # set query to be guide sequence 
                output_handle = open('query' + process_num + '.fasta', 'w')
                SeqIO.write(guide, output_handle, 'fasta')
                output_handle.close()

                # blast query against blast_db
                blastn_cmd = NcbiblastnCommandline(query='query' + process_num + '.fasta', subject='blast_db.fasta', \
                    out='blast' + process_num + '.xml', outfmt=5)
                stdout, stderr = blastn_cmd()

                # parse blast output
                blastn_xml = open('blast' + process_num + '.xml', 'r')
                blast_records = NCBIXML.parse(blastn_xml)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        # loop through each high-scoring segment pair (HSP)
                        for hsp in alignment.hsps:
                            length1 = len(guide.seq)
                            length2 = alignment.length
                            # check if length similarity threshold met
                            # and check to see if evalue_threshold is met
                            if (length2 < length1 * (1 + float(length_threshold))) and (length2 > length1 * (1 - float(length_threshold))) \
                                and (hsp.expect < evalue_threshold):
                                # blast hit found, add sequence to cluster
                                # get accession number from sequence title:
                                # Subject_145 AJ620515.1 Gaura parviflora ITS1, 5.8S rRNA gene, and ITS2
                                accession = alignment.title.split(" ")[1]
                                with lock:
                                    temp_cluster = clusters[guide_sequences.index(guide)]
                                    temp_cluster.append(accession)
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


