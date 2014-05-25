#! /usr/bin/python

import os
import sys
import multiprocessing
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from util import Color

class DistanceMatrixBuilder:
    """
    Builds distance matrix
    """
    
    distance_matrix = []

    def __init__(self, gb, seq_keys, length_threshold, gb_dir):
        """
        Takes as input a dictionary of SeqRecords gb and the keys to all sequences.
        length_threshold is the threshold of sequence length percent similarity to cluster taxa.
        For example if length_threshold = 0.25, and one sequence has
        length 100, the other sequence must have length 75 to 125. If the lengths are not similar
        enough the distance is set to 50 (which keeps them from being clustered).
        Generates a 2 dimensional list of distances. Distances are blastn e-values.
        """
        lock = multiprocessing.Lock()
        manager = multiprocessing.Manager()
        already_compared = manager.list()
        dist_matrix = manager.list()
        row = []
        for i in range(len(seq_keys)):
            row.append(99)
        for i in range(len(seq_keys)):
            dist_matrix.append(row)

        num_cores = multiprocessing.cpu_count()
        print(Color.blue + "Spawning " + Color.red + str(num_cores) + Color.blue + " processes to make distance matrix." + Color.done)
        processes = []

        for i in range(num_cores):
            p = multiprocessing.Process(target=self.distance_matrix_worker, args=(seq_keys, length_threshold, dist_matrix, already_compared, lock, i, gb_dir))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

        sys.stdout.write("\n")
        sys.stdout.flush()
        self.distance_matrix = dist_matrix


    def distance_matrix_worker(self, seq_keys, length_threshold, dist_matrix, already_compared, lock, process_num, gb_dir):
        """
        Worker process for make_distance_matrix(). Takes a list "already_compared" of sequences that have
        already had all pairwise comparisons. Each worker process will work making pairwise comparisons
        for a different sequence, adding them to the "already_compared" list as they are completed.
        """
        # each process must load its own sqlite gb
        gb = SeqIO.index_db(gb_dir + "/gb.idx")
        process_num = str(process_num)
        i = 0
        for key in seq_keys:
            # check whether another process is already comparing this row
            compare_row = False
            with lock:
                if key not in already_compared:
                    already_compared.append(key)
                    compare_row = True
            if compare_row:
                # get the sequence record to compare
                record1 = gb[key]
                output_handle = open('subject' + process_num + '.fasta', 'w')
                SeqIO.write(record1, output_handle, 'fasta')
                output_handle.close()
                j = 0
                for key2 in seq_keys:
                    # only calculate e-values for pairs that have not yet been compared
                    if dist_matrix[i][j] == 99:
                        if key == key2:
                            row = dist_matrix[i]
                            row[j] = 0.0
                            dist_matrix[i] = row
                        # check sequence lengths
                        else:
                            # print("proc # = "+process_num+" i = "+str(i)+ " j = "+str(j))
                            record2 = gb[key2]
                            length1 = len(record1.seq)
                            length2 = len(record2.seq)
                            # set distance to 50.0 if length similarity threshold not met
                            if (length1 > length2 * (1 + float(length_threshold))) or (length1 < length2 * (1 - float(length_threshold))):
                                row = dist_matrix[i]
                                row[j] = 50.0
                                dist_matrix[i] = row
                                row = dist_matrix[j]
                                row[i] = 50.0
                                dist_matrix[j] = row
                            else:
                                # do the blast comparison
                                output_handle = open('query' + process_num + '.fasta', 'w')
                                SeqIO.write(record2, output_handle, 'fasta')
                                output_handle.close()

                                blastn_cmd = NcbiblastnCommandline(query='query' + process_num + '.fasta', subject='subject' + process_num + \
                                    '.fasta', out='blast' + process_num + '.xml', outfmt=5)
                                stdout, stderr = blastn_cmd()
                                blastn_xml = open('blast' + process_num + '.xml', 'r')
                                blast_records = NCBIXML.parse(blastn_xml)

                                for blast_record in blast_records:
                                    if blast_record.alignments:
                                        if blast_record.alignments[0].hsps:
                                            # blast hit found, set distance to e-value
                                            row = dist_matrix[i]
                                            row[j] = blast_record.alignments[0].hsps[0].expect
                                            dist_matrix[i] = row
                                            row = dist_matrix[j]
                                            row[i] = blast_record.alignments[0].hsps[0].expect
                                            dist_matrix[j] = row
                                    else:
                                        # no blast hit found, set distance to default 10.0
                                        row = dist_matrix[i]
                                        row[j] = 10.0
                                        dist_matrix[i] = row
                                        row = dist_matrix[j]
                                        row[i] = 10.0
                                        dist_matrix[j] = row
                                blastn_xml.close()
                    j += 1
            i += 1
            # update status
            percent = str(round(100 * len(already_compared)/float(len(seq_keys)), 2))
            sys.stdout.write('\r' + Color.blue + 'Completed: ' + Color.red + str(len(already_compared)) + '/' + str(len(seq_keys)) + ' (' + percent + '%)' + Color.done)
            sys.stdout.flush()
        # done looping through all keys, now clean up
        os.remove("blast" + process_num + ".xml")
        os.remove("query" + process_num + ".fasta")
        os.remove("subject" + process_num + ".fasta")


