#! /usr/bin/python


import os
import sys
import argparse
import multiprocessing
import copy_reg
import types
from Bio import Entrez
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from util import Color



class Alignments(object):
    """
    Class responsible for managing sequence alignments
    """


    files = []


    def __init__(self, cluster_files):
        """
        Inputs a list of FASTA files, each file containing an unaligned sequence cluster.
        Creates new processes to align each sequence cluster.
        Generates a list of aligned FASTA files.
        """
        alignment_files = []
        if not os.path.exists("alignments"):
            os.makedirs("alignments")
        color = Color()
        print(color.blue + "Spawning " + color.red + str(multiprocessing.cpu_count()) + color.blue + " processes to align clusters." + color.done)
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        alignment_files = pool.map(self.align_cluster, cluster_files)
        pool.close()
        pool.join()
        self.files = alignment_files



    def align_cluster(self, cluster_file):
        """
        Worker fuction for align_clusters
        Inputs a FASTA file containing an unaligned sequence cluster.
        Uses MUSCLE to align the cluster.
        """
        # TODO:   
        # use MAFFT instead of MUSCLE to deal with forward/reverse strand issue
        alignment_file = "alignments" + cluster_file[cluster_file.index("/"):]
        muscle_cline = MuscleCommandline(input=cluster_file, out=alignment_file)
        color = Color()
        print(color.red + str(muscle_cline) + color.done)
        sys.stdout.flush()
        stdout, stderr = muscle_cline()
        return alignment_file



    def print_region_data(self):
        """
        Prints the name of each DNA region, the number of taxa, the aligned length,
        missing data (%), and taxon coverage density
        """
        # TODO: calculate the variable characters (and %), PI characters (and %)
        # first get list of all taxa
        taxa = []
        for alignment in self.files:
            records = SeqIO.parse(alignment, "fasta")
            for record in records:
                descriptors = record.description.split(" ")
                taxon = descriptors[1] + " " + descriptors[2]
                if taxon not in taxa:
                    taxa.append(taxon)

        # print data for each region
        i = 0
        color = Color()
        for alignment in self.files:
            records = list(SeqIO.parse(alignment, "fasta"))
            descriptors = records[0].description.split(" ")
            print(color.blue + "Aligned cluster #: " + color.red + str(i) + color.done)
            print(color.yellow + "DNA region: " + color.red + " ".join(descriptors[3:]) + color.done)
            print(color.yellow + "OTUs: " + color.red + str(len(records)) + color.done)
            print(color.yellow + "Aligned length: " + color.red + str(len(records[0].seq)) + color.done)
            print(color.yellow + "Missing data (%): " + color.red + str(round(100 - (100 * len(records)/float(len(taxa))), 1)) + color.done)
            print(color.yellow + "Taxon coverage density: "  + color.red + str(round(len(records)/float(len(taxa)), 2)) + color.done)
            i += 1



    def concatenate(self):
        """
        Returns a concatenated FASTA file
        """
        # first build up a dictionary of all OTUs
        otus = {}
        for alignment in self.files:
            records = SeqIO.parse(alignment, "fasta")
            for record in records:
                # sample record.description:
                # AF495760.1 Lythrum salicaria chloroplast ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit-like mRNA, partial sequence
                descriptors = record.description.split(" ")
                otu = descriptors[1] + " " + descriptors[2]
                if otu not in otus:
                    otus[otu] = ""

        # now concatenate the sequences
        total_length = 0
        for alignment in self.files: 
            records = SeqIO.parse(alignment, "fasta")
            # make sure to only add 1 sequence per cluster for each otu
            already_added = []
            for record in records:
                descriptors = record.description.split(" ")
                otu = descriptors[1] + " " + descriptors[2]
                if otu not in already_added:
                    otus[otu] = otus[otu] + record.seq
                    already_added.append(otu)
                loci_length = len(record.seq)
            total_length += loci_length
            # add gaps for any OTU that didn't have a sequence
            for otu in otus:
                if len(otus[otu]) < total_length:
                    otus[otu] = otus[otu] + self.make_gaps(loci_length)

        # write to FASTA file
        f = open("alignments/combined.fasta", "w")
        for otu in otus:
            # >otu
            # otus[otu]
            f.write("> " + otu + "\n")
            sequence = str(otus[otu])
            i = 0
            while i < len(sequence):
                f.write(sequence[i:i+80] + "\n")
                i += 80
        f.close()
        return "alignments/combined.fasta"



    def make_gaps(self, length):
        """
        Inputs an integer.
        Returns a string of '-' of length
        """
        gap = ""
        for i in range(length):
            gap = "-" + gap
        return gap



    def print_supermatrix_attributes(self):
        """
        Prints out details on the final aligned super matrix.
        """
        # TODO: make the output of this more useful
        color = Color()
        print(color.blue + "Supermatrix attributes:")
        records = SeqIO.parse("alignments/combined.fasta", "fasta")
        num_records = 0
        total_gap = 0
        for record in records:
            otu = record.description
            gap = 0
            for letter in record.seq:
                if letter == '-':
                    gap += 1
                    total_gap += 1
            print(color.yellow + "OTU: " + color.red + otu + color.yellow + " % gaps = " + color.red + str(round(gap/float(len(record.seq)), 2)))
            num_records += 1
            matrix_length = len(record.seq)
        print(color.blue + "Total number of OTUs = " + color.red + str(num_records))
        print(color.blue + "Total length of matrix = " + color.red + str(matrix_length))
        print(color.blue + "Total % gaps = " + color.red + str(round(total_gap/float(matrix_length * num_records), 2)) + color.done)



# pickle method recipe by Steven Bethard
# see http://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-pythons-multiprocessing-pool-ma/
def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)


copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
