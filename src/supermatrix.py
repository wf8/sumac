"""
SUMAC: supermatrix constructor

Copyright 2014 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


import os
import sys
import csv
from Bio import Entrez
from Bio import SeqIO
from util import Color



class Supermatrix(object):
    """
    Class responsible for managing the final supermatrix
    """


    file = ""           # FASTA file of final supermatrix
    otus = {}           # dictionary of Otu objects
    pd = None           # fraction of triples, a measure of the partial decisiveness of a supermatrix
    missing_data = None # % of sequence data missing


    def __init__(self, alignments=None):
        """
        Optionally accept an Alignment object.
        """
        self.file = ""
        self.otus = {}
        self.pd = None
        self.missing_data = None
        if alignments is not None:
            self.concatenate(alignments)



    def concatenate(self, alignments):
        """
        Builds a supermatrix from a set of alignments.
        """
        otus = {}
        for alignment in alignments.files:
            records = SeqIO.parse(alignment, "fasta")
            for record in records:
                # sample record.description:
                # AF495760.1 Lythrum salicaria chloroplast ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit-like mRNA, partial sequence
                descriptors = record.description.split(" ")
                otu = descriptors[1] + " " + descriptors[2]
                if otu not in otus.keys():
                    otus[otu] = Otu(otu)

        # now concatenate the sequences
        total_length = 0
        for alignment in alignments.files: 
            records = SeqIO.parse(alignment, "fasta")
            # make sure to only add 1 sequence per cluster for each otu
            already_added = []
            for record in records:
                descriptors = record.description.split(" ")
                otu = descriptors[1] + " " + descriptors[2]
                if otu not in already_added:
                    otus[otu].update(record.seq, descriptors[0], self.get_ungapped_length(record.seq))
                    already_added.append(otu)
                loci_length = len(record.seq)
            total_length += loci_length
            # add gaps for any OTU that didn't have a sequence
            for otu in otus:
                if len(otus[otu].sequence) < total_length:
                    otus[otu].update(self.make_gaps(loci_length), "-", 0)

        # write to FASTA file
        f = open("alignments/combined.fasta", "w")
        for otu in otus:
            # >otu
            # otus[otu]
            f.write("> " + otu + "\n")
            sequence = str(otus[otu].sequence)
            i = 0
            while i < len(sequence):
                f.write(sequence[i:i+80] + "\n")
                i += 80
        f.close()
        self.file = "alignments/combined.fasta"
        self.otus = otus



    def make_gaps(self, length):
        """
        Inputs an integer.
        Returns a string of '-' of length
        """
        gap = ""
        for i in range(length):
            gap = "-" + gap
        return gap



    def get_ungapped_length(self, sequence):
        """
        Inputs a sequence, and returns the length of the sequence minus any gaps ('-')
        """
        length = 0
        for i in sequence:
            if i != "-":
                length += 1
        return length



    def print_data(self):
        """
        Prints out details on the final aligned supermatrix.
        """
        # TODO: make the output of this more useful
        color = Color()
        print(color.blue + "Supermatrix attributes:")
        records = SeqIO.parse(self.file, "fasta")
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
        print(color.blue + "Total % missing data = " + color.red + str(self.get_missing_data()) + color.done)
        print(color.blue + "Partial decisiveness (fraction of triples) = " + color.red + str(self.get_PD()) + color.done) 
        #for otu in self.otus: 
        #    self.otus[otu].print_data()



    def get_missing_data(self):
        """
        Method to calculate the percent of sequence data missing from the supermatrix.
        """
        if self.missing_data != None:
            return self.missing_data
        else:
            total_seq = 0
            missing_seq = 0
            for otu in self.otus:
                for accession in self.otus[otu].accessions:
                    total_seq += 1
                    if accession == "-":
                        missing_seq += 1
            self.missing_data = round(missing_seq/float(total_seq), 2)
            return self.missing_data



    def make_genbank_csv(self):
        """
        Method to generate a CSV file with all GenBank accessions used in supermatrix.
        """
        with open('genbank_accessions.csv', 'wb') as csv_output:
            csvwriter = csv.writer(csv_output)
            header = ["Taxa"]
            i = 1
            for accession in self.otus[self.otus.keys()[0]].accessions:
                header.append("Gene Region " + str(i))
                i += 1
            csvwriter.writerow(header)
            for otu in self.otus:
                row = [self.otus[otu].name]
                for accession in self.otus[otu].accessions:
                    row.append(accession)
                csvwriter.writerow(row)



    def normalize(self):
        """
        function to normalize the sequence length data - this is necessary for plotting the data
        """
        # initialize the list
        longest_sequences = []
        for length in self.otus[self.otus.keys()[0]].sequence_lengths:
            longest_sequences.append(length)

        # now get maximum values
        for otu in self.otus:
            i = 0
            for length in self.otus[otu].sequence_lengths:
                if length > longest_sequences[i]:
                    longest_sequences[i] = length
                i += 1

        # save the normalized data in each OTU object
        for otu in self.otus:
            i = 0
            self.otus[otu].normalized_sequence_lengths = []
            for length in self.otus[otu].sequence_lengths:
                self.otus[otu].normalized_sequence_lengths.append((100 * length) / longest_sequences[i])
                i += 1



    def make_figure(self):
        
        # TODO:
        # test for numpy and matplotlib, show message if not present
        import numpy as np
        import matplotlib.pyplot as plt

        self.normalize()

        data = []
        otu_names = []
        genes = []
        i = 1
        for otu in self.otus:
            data.append(self.otus[otu].normalized_sequence_lengths)
            otu_names.append(otu)
            genes.append(str(i))
            i += 1

        supermatrix = np.array(data)

        # set up figure
        fig, ax = plt.subplots()
        ax.set_position([0.4, 0.1, .3, .7])
        # fig.set_size_inches(8.5, 11)

        # add data
        heatmap = ax.pcolor(supermatrix, cmap=plt.cm.Blues, alpha=0.8)

        # put the labels in the middle of each cell
        ax.set_yticks(np.arange(supermatrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(supermatrix.shape[1]) + 0.5, minor=False)

        # move x axis label stuff to top
        ax.invert_yaxis()
        ax.xaxis.tick_top()

        # add the labels
        ax.set_xticklabels(genes, minor=False, family="Arial", size=10)
        ax.set_yticklabels(otu_names, minor=False, family="Arial", size=8)

        # rotate the gene names
        #plt.xticks(rotation=90)

        # remove junk from axes
        ax.grid(False)
        ax.set_frame_on(False)

        # turn off all the ticks
        for t in ax.xaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        for t in ax.yaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False

        # add color legend
        axcolor = fig.add_axes([0.425, .90, 0.25, 0.015])
        cbar = fig.colorbar(heatmap, cax=axcolor, ticks=[0,50, 100], orientation='horizontal')
        cbar.ax.set_xticklabels(['0%', '50%', '100%'], family="Arial", size=10)

        plt.savefig("plot.pdf")



    def get_PD(self):
        """
        method to get pd of supermatrix, calculating it if need be.
        """
        if self.pd == None:
            self.pd = self.calculate_PD()
        return self.pd



    def calculate_PD(self):
        """
        Method to calculate the fraction of triples, a measure of partial decisiveness (PD).
        See: Sanderson, M.J., McMahon, M.M. & Steel, M., 2010. BMC evolutionary biology, 10. 
        """
        triplets = []
        i = 0
        for otu1 in self.otus:
            if i < (len(self.otus) - 2):
                j = 0
                for otu2 in self.otus:
                    if i < j:
                        k = 0
                        for otu3 in self.otus:
                            if j < k:
                                triplet = [otu1, otu2, otu3]
                                triplets.append(triplet)
                            k += 1
                    j += 1
            i += 1

        decisive_triples = 0
        total_triples = len(triplets)
        for triplet in triplets:
            i = 0
            for seq_length in self.otus[triplet[0]].sequence_lengths:
                if seq_length == 0:
                    s1 = 0
                else:
                    s1 = 1
                if self.otus[triplet[1]].sequence_lengths[i] == 0:
                    s2 = 0
                else:
                    s2 = 1
                if self.otus[triplet[2]].sequence_lengths[i] == 0:
                    s3 = 0
                else:
                    s3 = 1
                if (s1 * s2 * s3) == 1:
                    decisive_triples += 1
                    break
                i += 1
        return round(decisive_triples/float(total_triples), 2)




class Otu(object):
    """
    Class responsible for managing the data of each OTU in the supermatrix
    """

    
    name = ""               # the name of the OTU
    sequence = ""           # the full aligned sequence for this OTU in the supermatrix
    accessions = []         # list of each GenBank accession used where "-" means no sequence for that region
    sequence_lengths = []   # list of the length of each sequence



    def __init__(self, name):
        """
        Takes as input the name of the OTU
        """
        self.name = name
        self.sequence = ""
        self.accessions = []
        self.sequence_lengths = []



    def update(self, sequence, accession, sequence_length):
        """
        Inputs the aligned sequence (gaps already added), the accession #, and the unaligned sequence length.
        """
        self.sequence = self.sequence + sequence
        self.accessions.append(accession)
        self.sequence_lengths.append(sequence_length)



    def print_data(self):
        color = Color()
        print(color.blue + "Name = " + color.red + self.name)
        print(color.blue + "Sequence = " + color.red + self.sequence)
        print(color.blue + "Accessions = "  + color.red)
        print(self.accessions)
        print(color.blue + "Sequence_lengths = " + color.red)
        print(self.sequence_lengths)



