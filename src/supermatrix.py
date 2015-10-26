"""
SUMAC: supermatrix constructor

Copyright 2014 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


import os
import sys
import csv
import multiprocessing
from collections import OrderedDict
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from util import Color



class Supermatrix(object):
    """
    Class responsible for managing the final supermatrix
    """


    file = ""           # FASTA file of final supermatrix
    otus = {}           # dictionary of Otu objects
    pd = None           # fraction of triples, a measure of the partial decisiveness of a supermatrix
    missing_data = None # % of sequence data missing
    loci = None


    def __init__(self, alignments=None):
        """
        Optionally accept an Alignment object.
        """
        self.file = ""
        self.otus = {}
        self.pd = None
        self.coverage_density = None
        self.loci == None
        self.highest_OTU_decisiveness_score = 0
        self.highest_locus_decisiveness_score = 0
        self.lowest_OTU_decisiveness_score = 0
        self.lowest_locus_decisiveness_score = 0
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
                if alignments.user_provided:
                    otu = record.description
                else:
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
            records_for_loci = []
            for record in records:
                if alignments.user_provided:
                    otu = record.description
                else:
                    descriptors = record.description.split(" ")
                    otu = descriptors[1] + " " + descriptors[2]
                if otu not in already_added:
                    records_for_loci.append(SeqRecord(record.seq, id=otu, description=""))
                    if alignments.user_provided:
                        otus[otu].update(record.seq, alignment, self.get_ungapped_length(record.seq))
                    else:
                        otus[otu].update(record.seq, descriptors[0], self.get_ungapped_length(record.seq))
                    already_added.append(otu)
                loci_length = len(record.seq)
            total_length += loci_length
            
            # add '?' for any OTU that didn't have a sequence
            for otu in otus:
                if len(otus[otu].sequence) < total_length:
                    missing_seq = self.make_missing(loci_length)
                    records_for_loci.append(SeqRecord(Seq(missing_seq), id=otu, description=""))
                    otus[otu].update(missing_seq, "-", 0)
            
            # make fasta file of this loci
            alignment_out = alignment.split("/")
            f = open("alignments/supermatrix_" + alignment_out[len(alignment_out)-1], "w")
            sorted_records = sorted(records_for_loci, key=lambda record: record.id)
            SeqIO.write(sorted_records, f, "fasta")
            f.close()

        # order otus
        otus = OrderedDict(sorted(otus.items(), key=lambda t: t[0]))

        # write to FASTA file
        f = open("alignments/supermatrix_concatenated.fasta", "w")
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
        self.file = "alignments/supermatrix_concatenated.fasta"
        self.otus = otus



    def make_missing(self, length):
        """
        Inputs an integer.
        Returns a string of '-' of length
        """
        gap = ""
        for i in range(length):
            gap = "?" + gap
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
        total_missing = 0
        for record in records:
            otu = record.description
            missing = 0
            for letter in record.seq:
                if letter == '?':
                    missing += 1
                    total_missing += 1
            print(color.yellow + "OTU: " + color.red + otu + color.yellow + " % missing data = " + color.red + str(round(missing/float(len(record.seq)), 2)))
            num_records += 1
            matrix_length = len(record.seq)
        print(color.blue + "Total number of OTUs = " + color.red + str(num_records))
        print(color.blue + "Total length of matrix = " + color.red + str(matrix_length))
        print(color.blue + "Taxon coverage density = " + color.red + str(self.get_coverage_density()))
        print(color.blue + "Total % missing data = " + color.red + str(round(total_missing/float(matrix_length * num_records), 2)) + color.done)
        #for otu in self.otus: 
        #    self.otus[otu].print_data()



    def print_PD(self):
        """
        Prints partial decisiveness.
        """
        color = Color()
        print(color.blue + "Partial decisiveness (fraction of triples) = " + color.red + str(self.get_PD()) + color.done)



    def get_coverage_density(self):
        """
        Method to calculate the taxon coverage density.
        """
        if self.coverage_density != None:
            return self.coverage_density
        else:
            total_seq = 0
            missing_seq = 0
            for otu in self.otus:
                for accession in self.otus[otu].accessions:
                    total_seq += 1
                    if accession == "-":
                        missing_seq += 1
            self.coverage_density = round((total_seq - missing_seq)/float(total_seq), 2)
            return self.coverage_density



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
                row = [otu]
                for accession in self.otus[otu].accessions:
                    row.append(accession)
                csvwriter.writerow(row)
    
    
    
    def make_decisiveness_csv(self):
        """
        Method to generate a CSV file with the sequence decisiveness scores for all missing data.
        """
        with open('missing_sequence_decisiveness.csv', 'wb') as csv_output:
            csvwriter = csv.writer(csv_output)
            header = ["OTU", "Gene Region", "Sequence Decisiveness Score"]
            csvwriter.writerow(header)
            for otu in self.otus:
                i = 0
                for seq in self.otus[otu].sequence_lengths:
                    if seq == 0:
                        row = [self.otus[otu].name]
                        row.append(str(i+1))
                        score = self.calculate_sequence_decisiveness_score(self.otus[otu].decisiveness_score, self.loci[i][2])
                        row.append(str(score))
                        csvwriter.writerow(row)
                    i += 1



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



    def make_sequence_data_figure(self):
        
        # TODO:
        # test for numpy and matplotlib, show message if not present
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
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
        
        # set font size to scale with number of otus
        n = len(self.otus)
        if n < 25:
            font_size = 10
        else:
            font_size = round(11.8-(0.073*n), 1)
            if font_size < 0:
                font_size = 0

        supermatrix = np.array(data)
        
        # set up figure
        fig, ax = plt.subplots()
        ax.set_position([0.4, 0.1, .3, .7])
        # fig.set_size_inches(8.5, 11)

        # add data
        heatmap = ax.pcolor(supermatrix, cmap=plt.cm.Blues, alpha=0.8, vmin=0, vmax=100)

        # put the labels in the middle of each cell
        ax.set_yticks(np.arange(supermatrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(supermatrix.shape[1]) + 0.5, minor=False)

        # move x axis label stuff to top
        ax.invert_yaxis()
        ax.xaxis.tick_top()

        # add the labels
        ax.set_xticklabels(genes, minor=False, family="Arial", size=6)
        ax.set_yticklabels(otu_names, minor=False, family="Arial", size=font_size)

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

        plt.savefig("sequence_data_plot.pdf")



    def make_sequence_decisiveness_figure(self):
        """
        Method to make plot for the Sequence Decisiveness Scores.
        """

        data_all = []
        data_missing = []
        otu_names = []
        genes = []
        i = 1
        for otu in self.otus:
            seq_dec_scores_all = []
            seq_dec_scores_missing = []
            for locus in self.loci:
                score = self.calculate_sequence_decisiveness_score(self.otus[otu].decisiveness_score, self.loci[locus][2])
                # if sequence data present, make score 0
                if self.otus[otu].sequence_lengths[locus] != 0:
                    seq_dec_scores_missing.append(-0.25)
                    #seq_dec_scores_missing.append(-1)
                else:
                    seq_dec_scores_missing.append(score)
                seq_dec_scores_all.append(score)
            data_all.append(seq_dec_scores_all)
            data_missing.append(seq_dec_scores_missing)
            otu_names.append(otu)
            genes.append(str(i))
            i += 1
        
        # set font size to scale with number of otus
        n = len(self.otus)
        if n < 25:
            font_size = 10
        else:
            font_size = round(11.8-(0.073*n), 1)
            if font_size < 0:
                font_size = 0

        self.finish_sequence_decisiveness_figure(data_missing, font_size, otu_names, genes, "sequence_decisiveness_scores_plot.pdf")
        #self.finish_sequence_decisiveness_figure(data_all, font_size, otu_names, genes, "sequence_decisiveness_plot_all_scores.pdf")



    def finish_sequence_decisiveness_figure(self, data, font_size, otu_names, genes, file_name):
        """
        Method to finish plots for the Sequence Decisiveness Scores.
        """
        # TODO:
        # test for numpy and matplotlib, show message if not present
        import numpy as np
        import matplotlib
        #matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        supermatrix = np.array(data)

        # set up figure
        fig, ax = plt.subplots()
        ax.set_position([0.4, 0.1, .3, .7])
        # fig.set_size_inches(8.5, 11)

        # add data
        heatmap = ax.pcolor(supermatrix, cmap=plt.cm.YlOrRd, vmin=-0.25, vmax=1)
        #heatmap = ax.pcolor(supermatrix, cmap=plt.cm.Spectral_r, vmin=-1, vmax=1)

        # put the labels in the middle of each cell
        ax.set_yticks(np.arange(supermatrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(supermatrix.shape[1]) + 0.5, minor=False)

        # move x axis label stuff to top
        ax.invert_yaxis()
        ax.xaxis.tick_top()

        # add the labels
        ax.set_xticklabels(genes, minor=False, family="Arial", size=6)
        ax.set_yticklabels(otu_names, minor=False, family="Arial", size=font_size)

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
        cbar = fig.colorbar(heatmap, cax=axcolor, ticks=[-0.25, 0, 1], orientation='horizontal')
        cbar.ax.set_xticklabels(['-0.25', '0', '1'], family="Arial", size=6)
        #cbar = fig.colorbar(heatmap, cax=axcolor, ticks=[-1, 0, 1], orientation='horizontal')
        #cbar.ax.set_xticklabels(['-1', '0', '1'], family="Arial", size=6)

        plt.savefig(file_name)



    def get_PD(self):
        """
        method to get pd of supermatrix, calculating it if need be.
        """
        if self.pd == None:
            self.pd = self.calculate_PD()
            self.calculate_OTU_decisiveness_scores()
            self.calculate_locus_decisiveness_scores()
        return self.pd



    def calculate_PD(self):
        """
        Method to calculate the fraction of triples, a measure of partial decisiveness (PD).
        See: Sanderson, M.J., McMahon, M.M. & Steel, M., 2010. BMC evolutionary biology, 10. 
        """
        color = Color()
        decisive_triples = 0
        total_triples = 0
        total = self.binomial_coefficient(len(self.otus), 3)
        i = 0
        # nested loops to run through every possible triplet
        for otu1 in self.otus:
            if i < (len(self.otus) - 2):
                j = 0
                for otu2 in self.otus:
                    if i < j:
                        k = 0
                        for otu3 in self.otus:
                            if j < k:
                                sys.stdout.write("\r" + color.blue + "Calculating PD: " + color.red + str(round(100 * total_triples/float(total), 4)) + \
                                    "% " + color.blue + "finished" + color.done)
                                sys.stdout.flush()
                                # do PD calculations for this triplet
                                triplet = [otu1, otu2, otu3]
                                decisive, decisive_loci = self.calculate_triplet_PD(triplet)
                                decisive_triples += decisive
                                total_triples += 1
                                # triplet calculations for OTU and loci decisiveness scores
                                self.update_OTU_decisiveness(triplet, decisive)
                                self.update_locus_decisiveness(decisive_loci, decisive)
                            k += 1
                    j += 1
            i += 1
        sys.stdout.write("\r" + color.blue + "Calculating PD: " + color.red + "100.00% " + color.blue + "finished\n" + color.done)
        sys.stdout.flush()
        return round(decisive_triples/float(total_triples), 2)
    


    def calculate_PD_parallel(self, num_cores):
        """
        Method to calculate the fraction of triples, a measure of partial decisiveness (PD).
        See: Sanderson, M.J., McMahon, M.M. & Steel, M., 2010. BMC evolutionary biology, 10. 
        """
        color = Color()
        lock = multiprocessing.Lock()
        manager = multiprocessing.Manager()
        #already_compared = manager.list()
        #dist_matrix = manager.list()
        otus_shared = manager.dict()
        otus_shared = otus
        decisive_triples = manager.Value('i', 0)
        total_triples = manager.Value('i', 0)
        total = self.binomial_coefficient(len(self.otus), 3)

        for i in range(num_cores):
            p = multiprocessing.Process(target=calculate_PD_worker, args=(lock, i, num_cores, decisive_triples, total_triples, total, otus_shared))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()
        
        self.otus = otus_shared

        sys.stdout.write("\r" + color.blue + "Calculating PD: " + color.red + "100.00% " + color.blue + "finished\n" + color.done)
        sys.stdout.flush()
        return round(decisive_triples/float(total_triples), 2)



    def calculate_PD_worker(self, lock, process_num, num_cores, decisive_triples, total_triples, total, otus):
        """
        Worker function to help calculate fraction of triples
        """
        self.otus = otus
        # nested loops to run through every possible triplet
        i = process_num
        for otu1 in self.otus:
            if i < (len(self.otus) - 2):
                j = 0
                for otu2 in otus:
                    if i < j:
                        k = 0
                        for otu3 in otus:
                            if j < k:
                                sys.stdout.write("\r" + color.blue + "Calculating PD: " + color.red + str(round(100 * total_triples/float(total), 4)) + \
                                    "% " + color.blue + "finished" + color.done)
                                sys.stdout.flush()
                                # do PD calculations for this triplet
                                triplet = [otu1, otu2, otu3]
                                decisive, decisive_loci = self.calculate_triplet_PD(triplet)
                                with lock:
                                    decisive_triples += decisive
                                    total_triples += 1
                                # triplet calculations for OTU and loci decisiveness scores
                                self.update_OTU_decisiveness(triplet, decisive)
                                self.update_locus_decisiveness(decisive_loci, decisive)
                            k += 1
                    j += 1
            i += num_cores



    def calculate_triplet_PD(self, triplet):
        """
        Function that returns 1 if the triplet contains at least
        one gene region with sequence data for all three OTUs.
        Otherwise returns 0. Also return the list of loci that may be decisive.
        """
        decisive = 0
        i = 0
        decisive_loci = []
        while i < len(self.otus[triplet[0]].sequence_lengths):
            if self.otus[triplet[0]].sequence_lengths[i] == 0:
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
                decisive = 1
                decisive_loci.append(i)
            i += 1
        return decisive, decisive_loci

    # TODO: make parallel version of the following three functions:

    def update_OTU_decisiveness(self, triplet, decisive):
        """
        Function that loops through all OTUs not in the triplet, updating
        other_decisive_triples and other_total_triples for each OTU.
        PD without OTU = other_decisive_triples / other_total_triples
        OTU Decisiveness Score = PD / PD without OTU
        """
        for otu in self.otus:
            if self.otus[otu].other_decisive_triples == None:
                self.otus[otu].other_decisive_triples = 0
            if self.otus[otu].other_total_triples == None:
                self.otus[otu].other_total_triples = 0
            if otu not in triplet:
                self.otus[otu].other_decisive_triples += decisive
                self.otus[otu].other_total_triples += 1



    def update_locus_decisiveness(self, decisive_loci, decisive):
        """
        Function that loops through all loci not making this triplet
        decisive and updating other_decisive_triples and other_total_triples
        for each loci.
        PD without loci = other_decisive_triples / other_total_triples
        Locus Decisiveness Score = PD / PD without locus
        """
        if self.loci == None:
            self.setup_loci()
        for locus in self.loci:
            if locus not in decisive_loci:
                self.loci[locus][0] += decisive
                self.loci[locus][1] += 1



    def setup_loci(self):
        """
        Sets up a dictionary of the supermatrix's loci. Each dictionary value
        is a list [other_decisive_triples, other_total_triples]
        """
        self.loci = {}
        for otu in self.otus:
            i = 0
            while i < len(self.otus[otu].sequence_lengths):
                self.loci[i] = [0, 0]
                i += 1
            break



    def calculate_OTU_decisiveness_scores(self):
        """
        Method to calculate OTU Decisiveness Scores.
        """
        self.highest_OTU_decisiveness_score = 0
        self.lowest_OTU_decisiveness_score = 999
        for otu in self.otus:
            if self.otus[otu].other_total_triples != 0:
                PD_otu = self.otus[otu].other_decisive_triples/float(self.otus[otu].other_total_triples)
            else:
                PD_otu = 0
            if PD_otu != 0:
                score = self.pd/PD_otu
            else:
                PD_otu = 1/float(self.binomial_coefficient(len(self.otus),3))
                score = self.pd/PD_otu
            self.otus[otu].decisiveness_score = score
            if score > self.highest_OTU_decisiveness_score:
                self.highest_OTU_decisiveness_score = score
            if score < self.lowest_OTU_decisiveness_score:
                self.lowest_OTU_decisiveness_score = score



    def calculate_locus_decisiveness_scores(self):
        """
        Method to calculate Locus Decisiveness Scores.
        Adds the score to the list that is the self.loci dictionary value.
        """
        self.highest_locus_decisiveness_score = 0
        self.lowest_locus_decisiveness_score = 999
        for locus in self.loci:
            if self.loci[locus][1] != 0:
                PD_locus = self.loci[locus][0]/float(self.loci[locus][1])
            else:
                PD_locus = 0
            if PD_locus != 0:
                score = self.pd/PD_locus
            else:
                PD_locus = 1/float(self.binomial_coefficient(len(self.otus), 3))
                score = self.pd/PD_locus
            self.loci[locus].append(score)
            if score > self.highest_locus_decisiveness_score:
                self.highest_locus_decisiveness_score = score
            if score < self.lowest_locus_decisiveness_score:
                self.lowest_locus_decisiveness_score = score



    def calculate_sequence_decisiveness_score(self, otu_score, locus_score):
        """
        Method to calculate sequence decisiveness score.
        This is the [(OTU score - min OTU score)/(highest OTU score - min OTU score) + (locus score - min locus score)/(highest locus score - min locus)]/2
        """
        if self.highest_OTU_decisiveness_score != self.lowest_OTU_decisiveness_score:
            otu_relative_score = (otu_score - self.lowest_OTU_decisiveness_score)/float(self.highest_OTU_decisiveness_score - self.lowest_OTU_decisiveness_score)
        else:
            otu_relative_score = otu_score/float(self.highest_OTU_decisiveness_score)
        if self.highest_locus_decisiveness_score != self.lowest_locus_decisiveness_score:
            locus_relative_score = (locus_score - self.lowest_locus_decisiveness_score)/float(self.highest_locus_decisiveness_score - self.lowest_locus_decisiveness_score)
        else:
            locus_relative_score = locus_score/float(self.highest_locus_decisiveness_score)
        return round((otu_relative_score + locus_relative_score)/2, 3)



    def binomial_coefficient(self, n, k):
        """
        A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
        """
        if 0 <= k <= n:
            ntok = 1
            ktok = 1
            for t in xrange(1, min(k, n - k) + 1):
                ntok *= n
                ktok *= t
                n -= 1
            return ntok // ktok
        else:
            return 0



class Otu(object):
    """
    Class responsible for managing the data of each OTU in the supermatrix
    """

    
    name = ""               # the name of the OTU
    sequence = ""           # the full aligned sequence for this OTU in the supermatrix
    accessions = []         # list of each GenBank accession used where "-" means no sequence for that region
    sequence_lengths = []   # list of the length of each sequence
    other_decisive_triples = None
    other_total_triples = None



    def __init__(self, name):
        """
        Takes as input the name of the OTU
        """
        self.name = name
        self.sequence = ""
        self.accessions = []
        self.sequence_lengths = []
        self.other_decisive_triples = None
        self.other_total_triples = None



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



