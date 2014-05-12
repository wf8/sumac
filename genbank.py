#! /usr/bin/python

import os
import sys
import gzip
from Bio import Entrez
from Bio import SeqIO
from ftplib import FTP
from util import Color


class GenBankSetup():
    """
    Class responsible for downloading and indexing GenBank files.
    """

    @staticmethod
    def gettext(ftp, filename, outfile=None):
        """
        Fetch a text file.
        """
        if outfile is None:
            outfile = sys.stdout
        ftp.retrlines("RETR " + filename, lambda s, w=outfile.write: w(s+"\n"))

    @staticmethod
    def getbinary(ftp, filename, outfile=None):
        """
        Fetch a binary file.
        """
        if outfile is None:
            outfile = sys.stdout
        ftp.retrbinary("RETR " + filename, outfile.write)

    @classmethod
    def download(cls, division_input):
        """
        Downloads and uncompresses files for a GenBank division.
        """
        color = Color()
        division = str(division_input).lower()
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
            cls.getbinary(ftp, file_name, file)
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
        # check if any files were downloaded
        if i == 1:
            print(color.red + "GenBank division " + division_input \
                  + " not found. Please use a valid division name " \
                  + "(e.g. VRT, INV, PLN)." + color.done
                 )
            sys.exit(0)
        ftp.quit()

    @staticmethod
    def sqlite():
        """
        Sets up the SQLite db for the GenBank division.
        Returns a dictionary of SeqRecord objects.
        """
        gb_dir = os.path.abspath("genbank/")
        files = os.listdir(gb_dir)
        abs_path_files = []
        for file in files:
            abs_path_files.append(gb_dir + "/" + file)
        gb = SeqIO.index_db(gb_dir + "/gb.idx", abs_path_files, "genbank")
        return gb



class GenBankSearch:
    """
    Class responsible for searching GenBank and managing lists of keys 
    to all sequences in ingroup and outgroup.
    """

    ingroup_keys = []
    outgroup_keys = []

    def __init__(self, gb, ingroup, outgroup):
        """
        Takes as input a dictionary of SeqRecords gb and the names of ingroup 
        and outgroup clades.
        Finds lists of keys to SeqRecords for the ingroup and outgroup.
        """
        keys = gb.keys()
        self.ingroup_keys = []
        self.outgroup_keys = []
        total = len(keys)
        i = 0
        for key in keys:
            if ingroup in gb[key].annotations['taxonomy']:
                self.ingroup_keys.append(key)
            elif outgroup in gb[key].annotations['taxonomy']:
                self.outgroup_keys.append(key)
            self.print_search_status(i)
            i += 1
            ## FOR TESTING ONLY
            #if len(ingroup_keys) == 50:  ##
            #    sys.stdout.write("\n")   ## FOR TESTING ONLY
            #    sys.stdout.flush()       ## # FOR TESTING ONLY
            #    return ingroup_keys, outgroup_keys    # FOR TESTING ONLY
            ## FOR TESTING ONLY
            ## remove above
        sys.stdout.write("\n")
        sys.stdout.flush()
        return ingroup_keys, outgroup_keys
    
    def print_search_status(self, i):
        color = Color()
        sys.stdout.write('\r' + color.yellow + 'Ingroup sequences found: ' \
                          + color.red + str(len(self.ingroup_keys)) + color.yellow \
                          + '  Outgroup sequences found: ' + color.red \
                          + str(len(self.outgroup_keys)) + color.yellow \
                          + '  Percent searched: ' + color.red \
                          + str(round( 100 * float(i) / total , 1)) + color.done
                        )
        sys.stdout.flush()

    def write_keys_to_file():
        return True

    def read_keys_from_file():
        return True

