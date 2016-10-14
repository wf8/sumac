"""
SUMAC: supermatrix constructor

Copyright 2014 Will Freyman - freyman@berkeley.edu
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""


import os
import sys
import gzip
import pickle
from Bio import Entrez
from Bio import SeqIO
from ftplib import FTP
from util import Color


class GenBankSetup(object):
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
    def download(cls, divisions, path):
        """
        Downloads and uncompresses files for a GenBank division.
        Path should be the absolute path to save the GB files.
        """
        for division_input in divisions:
            color = Color()
            division = str(division_input).lower()
            print(color.purple + "Connecting to ftp.ncbi.nlm.nih.gov..." + color.done)
            ftp = FTP("ftp.ncbi.nlm.nih.gov")
            ftp.login()
            print(color.yellow + "Opening directory genbank..." + color.done)
            ftp.cwd("genbank")
            file_list = ftp.nlst()
            i = 1
            file_name = "gb" + division + str(i) + ".seq.gz"
            if not os.path.exists(path):
                os.makedirs(path)
            path = path + "/"
            while file_name in file_list:
                print(color.red + "Downloading file " + file_name + color.done)
                file = open(path + file_name, "wb")
                cls.getbinary(ftp, file_name, file)
                file.close()
                print(color.yellow + "Uncompressing file " + file_name + color.done)
                file = gzip.open(path + file_name, "rb")
                file_content = file.read()
                file.close()
                file = open(path + file_name[:-3], "wb")
                file.write(file_content)
                file.close()
                os.remove(path + file_name)
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
    def sqlite(path):
        """
        Sets up the SQLite db for the GenBank division.
        Path is the absolute path of the GB files.
        Returns a dictionary of SeqRecord objects.
        """
        color = Color()
        if os.path.exists(path + "/gb.idx"):
            print(color.purple + "Genbank database already downloaded. Indexing sequences..." + color.done)
            return SeqIO.index_db(path + "/gb.idx")
        else:
            files = os.listdir(path)
            path_files = []
            if len(files) == 0:
                print(color.red + "GenBank files not found. Re-download with the -d option. See --help for more details." + color.done)
                sys.exit(0)
            for file in files:
                path_files.append(path + "/" + file)
            print(color.purple + "Genbank database already downloaded. Indexing sequences..." + color.done)
            return SeqIO.index_db(path + "/gb.idx", path_files, "genbank")




class GenBankSearch(object):
    """
    Class responsible for searching GenBank and managing lists of keys 
    to all sequences in ingroup and outgroup.
    Path is the absolute path of the GB files.
    """

    ingroup = ''
    ingroup_keys = []
    outgroup = ''
    outgroup_keys = []
    path = ''

    def __init__(self, gb, ingroup, outgroup, max_ingroup=None):
        """
        Takes as input a dictionary of SeqRecords gb and the names of ingroup 
        and outgroup clades. Option parameter location of directory to save GB files.
        Finds lists of keys to SeqRecords for the ingroup and outgroup.
        """
        self.ingroup = ingroup
        self.outgroup = outgroup
        self.ingroup_keys = []
        self.outgroup_keys = []

        # check to see if this ingroup and outgroup have already been found
        if self.check_for_results(): 
            self.read_file()
        else:
            self.search(gb, max_ingroup)


    def search(self, gb, max_ingroup):
        """
        Perform search of all GB SeqRecords for ingroup/outgroup,
        and save results of search to file.
        """
        keys = gb.keys()
        total = len(keys)
        i = 0
        for key in keys:
            if self.ingroup in gb[key].annotations['taxonomy']:
                self.ingroup_keys.append(key)
            elif self.outgroup in gb[key].annotations['taxonomy']:
                self.outgroup_keys.append(key)
            self.print_search_status(i, total)
            i += 1
            if max_ingroup is not None and len(self.ingroup_keys) == max_ingroup:
                sys.stdout.write("\n")
                sys.stdout.flush()  
                self.write_file()
                return 
        sys.stdout.write("\n")
        sys.stdout.flush()
        self.write_file()


    def print_search_status(self, i, total):
        color = Color()
        sys.stdout.write('\r' + color.yellow + 'Ingroup sequences found: ' \
                          + color.red + str(len(self.ingroup_keys)) + color.yellow \
                          + '  Outgroup sequences found: ' + color.red \
                          + str(len(self.outgroup_keys)) + color.yellow \
                          + '  Percent searched: ' + color.red \
                          + str(round( 100 * float(i) / total , 1)) + color.done
                        )
        sys.stdout.flush()


    def write_file(self):
        """
        Saves results of GB search to file.
        """
        data = {"ingroup": self.ingroup, "outgroup": self.outgroup, "ingroup_keys": self.ingroup_keys, "outgroup_keys": self.outgroup_keys}
        out = open( "gb_search_results", "w" )
        pickle.dump(data, out)
        out.close()


    def read_file(self):
        """
        Loads results of GB search from file.
        """
        groups = pickle.load( open( "gb_search_results", "rb" ) )
        self.ingroup = groups["ingroup"]
        self.outgroup = groups["outgroup"]
        self.ingroup_keys = groups["ingroup_keys"]
        self.outgroup_keys = groups["outgroup_keys"]
        color = Color()
        print(color.yellow + 'This search was already performed. Loading previous results...' + color.done)
        print(color.yellow + 'Ingroup sequences found: ' \
              + color.red + str(len(self.ingroup_keys)) + color.yellow \
              + '  Outgroup sequences found: ' + color.red \
              + str(len(self.outgroup_keys)) + color.done \
             )


    def check_for_results(self):
        """
        Check to see if ingroup/outgroup sequences have already been found.
        """
        if not os.path.exists("gb_search_results"):
            return False
        else:
            groups = pickle.load( open( "gb_search_results", "rb" ) )
            if self.ingroup == groups["ingroup"] and self.outgroup == groups["outgroup"]:
                return True
