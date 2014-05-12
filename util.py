#! /usr/bin/python

import gzip
from ftplib import FTP

class GenBank():
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
    def download_gb_db(cls, division):
        """
        Downloads and uncompresses files for a GenBank division.
        """
        print(Color.purple + "Connecting to ftp.ncbi.nih.gov..." + Color.done)
        ftp = FTP("ftp.ncbi.nih.gov")
        ftp.login()
        print(Color.yellow + "Opening directory genbank..." + Color.done)
        ftp.cwd("genbank")
        file_list = ftp.nlst()
        i = 1
        file_name = "gb" + division + str(i) + ".seq.gz"
        if not os.path.exists("genbank"):
            os.makedirs("genbank")
        while file_name in file_list:
            print(Color.red + "Downloading file " + file_name + Color.done)
            file = open("genbank/" + file_name, "wb")
            cls.getbinary(ftp, file_name, file)
            file.close()
            print(Color.yellow + "Uncompressing file " + file_name + Color.done)
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

    @staticmethod
    def setup_sqlite():
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



class Color:
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
