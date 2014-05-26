#! /usr/bin/python



import sys
import re



class Color(object):
    """
    Class for printing pretty colors to the terminal.
    """
    
    purple = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    done = '\033[0m'

    def __init__(self):
        # turn off colors if not running in terminal
        if not sys.stdout.isatty():
            self.disable()

    def disable(self):
        self.purple = ''
        self.blue = ''
        self.green = ''
        self.yellow = ''
        self.red = ''
        self.done = ''



class Logger(object):
    """
    Class for logging both to file and to terminal. Modified from:
    http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python/
    """

    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("sumac_log", "w")


    def write(self, message):
        self.terminal.write(message)
        trimmed_message = self.remove_ansi_colors(message)
        # remove updates flushed to terminal from log
        if len(trimmed_message.splitlines()) > 1:
            trimmed_message = ""
        self.log.write(trimmed_message)  


    def flush(self):
        self.terminal.flush()
        self.log.flush()


    def isatty(self):
        return self.terminal.isatty()


    def remove_ansi_colors(self, message):
        ansi_escape = re.compile(r'\x1b[^m]*m')
        return ansi_escape.sub('', message)
