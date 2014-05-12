#! /usr/bin/python

import sys

class Color:
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
