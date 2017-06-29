#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
    Samir Jain, Eric Epstein, Trevor Klemp, Maggie Gray, Selman Jawed, Derek 
    Braun* (*derek.braun@gallaudet.edu)
    
    Performs statistical analyses comparing data files created by Simulator.py.
    Last updated: 23-Jun-2017 by Derek Braun
'''

import argparse
import numpy
import fileio
import os


def foo:
    print 'bar'
    
    
#
#   MAIN ROUTINE
#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename',
                        help = 'filename for data file')
    parser.add_argument('-s','--sample arg',action='store_true',
                        help = 'sample argument')
    args=parser.parse_args()
    
    if os.path.isfile(args.filename):
        e = fileio.Experiment(args.filename)
        print '  Reading {}'.format(args.filename)
    else:
        print '  File not found.'