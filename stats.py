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
    parser.add_argument('filename_one',
                        help = 'filename for first data file')
    parser.add_argument('filename_two',
                        help = 'filename for second data file')
    parser.add_argument('-s','--mann_whitney',action='store_true',
    					
                        help = 'perform the mann-whitney tests and get stats')
    args=parser.parse_args()
    
    if os.path.isfile(args.filename_one) and os.path.isfile(args.filename_two):
        eOne = fileio.Experiment(args.filename_one)
        eTwo = fileio.Experiment(args.filename_two)
        print '  Reading {}'.format(args.filename_one) + '& Reading {}'.format(args.filename_one)
    elif os.path.isfile(args.filename_two) and not(os.path.isfile(args.filename_one)):
        print '  File one not found.'	
    elif os.path.isfile(args.filename_one) and not(os.path.isfile(args.filename_two)):
        print '  File two not found.'
    else:
    	print '  Both files not found.'
        
    if args.mann_whitney
    	