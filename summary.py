#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
Samir Jain, Eric Epstein, Maggie Gray, Derek Braun*
(*derek.braun@gallaudet.edu)

Display a X-Y table showing final median values from data files created by
simulator.py.

Last updated: 15-May-2019 by Derek Braun
'''

import sys
import os
import argparse
import numpy
import random
import fileio
from scipy import stats

#
#   MAIN ROUTINE
#
if __name__ == '__main__':
    # reading from the files
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filenames', nargs='+',
                        help = 'filenames for data files.')
    parser.add_argument('-f', '--field', action='store',
                        default = 'a',
                        help = 'the field to compare across simulations.')
    parser.add_argument('-x', action='store',
                        default = 'aa_fitness',
                        help = 'the x variable to compare across simulations.')
    parser.add_argument('-y', action='store',
                        default = 'aa_homogamy',
                        help = 'the y variable to compare across simulations.')
    args=parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    experiments = []
    for filename in args.filenames:
        # Check to see if each individual file exists
        if os.path.isfile(filename):
            experiments.append(fileio.Experiment(filename))
            print 'Reading {}'.format(filename)
        else:
            print "File {} not found.".format(filename)
            exit()

    # Check to see if the field is a valid one
    for e in experiments:
        if not hasattr(e, args.x):
            print '"{}" is not a metadata header in {}.'\
                  ''.format(args.x, e.filename)
            exit()
        if not hasattr(e, args.y):
            print '"{}" is not a metadata header in {}.'\
                  ''.format(args.y, e.filename)
            exit()

    l = []
    for e in experiments:
        l += [(getattr(e, args.y), getattr(e, args.x),
              numpy.median(e.select(args.field)[-1]))]
    table = sorted(l)

    # now print the damn table
    print
    print '** Medians for "{}" as a function of "{}" and "{}" **'.format(args.field, args.x, args.y)
    # print X headers
    #print '{:30}   '
    #for t in table:
    #    print '{:^8}   '.format(t[0]),
    #print

    y0 = None
    for (y, x, M) in table:
        if y == y0:
            print '{:0.6f}  '.format(M),
        else:
            y0 = y
            print
            print '{:^6}  '.format(y),
            continue
    print
    print 'Done.'
    exit()
