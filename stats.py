#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
    Samir Jain, Eric Epstein, Maggie Gray, Derek Braun*
    (*derek.braun@gallaudet.edu)

    Performs statistical analyses comparing data files created by simulator.py.

    Last updated: 11-Jul-2017 by Maggie Gray
'''

import argparse
import numpy
import random
import fileio
import os
from scipy import stats

#
#   MAIN ROUTINE
#
if __name__ == '__main__':
    # reading from the files
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filenames', nargs='+',
                        help = 'filenames for data files')
    parser.add_argument('-f', '--field', action='store',
                        help = 'the variables to compare among populations. a, aa, or F.')
    args=parser.parse_args()

    # Check to see if the field is a valid one
    if args.field not in ['a','aa']:
        print 'Field {} is not valid. Please try again with "a" or "aa".'\
              ''.format(args.field)
        exit()

    experiments = []
    for filename in args.filenames:
        # Check to see if each individual file exists
        if os.path.isfile(filename):
            experiments.append(fileio.Experiment(filename))
            print 'Reading {}'.format(filename)
        else:
            print "File {} not found.".format(filename)
            exit()

    print
    print '{:30}   {:^8}   {:^8}  {:^21}'.format('** Summary Statistics **', 'start', \
                                                 'end', '')
    print '{:30}   {:^8}   {:^8}  {:^21}'.format('filename','(median)','(median)','95% CI')
    for e in experiments:
        # select first and last set of values
        Xo = e.select(args.field)[0]
        X = e.select(args.field)[-1]
        # Find the mean and stdev
        start_median = numpy.median(Xo)
        end_median = numpy.median(X)
        X.sort()
        ci = '({:0.6f} - {:0.6f})'.format(X[int(0.025*len(X))],
                                      X[int(0.975*len(X))])
        print '{:30}   {:0.6f}   {:0.6f}  {:^21}'.format(e.filename, start_median, \
                                             end_median, ci)
    print
    print '** Shapiro-Wilk test of normality **'.format(filename)
    print '{:30}   {:^5}   {:^14}'.format('filename', 'p', 'interpretation')
    data_array = []
    for e in experiments:
        # select last set of values
        X = e.select(args.field)[-1]
        # Find the mean and stdev,
        data_array.append(X)
        # Thin X to 5,000 values if needed; the maximum for the Shapiro-Wilk
        # Then run the Shapiro-Wilk
        if len(X) > 5000:
            w, p = stats.shapiro(random.sample(X, 5000))
        else:
            w, p = stats.shapiro(X)
        print '{:30}   {:0.3f}   {:^14}'.format (e.filename, p, 'not normal' if p < 0.05 else 'normal')

    # Running the tests
    if len(experiments) == 2:
        print ''
        print '** Mann-Whitney U test **'
        print '{:5}   {:3}'.format('U', 'p-value')
        U, p = stats.mstats.mannwhitneyu(data_array[0], data_array[1])
        print "{:<5.1f}   {:.3f}".format(U, p)
        print
    else:
        print ''
        print '** Kruskal-Wallis one-way analysis of variance **'
        print '{:^5}   {:^5}'.format('H', 'p')
        H, p = stats.mstats.kruskal(*data_array)
        print "{:<5.1f}   {:.3f}".format(H, p)
        print
        if p <= 0.5:
            print '** post-hoc pairwise Mann-Whitney U test **'
            matrix = ''
            for i in range(len(data_array)):
                for j in range(len(data_array)):
                    if j >= i:
                        matrix += '{:^7}   '.format('-')
                    else:
                        U, p = stats.mstats.mannwhitneyu(data_array[i], data_array[j])
                        matrix += '{:.5f}   '.format(p)
                matrix += '\n'
            print matrix

    print 'Done. '
