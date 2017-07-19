#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
    Samir Jain, Eric Epstein, Trevor Klemp, Maggie Gray, Selman Jawed, Derek
    Braun* (*derek.braun@gallaudet.edu)
    
    Performs statistical analyses comparing data files created by simulator.py.
    Last updated: 11-Jul-2017 by Maggie Gray
'''

import argparse
import numpy
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
    print '** Shapiro-Wilk test of normality **'.format(filename)
    print '{:35}   {:^5}    {:^15}'.format('filename', 'p', 'mean ± stdev')
    data_array = []
    for e in experiments:
        # Find the mean and stdev, and run the shapiro-wilk test
        X = e.select(args.field)[-1]
        data_array.append(X)
        mean_stdev = "{:.4f} ± {:.4f}".format(numpy.mean(X), numpy.std(X))
        w, p = stats.shapiro(X)
        print '{:35}   {:.3f}   {}'.format (e.filename, p, mean_stdev)
        
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
        print '** Kruskal-Wallis H test **'
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
                        matrix += '{:^5}   '.format('-')
                    else:
                        U, p = stats.mstats.mannwhitneyu(data_array[i], data_array[j])
                        matrix += '{:.3f}   '.format(p)
                matrix += '\n'
            print matrix
        
    print 'Done. '



