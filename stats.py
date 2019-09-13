#!/usr/local/bin/python3 -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
Samir Jain, Eric Epstein, Maggie Gray, Derek Braun*
(*derek.braun@gallaudet.edu)

Performs statistical analyses comparing data files created by simulator.py.

Last updated: 2-May-2019 by Derek Braun
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
                        required = True,
                        help = 'the variable to compare among populations. a, aa, or F.')
    args=parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    experiments = []
    print('Reading file(s)...')
    for filename in args.filenames:
        # Check to see if each individual file exists
        if os.path.isfile(filename):
            experiments.append(fileio.Experiment(filename))
            print('   {}'.format(filename))
        else:
            print('   File {} not found.'.format(filename))
            exit()

    # Check to see if the field is a valid one
    for e in experiments:
        if args.field not in e.headers:
            print('Field ""{}"" is not a header in {}.'\
                  ''.format(args.field, e.filename))
            exit()

    print()
    print('** Summary Statistics for "{}" **'.format(args.field))
    print('{:30}   {:^8}   {:^8}  {:^21}'.format('', 'start', 'end', ''))
    print('{:30}   {:^8}   {:^8}  {:^21}'.format('filename','(median)','(median)','95% CI'))
    for e in experiments:
        # select first and last set of values
        Xo = e.select(args.field)[0]
        X = e.select_endpoint(args.field)
        # Find the mean and stdev
        start_median = numpy.median(Xo)
        end_median = numpy.median(X)
        X.sort()
        ci = '({:0.6f} - {:0.6f})'.format(X[int(0.025*len(X))],
                                      X[int(0.975*len(X))])
        print('{:30}   {:0.6f}   {:0.6f}  {:^21}'.format(e.filename, start_median, \
                                             end_median, ci))
    print()
    print('** Shapiro-Wilk test of normality for "{}" **'.format(args.field, \
                                                                 filename))
    print('{:30}   {:^8}   {:^14}'.format('filename', 'p-value', 'interpretation'))
    data_array = []
    for e in experiments:
        # select last set of values
        X = e.select_endpoint(args.field)
        # Find the mean and stdev,
        data_array.append(X)
        # Thin X to 5,000 values if needed; the maximum for the Shapiro-Wilk
        # Then run the Shapiro-Wilk
        if len(X) > 5000:
            w, p = stats.shapiro(random.sample(X, 5000))
        else:
            w, p = stats.shapiro(X)
        print('{:30}   {:^8.3g}   {:^14}'.format (e.filename, p, 'not normal' \
                                                  if p < 0.05 else 'normal'))

    # Running the tests
    if len(experiments) == 2:
        print('')
        print('** Mann-Whitney U test for "{}" **'.format(args.field))
        print('{:5}   {:^8}'.format('U', 'p-value'))
        U, p = stats.mstats.mannwhitneyu(data_array[0], data_array[1])
        print("{:<5.1f}   {:^8.5g}  {:e}".format(U, p, p))
        print()
    else:
        print('')
        print('** Kruskal-Wallis one-way analysis of variance for "{}" **' \
              ''.format(args.field))
        print('{:^8}   {:^8}'.format('H', 'p-value'))
        H, p = stats.mstats.kruskal(*data_array)
        print("{:<8.1f}   {:^8.5g}".format(H, p))
        print()
        if p <= 0.5:
            print('** post-hoc pairwise Mann-Whitney U test for "{}" **' \
                  ''.format(args.field))
            matrix = ''
            for i in range(len(data_array)):
                for j in range(len(data_array)):
                    if j >= i:
                        matrix += '{:^8}   '.format('-')
                    else:
                        U, p = stats.mstats.mannwhitneyu(data_array[i],
                                                         data_array[j])
                        matrix += '{:^8.3g}   '.format(p)
                matrix += '\n'
            print(matrix)
    print('Done. ')
