#!/usr/local/bin/python3 -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
*Derek C. Braun, Brian H. Greenwald, Samir Jain, Eric Epstein, Brienna Herold, Maggie Gray
*derek.braun@gallaudet.edu
Performs statistical analyses comparing data files created by simulator.py.
'''


import sys
import os
import argparse
import numpy
import random
import fileio
from scipy import stats


def interpret(p):
    if p > 0.05:
        return 'not significant'
    elif p > 0.01:
        return 'significant'
    else:
        return 'highly significant'
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
    print('{:30}   {:^8}  {:^21}'.format('', 'end', ''))
    print('{:30}   {:^8}  {:^21}'.format('filename','(median)','95% interval'))
    for e in experiments:
        # select first and last set of values
        X = e.select_endpoint(args.field)
        end_median = numpy.median(X)
        X.sort()
        ci = '({:0.6f} - {:0.6f})'.format(X[int(0.025*len(X))],
                                      X[int(0.975*len(X))])
        print('{:30}   {:0.6f}  {:^21}'.format(e.filename, end_median, ci))
    print()
    print('** Shapiro-Wilk test of normality for "{}" **'.format(args.field, \
                                                                 filename))
    print('{:30}   {:^7}   {:^18}'.format('filename', 'p', 'interpretation'))
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
        print('{:30}   {:^5.3G}   {:^18}'.format (e.filename, p, interpret(p)))

    # Running the tests
    if len(experiments) == 2:
        print('')
        print('** Mann-Whitney U test for "{}" **'.format(args.field))
        print('{:8}   {:8}   {:8}   {:^5}   {:^6}   {:^18}'.format('U1', 'U2', 'U', 'p', 'f1', 'interpretation'))
        n1, n2 = len(data_array[0]), len(data_array[1])
        U1, p = stats.mannwhitneyu(data_array[0], data_array[1], alternative = 'two-sided')  # U is U for y
        U2 = (n1 * n2) - U1
        U = min(U1, U2)
        f1 = U1/(n1 * n2)
        print("{:<8.3G}   {:<8.3G}   {:<8.3G}   {:^5.3G}   {:^5.2%}   {:^18}".format(U1, U2, U, p, f1, interpret(p)))
        print()
    else:
        print('')
        print('** Kruskal-Wallis one-way analysis of variance for "{}" **' \
              ''.format(args.field))
        print('{:^8}   {:^8}'.format('H', 'p-value'))
        H, p = stats.kruskal(*data_array)
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
                        U2, p = stats.mannwhitneyu(data_array[i],
                                                   data_array[j],
                                                   alternative = 'two-sided')
                        matrix += '{:^8.3g}   '.format(p)
                matrix += '\n'
            print(matrix)
    print('Done. ')
