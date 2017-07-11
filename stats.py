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
    parser.add_argument('--fields', action='store',
                        help = 'the variables to compare among populations. a, aa, or F.')
    parser.add_argument('-c','--shapiro_wilk', action='store_true',
                        default=True,
                        help = 'Run the Shapiro-Wilk test for normality')
    args=parser.parse_args()
                        
    # Checking to see if the field is a valid one
    if not(args.fields == "a" or args.fields == "aa" or args.fields == "F"):
        print '  field not valid. Please try again with "a", "aa", or "F".\n'
        exit()                        

    for file in args.filenames:
    	
    	# Checking to see if each individual file is actually an existing file
    	if os.path.isfile(file):
    		e = fileio.Experiment(file)
    		x = e.select(args.fields)[-1]	# The array of data at the final point of time
    		print '\n  Reading {}'.format(file)
    	else:
    		print "\n  File {} not found.".format(file)
    		exit()
    		     
        # Finding the mean and stdev
        f_mean_stdev = "{} ± {}".format(numpy.mean(x), numpy.std(x))
        
        # Running the shapiro wilk test for each individual file
        if args.shapiro_wilk:
            w, normality_pval = stats.shapiro(x)
            print "  ** Ran test of normality for [{}] **".format(file) # check
            print "  Shapiro-Wilk p-value:  {}".format(normality_pval)
            print "  mean ± standard deviation:  {} ".format(f_mean_stdev)
    
    # Creating a list of the numbers to compare (in preparation for the statistical test)
    fieldsList = []
    i = 0
    while i < len(args.filenames):
    	f = fileio.Experiment(args.filenames[i]).select(args.fields)[-1]
        fieldsList.append(f)
        i += 1
        
    # Running the tests
    if len(args.filenames) == 2:
        statistic, pval = stats.mstats.mannwhitneyu(fieldsList[0], fieldsList[1])
        print "\n  ** Ran comparisons **"
        print "  filenames:  {}".format(args.filenames)
        print "  Mann-Whitney U:  {}".format(statistic)
        print "  p-value:  {}".format(pval)
    elif len(args.filenames) > 2:
    	
        statistic, pval = stats.mstats.kruskal(fieldsList)
        print "\n  ** Ran comparisons **"
        print "  filenames:  {}".format(args.filenames)
        print "  Kruskal-Wallis H:  {}".format(statistic)
        print "  p-value:  {}\n".format(pval)
    
    print 'Done. '

# To run: write $ ./stats.py *.tsv --fields "aa"

