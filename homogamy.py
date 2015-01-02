#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-

'''
    See: README.md for To-Do List.

    Samir Jain, Eric Epstein, Derek Braun*
    *derek.braun@gallaudet.edu
    Summer of '14
    
    We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/
    
    Because this code is intended for archival and may be re-run years down
    the road, we purposely wrote this code to be as self-sufficient as 
    possible. We have contained all our code in one script and avoided
    spilling over to external modules wherever reasonable.
'''

DEBUG_MODE = True

import sys
import os
import time
import random
import argparse
import shutil
import time
import csv
import multiprocessing
import numpy 
import matplotlib
# The matplotlib.use command, which specifies the backend, *must* appear
# before importing pyplot. For choices of backends, see:
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
matplotlib.use('Agg') 
from matplotlib import pyplot as plt
import simuOpt
if not DEBUG_MODE:
    simuOpt.setOptions(optimized=True, numThreads=0)
from simuPOP import *
print

PROPOSALS = 1000
#A_FREQ = 0.01304
A_FREQ = 0.1304
EXPERIMENTS = [ # small pop, equal fitness, random mating
                {'constant_pop_size': 10000,      
                'aa_fitness'        : 1,
                'aa_homogamy'       : 0.},
                # small pop, equal fitness, assortative mating
               {'constant_pop_size' : 10000,     
                'aa_fitness'        : 1,
                'aa_homogamy'       : 0.9},
                # small pop, aa homozygotes have 2x fitness, random mating
               {'constant_pop_size' : 10000,      
                'aa_fitness'        : 2,
                'aa_homogamy'       : 0.},
                # small pop, aa homozygotes have 2x fitness, assortative mating
               {'constant_pop_size' : 10000,      
                'aa_fitness'        : 2,
                'aa_homogamy'       : 0.9}]

GALLAUDET_BLUE = '#003b65'
GALLAUDET_BUFF = '#e5d19e'
PEN_COLOR = 'black'


def varyAssort(constant_pop_size, aa_fitness, aa_homogamy):
    '''
        Accepts:
        constant_pop_size   population size, which remains constant throughout
        aa_fitness          _relative_ fitness of deaf individuals
        aa_homogamy         the percent of assortative mating between
                            deaf individuals

        Returns a dict containing the results from the simulation.
    '''        
    setRNG(random.seed(getRNG().seed()))
    pop = Population(constant_pop_size, loci=1, infoFields='fitness')
    # The next line creates virtual subpopulations that are needed for 
    # defining the mating scheme: 
    #   The (0,0) subPop represents AA/Aa
    #   The (0,1) subPop represents aa
    vsps = GenotypeSplitter(loci=0, alleles=[[0,0,0,1],[1,1]])  
    pop.setVirtualSplitter(vsps)
    pop.evolve(
        initOps = [InitSex(),
                   InitGenotype(freq=[1-A_FREQ, A_FREQ]),
                   PyExec('header=[]'),
                   PyExec('row=[]')],
        preOps = Stat(popSize=True, alleleFreq=[0], subPops=[(0,0), (0,1)]),
        matingScheme = RandomMating(subPops=[(0,0), (0,1)]),
        postOps = [PyExec(r"header += ['gen','AA/Aa_size','aa_size','A_freq',"\
                                      "'a_freq']"),
                   PyExec(r"row += [gen,"\
                                   "subPopSize[0],"\
                                   "subPopSize[1],"\
                                   "alleleFreq[0][0], alleleFreq[0][1]]")],
        gen = 10,
    )
    pop.evolve(
        preOps = [Stat(popSize=True, alleleFreq=[0], subPops=[(0,0), (0,1)]),
                  # Because simuPop uses _absolute_ fitness and my logic uses
                  # relative fitness, I need to convert relative fitness into
                  # absolute fitness before passing to simuPOP.
                  MapSelector(loci=0, fitness={(0,0):1./aa_fitness,
                                               (0,1):1./aa_fitness,
                                               (1,1):1.})],
        # from documentation:        
        # If multiple mating schemes are applied to the same subpopulation, 
        # a weight (parameter weight) can be given to each mating scheme to 
        # determine how many offspring it will produce. The default for all 
        # mating schemes are 0. In this case, the number of offspring each 
        # mating scheme produces is proportional to the size of its parental 
        # (virtual) subpopulation. If all weights are negative, the numbers of 
        # offspring are determined by the multiplication of the absolute values 
        # of the weights and their respective parental (virtual) subpopulation 
        # sizes. If all weights are positive, the number of offspring produced 
        # by each mating scheme is proportional to these weights. Mating
        # schemes with zero weight in this case will produce no offspring. 
        # If both negative and positive weights are present, 
        # negative weights are processed before positive ones.
        matingScheme = HeteroMating([RandomMating(subPops=[(0,0), (0,1)],
                                                  weight = -1+aa_homogamy),
                                     RandomMating(subPops=[(0,0)],
                                                  weight = 0),
                                     RandomMating(subPops=[(0,1)],
                                                  weight = 0)]),
        postOps = [PyExec(r"header += ['gen','AA/Aa_size','aa_size','A_freq',"\
                                      "'a_freq']"),
                   PyExec(r"row += [gen,"\
                                   "subPopSize[0],"\
                                   "subPopSize[1],"\
                                   "alleleFreq[0][0], alleleFreq[0][1]]")],
        gen = 32,
    )
    return {'header':pop.dvars().header,
            'row':pop.dvars().row}


def contour_plot(ax, X, Yt, title=None, xlabel=None, ylabel=None,
                 titlesize=20, labelsize=20, ticklabelsize=16, scaling=1):
    '''
        Produces a type of line chart, where for each x, the median value is 
        shown as a line, and the area between the 5% and 95% CI are shaded.
        
        Accepts:
            ax              a matplotlib.pyplot axis instance
            X               an array of x values
            Yt              an array of y tuples (95%, median, 5%)
            xlabel          x axis label string
            ylabel          y axis label string
            title           axis title
            titlesize       font size (in points) for titles
            labelsize       font size (in points) for x and y labels
            ticklabelsize   font size (in points) for xtick and ytick labels
            scaling         a float which scales the graph, including line widths
                            and all font sizes.
    '''
    # unpack Y tuple
    Y_uppers = [t[0] for t in Yt]
    Y_medians = [t[1] for t in Yt]
    Y_lowers = [t[2] for t in Yt]
    ax.set_xlim(min(X),max(X))
    if title is not None:
        ax.set_title(title, fontsize=titlesize*scaling)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=labelsize*scaling)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=labelsize*scaling)
    ax.fill_between(X, Y_uppers, Y_lowers, color=GALLAUDET_BUFF)
    ax.text(max(X)*0.98, max(Y_uppers), '{:.2f}'.format(Y_uppers[-1]), 
            va='bottom', ha='right', fontsize=ticklabelsize*scaling)
    ax.plot(X, Y_medians, color=GALLAUDET_BLUE, lw=3*scaling)
    ax.text(max(X)*0.98, max(Y_medians), '{:.2f}'.format(Y_medians[-1]), 
            va='bottom', ha='right', fontsize=ticklabelsize*scaling)
    ax.text(max(X)*0.98, min(Y_lowers), '{:.2f}'.format(Y_lowers[-1]), 
            va='top', ha='right', fontsize=ticklabelsize*scaling)
    # turn off gridlines
    ax.grid(False)
    # set tick label sizes
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(ticklabelsize*scaling) 
    return ax   
    

def write_summary_contour_plot(filename, Xarr, Yarr, nrows, ncols, titles=[], 
                               title=None, xlabel=None, ylabel=None,
                               titlesize=24, labelsize=20, ticklabelsize=14):
    '''
        Produces a type of line chart, where for each x, the median value is 
        shown as a line, and the area between the 5% and 95% CI are shaded.
        
        Accepts:
            filename    the filename (including path and extension) for the
                        new plot to be created
            Xarr        an array of arrays of x values
            Yarr        an array of arrays consisting of a tuple of three 
                        Y values:
                            (5%, median, 95%)
            titles      an array of plot titles
            title       main plot title
            xlabel      x axis label string
            ylabel      y axis label string
            
        Doesn't return anything
        Writes a PDF file to filename.
        
        Examples for setting up multiple charts on shared axes are at:
        http://matplotlib.org/examples/pylab_examples/subplots_demo.html
    '''
    
    def adjustFigAspect(fig, aspect=1):
        '''
            Adjusts the whitespace around a figure so that each subplot 
            achieves the desired aspect ratio (square by default).
            Accepts a matplotlib figure object.
            Doesn't need to return anything because it directly modifies the 
            figure object.
        '''
        xsize, ysize = fig.get_size_inches()
        minsize = min(xsize, ysize)
        xlim = .4*minsize/xsize
        ylim = .4*minsize/ysize
        if aspect < 1:
            xlim *= aspect
        else:
            ylim /= aspect
        fig.subplots_adjust(left=.5-xlim,
                            right=.5+xlim,
                            bottom=.5-ylim,
                            top=.5+ylim)

    plt.clf()
    fig, axarr = plt.subplots(nrows, ncols, sharex=True, sharey=True)
    fig.suptitle(title, fontsize=titlesize)
    for ax, X, Y, title in zip(axarr.flat, Xarr, Yarr, titles):
        ax = contour_plot(ax, X, Y,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          titlesize=titlesize,
                          labelsize=labelsize,
                          ticklabelsize=ticklabelsize,
                          scaling=1./len(axarr.flat))
    adjustFigAspect(fig)
    plt.savefig(filename, transparent=True)


def write_contour_plot(filename, X, Y, title=None, xlabel=None, ylabel=None,
                       titlesize=24, labelsize=20, ticklabelsize=14):
    '''
        Produces a type of line chart, where for each x, the median value is 
        shown as a line, and the area between the 5% and 95% CI are shaded.
        
        Accepts:
            filename    the filename (including path and extension) for the
                        new plot to be created
            X           an array of x values
            Yt              an array of y tuples (95%, median, 5%)
            xlabel          x axis label string
            ylabel          y axis label string
            title           axis title
            titlesize       font size (in points) for titles
            labelsize       font size (in points) for x and y labels
            ticklabelsize   font size (in points) for xtick and ytick labels
            scaling         a float which scales the graph, including line widths
                            and all font sizes.
            
        Doesn't return anything
        Writes a PDF file to filename.
    '''
    plt.clf()
    fig = plt.figure()
    if title is not None:
        plt.title(title, fontsize=titlesize)
    ax = fig.add_subplot(111)
    ax = contour_plot(ax, X, Y, xlabel=xlabel, ylabel=ylabel,
                      labelsize=labelsize, ticklabelsize=ticklabelsize)
    plt.tight_layout()     # eliminates whitespace around the plot
    plt.savefig(filename, transparent=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path',
                        help = 'results folder path. If it does not exist, '\
                               'it will be created.',
                        nargs = '?',  # makes this argument optional
                        default = os.path.dirname(__file__))
    parser.add_argument('-o','--overwrite',action='store_true',
                        help = 'overwrite old csv files.')
    parser.add_argument('-g','--graph_only',action='store_true',
                        help = 'skip simulations and only generate graphs.')
    parser.add_argument('-v','--verbose',action='store_true',
                        help = 'also outputs the sample runs.')
        
    args=parser.parse_args()
        
    # Clean up the path name. Create a new directory if it does not yet exist.
    # Copy the source code into the directory, if it is not already there.
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
        print "Created folder '{}'.".format(args.path)
    else:
        print "Using folder '{}'.".format(args.path)
    source_fn = os.path.split(__file__)[-1].replace('.pyc','.py')
    shutil.copyfile(source_fn, os.path.join(args.path,source_fn))
    print "Copied source code to '{}'.".format(args.path)

    if not args.graph_only:
        for experiment in EXPERIMENTS:
            # This quick sample run obtains the headers for the data file.
            sample_run = varyAssort(**experiment)
            sample_run['header'][0] = '# ' + sample_run['header'][0]
            headers =  [['# experiment date = {timestamp}' \
                         ''.format(timestamp=time.strftime('%Y %b %d'))],
                        ['# source code = ' + source_fn],
                        ['# constant_pop_size = {constant_pop_size}'\
                         ''.format(**experiment)],
                        ['# aa_fitness = {aa_fitness}'\
                         ''.format(**experiment)],
                        ['# aa_homogamy = {aa_homogamy}'\
                         ''.format(**experiment)],
                        ['# mutation allele start freq = {A_FREQ}'\
                         ''.format(A_FREQ=A_FREQ)],
                        sample_run['header']]
            fn = os.path.join(args.path, 
                              'pop{constant_pop_size}_fitness{aa_fitness}'\
                              '_homogamy{aa_homogamy:.2}.tsv'\
                              ''.format(**experiment))
            if os.path.isfile(fn) and not args.overwrite:
                print "File '{}' exists.".format(fn)
                print "  Use --overwrite to re-do the experiment."
                continue
            else:
                f = open(fn,'wb')
                o = csv.writer(f, dialect=csv.excel_tab)
                o.writerows(headers)
                f.close()
                print "Created '{}'.".format(fn)
                 
            def worker():
                '''
                    The worker function exists as a convenient way of passing
                    varyAssort with its parameters to the multiprocessing pool.
                '''
                return varyAssort(**experiment)
                
            proposals = 0
            mp_chunk_size = cpu_count = multiprocessing.cpu_count()
            pool = multiprocessing.Pool()
            while proposals < PROPOSALS:
                start_time = time.time()
                p = [pool.apply_async(worker) for i in range(mp_chunk_size)]
                table = [item.get()['row'] for item in p]
                f = open(fn,'ab')
                o = csv.writer(f, dialect=csv.excel_tab)
                o.writerows(table)
                f.close()
                proposals += mp_chunk_size
                rate = mp_chunk_size*60./(time.time()-start_time)
                print '{proposals:,} proposals completed ' \
                      '({cpu_count} CPUs; {rate:,.0f} proposals/min)'\
                      ''.format(proposals=proposals,
                                rate=rate,
                                cpu_count=cpu_count)
                # mp_chunk_size is dynamically adjusted based on actual
                # execution speed such that file writes occur once per minute.
                mp_chunk_size = int(rate - rate%cpu_count)
                if proposals + mp_chunk_size > PROPOSALS:
                    mp_chunk_size = PROPOSALS-proposals

    # This graphing routine retrieves data from all the saved tsv files in
    # the directory (which would have been generated by simulation runs).
    # This allows tweaking of graphing routines without re-running simulations. 
    print "Looking for .tsv files in '{}'.".format(args.path)
    files = [f for f in os.listdir(args.path) \
                if os.path.isfile(os.path.join(args.path, f)) and '.tsv' in f]
    a_freqs_arr = []
    params_arr= []
    titles = []
    for file in files:
        f = open(os.path.join(args.path, file),'r')
        rows = csv.reader(f, dialect=csv.excel_tab)
        params = {}
        headers = []
        data = []
        for row in rows:
            if '#' in row[0] and '=' in row[0]:
                key = row[0].replace('#','').split('=')[0].strip()
                value = row[0].replace('#','').split('=')[1].strip()
                params[key] = value
                continue
            elif '#' in row[0]:
                headers += row
                continue
            else:
                data += [row]
        # This clever Python shorthand transposes a table. It will cause data
        # truncation if the table does not have consistent dimensions.
        data = zip(*data)
        
        params_arr.append(params)
        # Select and type convert X, which should be the same for all data.
        X = []
        for h, col in zip(headers, data):
            if 'gen' in h:
                X.append(int(col[0]))
        # Select a_freq.
        a_freqs = []
        for h, col in zip(headers, data):
            if 'a_freq' in h:
                col = map(float, col)
                col.sort()
                a_freqs.append((col[int(0.975*len(col))], 
                                numpy.median(col), 
                                col[int(0.025*len(col))]))
        a_freqs_arr.append(a_freqs)
        
        title='popsize={}, '\
              'aa_fitness={}, '\
              'aa_homogamy={}'\
              ''.format(params['constant_pop_size'],
                        params['aa_fitness'],
                        params['aa_homogamy'])
        titles.append(title)
        
        # Create individual contour charts
        filename = os.path.join(args.path,file.replace('.tsv','.pdf'))
        print "Saving chart to '{}'.".format(filename)
        write_contour_plot(filename, X, a_freqs,
                           title=title,
                           xlabel='Generation',
                           ylabel='Recessive allele frequency')
    filename = os.path.join(args.path, 'summary.pdf')
    print "Saving summary chart to '{}'.".format(filename)
    write_summary_contour_plot(filename, 
                               Xarr=[X for i in range(len(a_freqs))], 
                               Yarr=a_freqs_arr,
                               nrows=2,
                               ncols=2,
                               title='homogamy.py summary data',
                               titles=titles,
                               xlabel='Generation',
                               ylabel='Recessive allele frequency')