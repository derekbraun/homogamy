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

#   Needed changes:
#   1.  We need to be able to run this simulation with a smaller allele
#       frequency without running out of males or females.
#
#       To do this, change the Mating Scheme to generate more than 1 offspring 
#       per mating, which is not realistic anyway. Refer to Arnos paper on
#       fitness for what this number of offspring should be.
#
#       Bo Peng suggested starting with:
#       matingScheme = HeteroMating([HomoMating(subPops=[(0,1)], weight=aa_homogamy),
#                                    HomoMating(weight=1-aa_homogamy)])
#
#       This is close to the predefined RandomMating scheme. The above 
#       suggestion did not work "out of the box" because additional
#       parameters need to be specified. From looking at the _reference manual_
#       This looked easy enough and would not take long, but I am currently 
#       in a time crunch for my presentation at NIDCD.  
#
#   2.  Change the plot into a gradient density plot. This can be done by
#       creating many plots on the same axis with alpha channels to achieve
#       a gradient effect (or alternatively, creating overlapping plots using
#       zorders). Write this routine for a variable number of gradients. I 
#       think a default of 100 gradients would be more than enough. Too
#       many gradients would make the vector file too large.

DEBUG_MODE = False

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
import simuOpt
if DEBUG_MODE:
    PROPOSALS = 100
else:
    PROPOSALS = 1000
    simuOpt.setOptions(optimized=True, numThreads=0)
from simuPOP import *
print

A_FREQ = 0.1304   # should be 0.01304 but this currently causes problems with 
                  # running out of females
GEN = 40
EXPERIMENTS = [ # small pop, equal fitness, random mating
                {'constant_pop_size': 10000,      
                'aa_fitness'        : 1,
                'aa_homogamy'       : 0.,
                'a_freq'            : A_FREQ,
                'gen'               : GEN},
                # small pop, equal fitness, assortative mating
               {'constant_pop_size' : 10000,     
                'aa_fitness'        : 1,
                'aa_homogamy'       : 0.9,
                'a_freq'            : A_FREQ,
                'gen'               : GEN},
                # small pop, aa homozygotes have 2x fitness, random mating
               {'constant_pop_size' : 10000,      
                'aa_fitness'        : 1.2,
                'aa_homogamy'       : 0.,
                'a_freq'            : A_FREQ,
                'gen'               : GEN},
                # small pop, aa homozygotes have 2x fitness, assortative mating
               {'constant_pop_size' : 10000,      
                'aa_fitness'        : 1.2,
                'aa_homogamy'       : 0.9,
                'a_freq'            : A_FREQ,
                'gen'               : GEN}]


def simuAssortativeMatingWithFitness(constant_pop_size, gen, a_freq, 
                                     aa_fitness, aa_homogamy):
    '''
        Accepts:
        constant_pop_size   population size, which remains constant throughout
        gen                 number of generations
        a_freq              starting frequency of the a allele
        aa_fitness          _relative_ fitness of deaf (aa) individuals
        aa_homogamy         the percent of assortative mating between
                            deaf individuals

        Returns a dict containing the results from the simulation:
        gen                 generation number
        AA/Aa_size          size of the AA/aa population
        aa_size             size of the aa population
        A_freq              frequency of the A allele
        a_freq              frequency of the a allele
    '''        
    setRNG(random.seed(getRNG().seed()))
    pop = Population(constant_pop_size, loci=[1], infoFields='fitness')
    pop.dvars().header = [] 
    pop.dvars().row = []
    pop.setVirtualSplitter(GenotypeSplitter(loci=[0], alleles=[[0,0,0,1],[1,1]]))
    # Creates two virtual subpopulations needed in order to define a mating 
    # scheme: 
    #   Note: allele definitions are _unphased_ so (0,1) and (1,0) are equivalent
    #   alleles=[0,0,0,1] are individuals with 00 or 01 (AA/Aa)
    #   alleles=[1,1] are individuals with 11 (aa)
    pop.evolve(
        initOps= [InitSex(),
                  # assigns individuals randomly to be male or female.
                  # This can result in slightly more males or females, 
                  # which can cause errors if the wrong mating scheme is 
                  # selected.
                  InitGenotype(freq=[1-a_freq, a_freq])
                  ],
        preOps = [MapSelector(loci=[0], fitness={(0,0):1,
                                                 (0,1):1,
                                                 (1,1):aa_fitness})
                  # Assigns fitness values to individuals with different
                  # genotypes. This is stored in a field called 'fitness'
                  # by default, which is then applied by the appropriate
                  # mating scheme, also by default.
                  # Fitness in this case is relative fitness and is applied 
                  # during the mating scheme. Some mating schemes do not 
                  # support fitness, so check the documentation!
                  # If fitness is used, _all_ individuals must be assigned a 
                  # fitness or those with no assignment will be calculated as
                  # having zero fitness and will be discarded during the mating
                  # scheme.
                  ],
        matingScheme = HeteroMating([RandomMating(subPops=[(1,1)], weight=aa_homogamy),
                                     RandomMating(weight=1-aa_homogamy)]),
        
        postOps = [Stat(popSize=True, alleleFreq=[0], subPops=[(0,0), (0,1)]),
                   PyExec(r"header += ['gen','AA/Aa_size','aa_size','A_freq','a_freq']"),
                   PyExec(r"row += [gen,"\
                                   "subPopSize[0],"\
                                   "subPopSize[1],"\
                                   "alleleFreq[0][0], alleleFreq[0][1]]")
                   ],
        gen = gen
    )
        
    
    return {'header':pop.dvars().header,
            'row':pop.dvars().row}


#
#   GRAPHING
#   

import matplotlib
# The matplotlib.use command, which specifies the backend, *must* appear
# before importing pyplot. For choices of backends, see:
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
matplotlib.use('pdf') 
from matplotlib import pyplot as plt, lines

#see: http://matplotlib.sourceforge.net/users/customizing.html
ASPECT_RATIO = (16,9)
GALLAUDET_BLUE = '#003b65'
GALLAUDET_BUFF = '#e5d19e'
PLOT_PARAMS = { 'print': 
                {   'aspect'    : 1.0,
                    'colors'    : ['black', 'white', GALLAUDET_BLUE,
                                    GALLAUDET_BUFF, 'LightSteelBlue',
                                    'LightGoldenRodYellow','LightSkyBlue',
                                    'LightPink','LightGreen','LightSalmon'],
                    'font'      : 12,
                    'title'     : 14,
                    'axes'      : 12,
                    'ticklabel' : 10,
                    'minsize'   : 6,
                    'linewidth' : 1,
                    'minwidth'  : 0.5},
                'slides_light_bg': 
                {   'aspect'    : 1.6,
                    'colors'    : ['black', 'white', GALLAUDET_BLUE,
                                    GALLAUDET_BUFF, 'LightSteelBlue',
                                    'LightGoldenRodYellow','LightSkyBlue',
                                    'LightPink','LightGreen','LightSalmon'],
                    'font'      : 12,
                    'title'     : 14,
                    'axes'      : 12,
                    'ticklabel' : 10,
                    'minsize'   : 6,
                    'linewidth' : 1,
                    'minwidth'  : 0.5},
                'slides_dark_bg': 
                {   'aspect'    : 1.6,
                    'colors'    : ['white', 'black', GALLAUDET_BUFF,
                                    GALLAUDET_BLUE, 'LightSteelBlue',
                                    'LightGoldenRodYellow','LightSkyBlue',
                                    'LightPink','LightGreen','LightSalmon'],
                    'font'      : 12,
                    'title'     : 14,
                    'axes'      : 12,
                    'ticklabel' : 10,
                    'minsize'   : 6,
                    'linewidth' : 1,
                    'minwidth'  : 0.5}}
    
def _set_plot_params(use='print', scaling=1.0):
    '''
        This routine sets all the plot parameters globally, which takes effect
        with the next plot. Setting these parameters globally using 
        matplotlib's rcParams simplifies the downstream charting code because
        it is no longer necessary specify the size and color parameters
        in every charting command.
        
        In the future, consider: 
            (a) using matplotlib.rc_params_from_file, and setting up a 
                different file for each use.
            (b) using the same rcParams name for each dict field, then
                the rcparams can just be loaded iteratively, although
                this would complicate scaling. Scaling would need to be done
                by reading the appropriately scalable rcParams, scaling,
                then rewriting (which could also be done iteratively).
    '''
    d = PLOT_PARAMS[use]
    matplotlib.rc('font',  size=max(d['font']*scaling, d['minsize']))
    #matplotlib.rc('figure',titlesize=d['title']*scaling)
    matplotlib.rc('font',  size=max(d['font']*scaling, d['minsize']))
    matplotlib.rc('text',  color=d['colors'][0],
                           usetex=False)
    matplotlib.rc('axes',  titlesize=max(d['title']*scaling, d['minsize']),
                           labelsize=max(d['axes']*scaling, d['minsize']),
                           labelcolor=d['colors'][0],
                           linewidth=max(d['linewidth']*scaling, d['minwidth']),
                           edgecolor='#cccccc')                           
    matplotlib.rc('xtick', labelsize=max(d['ticklabel']*scaling, d['minsize']),
                           color=d['colors'][0])
    matplotlib.rc('ytick', labelsize=max(d['ticklabel']*scaling, d['minsize']),
                           color=d['colors'][0])
    matplotlib.rc('lines', linewidth = max(d['linewidth']*scaling, d['minwidth']),
                           color=d['colors'][0])
    matplotlib.rc('figure.subplot', wspace=0.3,
                                    hspace=0.3)
    return d['colors']

def contour_plot(ax, X, Yt, title=None, xlabel=None, ylabel=None, 
                            use='print', scaling=1.0):
    '''
        Produces a type of line chart, where for each x, the median value is 
        shown as a line, and the area between the 5% and 95% CI are shaded.
        
        Accepts:
            ax              a matplotlib.pyplot axis instance
            X               an array of x values
            Yt              an array of y tuples (95%, median, 5%)
            xlabel          x axis label string
            ylabel          y axis label string
            use         determines the font size and line width according
                            to the global PLOT_PARAMS dict.
    '''
    _set_plot_params(use, scaling)
    # unpack Y tuple
    Y_uppers = [t[0] for t in Yt]
    Y_medians = [t[1] for t in Yt]
    Y_lowers = [t[2] for t in Yt]
    ax.set_xlim(min(X),max(X))
    if title is not None:
        ax.set_title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.fill_between(X, Y_uppers, Y_lowers, color=GALLAUDET_BUFF)
    ax.text(max(X)*1.02, max(Y_uppers), '{:.2f}'.format(Y_uppers[-1]), 
            va='bottom', ha='left')
    ax.plot(X, Y_medians, color=GALLAUDET_BLUE, lw=3*scaling)
    ax.text(max(X)*1.02, max(Y_medians), '{:.2f}'.format(Y_medians[-1]), 
            va='bottom', ha='left')
    ax.text(max(X)*1.02, min(Y_lowers), '{:.2f}'.format(Y_lowers[-1]), 
            va='top', ha='left')
    
    ax.grid(False)      # turns off gridlines
    return ax
    

def write_summary_contour_plot(filename, Xarr, Yarr, nrows, ncols, titles=[], 
                               title=None, xlabel=None, ylabel=None,
                               use='print'):
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
    
    def adjustFigAspect(fig, aspect=1.5):
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

    _set_plot_params(use, scaling=1./nrows)
    plt.clf()
    fig, axarr = plt.subplots(nrows, ncols, sharex=True, sharey=True)
    fig.suptitle(title)
    for ax, X, Y, title in zip(axarr.flat, Xarr, Yarr, titles):
        ax = contour_plot(ax, X, Y,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          use=use,
                          scaling=1./nrows)
    adjustFigAspect(fig)
    plt.savefig(filename, transparent=True)
    plt.close()


def write_contour_plot(filename, X, Y, title=None, xlabel=None, ylabel=None,
                       use='print'):
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
            
        Doesn't return anything
        Writes a PDF file to filename.
    '''
    _set_plot_params(use)
    plt.clf()
    fig = plt.figure()
    if title is not None:
        plt.title(title)
    ax = fig.add_subplot(111)
    ax = contour_plot(ax, X, Y, 
                      xlabel=xlabel,
                      ylabel=ylabel,
                      use=use)
    plt.savefig(filename, transparent=True)
    plt.close()


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
            sample_run = simuAssortativeMatingWithFitness(**experiment)
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
                        ['# a_freq = {a_freq}'\
                         ''.format(**experiment)],
                        ['# gen = {gen}'\
                         ''.format(**experiment)],
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
                    simuAssortativeMatingWithFitness with its parameters to the multiprocessing pool.
                '''
                return simuAssortativeMatingWithFitness(**experiment)
                
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
        # note: the transposed table is a list of tuples, not a list of lists!
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
        title='N={:,}   '\
              'fitness={:.1f}    '\
              'homogamy={:.1f}'\
              ''.format(int(params['constant_pop_size']),
                        float(params['aa_fitness']),
                        float(params['aa_homogamy']))
        titles.append(title)
        
        # Create individual contour charts
        for use in ['print','slides_light_bg', 'slides_dark_bg']:
            new_ext = '.{}.pdf'.format(use)
            filename = os.path.join(args.path, file.replace('.tsv', new_ext))
            print "Saving chart to '{}'.".format(filename)
            write_contour_plot(filename, X, a_freqs,
                               title=title,
                               xlabel='Generations',
                               ylabel='Recessive Allele Frequency',
                               use=use)
    
    # Create summary contour charts
    for use in ['print','slides_light_bg', 'slides_dark_bg']:
        bn = 'summary.{}.pdf'.format(use)
        filename = os.path.join(args.path, bn)
        print "Saving summary chart to '{}'.".format(filename)
        write_summary_contour_plot(filename, 
                                   Xarr=[X for i in range(len(a_freqs))], 
                                   Yarr=a_freqs_arr,
                                   nrows=2,
                                   ncols=2,
                                   title='',
                                   titles=titles,
                                   xlabel='Generations',
                                   ylabel='Recessive Allele Frequency',
                                   use = use)