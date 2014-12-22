#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-

'''

    TO DO LIST:
    
    1.  Improve plot: outside ticks, axes slightly offset, no square box,
        charcoal pen color
        
    1.  Change plotting (and data collection?) to AA_freq not AA_size
        (Need to go through simuPOP examples to make this work)
        
    2.  Also record (hopefully) the inbreeding coefficient, F, after each
        generation. We might need to make the simupop object into an iterable
        to accomplish this.
    
    4.  Generate a final table that gives final medians and HPDs for
        AA_freq, Aa_freq, aa_freq, A_freq, a_freq, F.
    
    3.  Create a stand-alone routine to produce the summary chart:
    
        a.  We want to plot one row with a_freq and a second row with aa_freq.
            (4 experiments across), with one page (or a new figure?) per
            popsize.
    
        b.  Shared y-axis. Shared x-axis label only. Each experiment clearly
            identified.
            
    6. Plot the inbreeding coefficient, F, calculated at each time point,
       on top of the a_freq figures (and I guess, the AA_freq figures as well).



    Samir Jain, Eric Epstein, Derek Braun*
    *derek.braun@gallaudet.edu
    Summer of '14
    
    We follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/
    
    Because this code is intended to be archived, this doesn't follow the 
    don't-repeat-yourself (DRY) philosophy of object-oriented programming.
    All of our code is contained in this one script, and is not intended to
    be modular.
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
                # small pop, aa have 2x fitness, random mating
               {'constant_pop_size' : 10000,      
                'aa_fitness'        : 2,
                'aa_homogamy'       : 0.},
                # small pop, aa have 2x fitness, assortative mating
               {'constant_pop_size' : 10000,      
                'aa_fitness'        : 2,
                'aa_homogamy'       : 0.9}]

ASPECT_RATIO = (16,9)
GALLAUDET_BLUE = '#003b65'
GALLAUDET_BUFF = '#e5d19e'
CHART_COLORS = [GALLAUDET_BLUE,GALLAUDET_BUFF,'LightSkyBlue',
                'LightPink','LightGreen','LightSalmon']
PEN_COLOR = 'black'


def varyAssort(constant_pop_size, aa_fitness, aa_homogamy):
    '''
        This is the actual simulation routine.
        It calls simuPOP and passes on the appropriate parameters.
        
        Accepts:
        constant_pop_size   population size, which remains constant throughout
        aa_fitness          _relative_ fitness of deaf individuals
        aa_homogamy         the percent of assortative mating between
                            deaf individuals
        Returns: a dict containing the results from the simulation.
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


def contour_plot(ax, X, Y, xlabel='', ylabel='', title='', scaling=1.0):
    '''
        Produces a form of line chart, where the median is shown as a line
        and the area between the 5% and 95% CI are shaded.
        Accepts:
            ax          a matplotlib.pyplot axis instance
            X           an array of float or int x values
            Y           an array of tuples of float y values (95%, median, 5%)
            xlabel      x axis label string
            ylabel      y axis label string
            title       axis title
            scaling     a float which scales the graph. The default scaling
                        (1.0) is appropriate for poster size.
    '''
    # unpack Y tuple
    upper_CIs = [y[0] for y in Y]
    medians = [y[1] for y in Y]
    lower_CIs = [y[2] for y in Y]
    ax.set_xlim(min(X),max(X))
    if title <> '':
        ax.set_title(title, fontsize=20*scaling)
    ax.set_xlabel(xlabel, fontsize=20*scaling)
    ax.set_ylabel(ylabel, fontsize=20*scaling)
    ax.fill_between(X, upper_CIs, lower_CIs, color=GALLAUDET_BUFF)
    ax.text(max(X)*0.98, max(upper_CIs), '{:.2f}'.format(upper_CIs[-1]), 
            va='bottom', ha='right', fontsize=16*scaling)
    ax.plot(X, medians, color=GALLAUDET_BLUE, lw=3*scaling)
    ax.text(max(X)*0.98, max(medians), '{:.2f}'.format(medians[-1]), 
            va='bottom', ha='right', fontsize=16*scaling)
    ax.text(max(X)*0.98, min(lower_CIs), '{:.2f}'.format(lower_CIs[-1]), 
            va='top', ha='right', fontsize=16*scaling)
    ax.grid(False)         # turns off gridlines
    return ax
    
    
def write_summary_contour_plot(filename, Xarr, Yarr, cols, titles=[], title='',
                  xlabel='', ylabel=''):
    '''
        Writes a 2D array of scatter plots to filename.
        from http://matplotlib.org/examples/pylab_examples/subplots_demo.html
        
        The median is shown as a line and the area between the 5% and 95% CI 
        are shaded.
        Accepts:
            filename    the filename (including path and extension) for the
                        new plot to be created
            Xarr        an array of arrays of x values
            Yarr        an array of arrays of a tuple of three Y values:
                            (5%, median, 95%)
            titles      an array of plot titles
            title       main plot title
            xlabel      x axis label string
            ylabel      y axis label string
    '''
    rows = 1            # we might add the ability for more rows later
    scaling = 0.5       # controls the size of text, which needs to be made
                        # smaller as more subplots are added
    
    # The link below gives examples of setting up multiple charts on shared 
    # axes:
    # 
    # http://matplotlib.org/mpl_examples/pylab_examples/subplots_demo.py
    
    def adjustFigAspect(fig,aspect=1):
        '''
            Adjust the subplot parameters so that the figure has the correct
            aspect ratio.
        '''
        xsize,ysize = fig.get_size_inches()
        minsize = min(xsize,ysize)
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
    plt.title(title, y=1.10, fontsize=20*scaling)
    fig = plt.figure(figsize=ASPECT_RATIO)
    fig, axarr = plt.subplots(rows, cols, sharex=True, sharey=True)
    for ax, X, Y, title in zip(axarr, Xarr, Yarr, titles):
        ax = contour_plot(ax, X, Y,
                          xlabel=xlabel,
                          ylabel=ylabel, 
                          title=title.replace(',','\n'),
                          scaling=scaling)
    plt.tight_layout()
    # hide ytick labels in all plots but the leftmost
    plt.setp([ax.get_yticklabels() for ax in axarr[0:]], visible=False) 
    plt.setp([ax.get_yticklines() for ax in axarr[0:]], visible=False)
    plt.setp([ax.get_xticklabels() for ax in axarr[0:]], fontsize=16*scaling) 
    plt.savefig(filename, transparent=True)


def write_contour_plot(filename, X, Y, title='', xlabel='', ylabel=''):
    '''
        Writes a scatter plot to filename.
        The median is shown as a line and the area between the 5% and 95% CI 
        are shaded.
        Accepts:
            filename    the filename (including path and extension) for the
                        new plot to be created
            X           an array of x values
            Y           an array consisting of a tuple of three Y values:
                            (5%, median, 95%)
            xlabel      x axis label string
            ylabel      y axis label string
    '''
    plt.clf()
    fig = plt.figure(figsize=ASPECT_RATIO)
    plt.title(title, y=1.10, fontsize=20)
    ax = fig.add_subplot(111)
    ax = contour_plot(ax, X, Y, xlabel=xlabel, ylabel=ylabel)
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
        print "Saving chart to '{filename}'.".format(filename=filename)
        write_contour_plot(filename, X, a_freqs,
                           title=title,
                           xlabel='Generation',
                           ylabel='Recessive allele frequency')
    
    write_summary_contour_plot(args.path+'summary.pdf', 
                               Xarr=[X for i in range(len(a_freqs))], 
                               Yarr=a_freqs_arr,
                               cols=4,
                               titles=titles,
                               xlabel='Generation',
                               ylabel='Recessive allele frequency')