#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
    Samir Jain, Eric Epstein, Trevor Klemp, Maggie Gray, Selman Jawed, Derek 
    Braun* (*derek.braun@gallaudet.edu)
    
    Produces graphs for the data files created by simulator.py.
    Last updated: 23-Jun-2017 by Derek Braun
'''

#   To do list:
#   2a. Add the histogram code
#   3.  Streamline the copious rcparams code somehow. Either use existing 
#       styles, or move these into a local rcparams file, or dynamically use
#       the matplotlib.rc_params_from_file method.

import os
import argparse
import numpy
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt, lines
import fileio

#see: http://matplotlib.sourceforge.net/users/customizing.html
ASPECT_RATIO = (16,9)
BLUE = '#00457c'
BUFF = '#e8d4a2'
PLOT_PARAMS = { 'print':
                { 'aspect'      : 1.0,
                  'colors'      : ['black', 'white', BLUE,
                                    BUFF, 'LightSteelBlue',
                               'LightGoldenRodYellow','LightSkyBlue',
                               'LightPink','LightGreen','LightSalmon'],
                  'font'        : 12,
                  'title'       : 14,
                  'axes'        : 12,
                  'ticklabel'   : 10,
                  'minsize'     : 6,
                  'linewidth'   : 1,
                  'minwidth'    : 0.5},
                'slides_light_bg':
                { 'aspect'      : 1.6,
                  'colors'      : ['black', 'white', BLUE,
                                   BUFF, 'LightSteelBlue',
                                   'LightGoldenRodYellow','LightSkyBlue',
                                   'LightPink','LightGreen','LightSalmon'],
                  'font'        : 20,
                  'title'       : 14,
                  'axes'        : 20,
                  'ticklabel'   : 20,
                  'minsize'     : 8,
                  'linewidth'   : 1,
                  'minwidth'    : 0.5},
                'slides_dark_bg':
                {   'aspect'    : 1.6,
                    'colors'    : ['white', 'black', BUFF,
                                   BLUE, 'LightSteelBlue',
                                   'LightGoldenRodYellow','LightSkyBlue',
                                   'LightPink','LightGreen','LightSalmon'],
                    'font'      : 20,
                    'title'     : 14,
                    'axes'      : 20,
                    'ticklabel' : 20,
                    'minsize'   : 8,
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


def _adjustFigAspect(fig, aspect=1.5):
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
        fig.subplots_adjust(left=.375-xlim,
                            right=.625+xlim,
                            bottom=.375-ylim,
                            top=.625+ylim)


#
#   Put the histogram plot code here.
#


def hist_with_hpd(filename, data, title=None, xlabel=None, use='print'):
    '''
        Produces and writes a histogram where where a two-tailed HPD is 
        color-coded. This coloring sets this routine apart from pyplot.hist().
        
        Accepts:
            filename        the file to be written
            data            a list of values to be plotted
            title           the plot title
            xlabel          x axis label string
            use             determines the font size and line width according
                            to the global PLOT_PARAMS dict.
    '''
    CROP_HPD = .998
    hpd = 0.95
    num_of_bins = min([len(data)/100,100])
    hist, bin_edges = numpy.histogram(data,bins=25,range =(-0.006, 0.006))
    
    # find exact HPD values
    data.sort()
    hpd_min_index = int((1-hpd)/2.*len(data))
    hpd_min = data[hpd_min_index]
    hpd_max_index = int((1+hpd)/2.*len(data))
    hpd_max = data[hpd_max_index]
    crop_min_index = int((1-CROP_HPD)/2.*len(data))
    crop_max_index = int((1+CROP_HPD)/2.*len(data))
    
    # organize bins based on HPD (for coloring) and
    # find bin indexes for cropping (xmin and xmax) and HPD labeling 
    # (hpd_min_bin_index and hpd_max_bin_index)
    xmin = False
    xmax = False
    hpd_min_bin_index = False
    hpd_max_bin_index = False
    sig_bars = []
    insig_bars = []
    n = 0
    for i, bin in enumerate(hist):
        n += bin
        if hpd_min_index < n < hpd_max_index:
            if not hpd_min_bin_index:
                hpd_min_bin_index = i
            hpd_max_bin_index = i+1
            sig_bars.append(bin)
            insig_bars.append(0)
        else:
            sig_bars.append(0)
            insig_bars.append(bin)
        if not xmax and n >= crop_max_index:
            xmax = len(hist)
    xmin = 0
    #changes label to percentage
    hpd_min = hpd_min * 100
    hpd_max = hpd_max * 100
     
    _set_plot_params(use)
    plt.clf()
    fig = plt.figure(figsize=ASPECT_RATIO)
    ax = fig.add_subplot(111)
    if title is not None:
        plt.title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    ax.bar(numpy.arange(len(sig_bars)),sig_bars,width=1,color=BLUE)
    ax.bar(numpy.arange(len(insig_bars)),insig_bars,width=1,color=BLUE)
    ax.bar(numpy.arange(len(insig_bars)),insig_bars,width=1,color=BUFF)
    
    # calculate and print ticks
    tick_pos = [i for i, bin_edge in enumerate(bin_edges) \
        if i % 10 == 0]
    tick_labels = ['{:,.4f}'.format(bin_edge) for i, bin_edge \
        in enumerate(bin_edges) if i % 10 == 0]
         
    # turn off axis box, turn off y-axis, crop axes, and set ticks
    #ax.box(on=None)
    ax.yaxis.set_visible(False)
    ax.set_xlim(xmin,xmax)
    #ax.setp(ax.get_xticklabels(), visible=False)
    # label hpd lines
    ax.text(hpd_min_bin_index-1,max(hist) + 2000,'{:,.4f}'.format(hpd_min),
                ha='right')
    ax.axvline(hpd_min_bin_index,ls=':',lw=2)
    ax.text(hpd_max_bin_index+1,max(hist) + 2000,'{:,.4f}'.format(hpd_max),
                ha='left')
    ax.axvline(hpd_max_bin_index,ls=':',lw=2)
    #label zero line
    zero_index = i / 2.0
    ax.axvline(zero_index,color='black',ls=':',lw=2)
    ax.text(zero_index+0.75,max(hist) + 2000,'0',
                ha='right')
    ax.grid(False)
    plt.savefig(filename, transparent=True)
    plt.close()


def multiline_plot(filename, X, Ya, title=None, xlabel=None, ylabel=None, 
                   use='print'):
    '''
        Produces and writes a plot where each proposal is represented by its 
        own line.
        
        This type of plot is generally deprecated because it creates
        enormously large vector files that may not be printable.
        
        Accepts:
            filename        the file to be written
            X               an array of x values
            Ya              an array of an array of y values
            title           the plot title
            xlabel          x axis label string
            ylabel          y axis label string
            use             determines the font size and line width according
                            to the global PLOT_PARAMS dict.
    '''
    _set_plot_params(use)
    plt.clf()
    fig = plt.figure()
    if title is not None:
        plt.title(title)
    ax = fig.add_subplot(111)
    ax.set_xlim(min(X),max(X))
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        if len(ylabel) < 4:
            ax.set_ylabel(ylabel, rotation=0)
        else:
            ax.set_ylabel(ylabel)
    ax.plot(X, Ya, color=BLUE, alpha=20./len(Ya[0]))  # adjust the numerator to 
                                                      # get the appropriate line
                                                      # darkness
    #plt.yticks(numpy.arange(0, 2.6, 0.5))
    #plt.xticks(numpy.arange(start_year, start_year + 20 * gen + 1, 80))        
    #_adjustFigAspect(fig)
    plt.savefig(filename, transparent=True)
    plt.close()   



def contour_plot(filename, X, Ya, title=None, xlabel=None, ylabel=None, 
                 use='print', gradients=32, scaling=1.0):
    '''
        Produces a contour plot.
        The median values and 95% credible intervals are also represented 
        with lines.
        
        Accepts:
            filename        the file to be written
            X               an array of x values
            Ya              an array of an array of y values
            title           the plot title
            xlabel          x axis label string
            ylabel          y axis label string
            use             determines the font size and line width according
                            to the global PLOT_PARAMS dict.
            gradients       the number of gradients. More gradients means
                            a smoother density plot at the cost of a larger 
                            vector file.
            scaling         scaling factor used for creating subplots,
                            where the text and lines need to be scaled down.
    '''
    _set_plot_params(use, scaling)
    # calculate gradient
    Y_upper_cis = []
    Y_medians = []
    Y_lower_cis = []
    Ygrads = []
    for Y in Ya:
        Y = list(Y)
        Y = map(float, Y)
        Y.sort()
        Y_upper_cis.append(Y[int(0.975*len(Y))])
        Y_medians.append(numpy.median(Y))
        Y_lower_cis.append(Y[int(0.025*len(Y))])
    Ygrads = []
    for i in range(gradients):
        Yugs=[]
        Ylgs=[]
        for Y in Ya:
            # horribly inefficient
            Y = list(Y)
            Y = map(float, Y)
            Y.sort()
            step = len(Y)/(2*gradients)
            uidx = len(Y) - i*step - 1
            lidx = i*step
            Yugs.append(Y[uidx])
            Ylgs.append(Y[lidx])
        Ygrads.append((Yugs,Ylgs))

    _set_plot_params(use)
    plt.clf()
    fig = plt.figure()
    if title is not None:
        plt.title(title)
    ax = fig.add_subplot(111)
    ax.set_xlim(min(X),max(X))
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        if len(ylabel) < 4:
            ax.set_ylabel(ylabel, rotation=0)
        else:
            ax.set_ylabel(ylabel)
    for i, Yg in enumerate(Ygrads):
        ax.fill_between(X, Yg[0], Yg[1], color=BUFF, alpha=1./(gradients/4.))
    ax.plot(X, Y_upper_cis, color=BLUE)
    ax.text(max(X)*1.02, Y_upper_cis[-1], '{:.2f}'.format(Y_upper_cis[-1]), 
            va='bottom', ha='left')
    ax.plot(X, Y_medians, color=BLUE)
    ax.text(max(X)*1.02, max(Y_medians), '{:.2f}'.format(Y_medians[-1]), 
            va='bottom', ha='left')
    ax.plot(X, Y_lower_cis, color=BLUE)
    ax.text(max(X)*1.02, Y_lower_cis[-1], '{:.2f}'.format(Y_lower_cis[-1]), 
            va='top', ha='left')
    ax.grid(False)  
    #_adjustFigAspect(fig)
    plt.savefig(filename, transparent=True)
    plt.close()


#
#   MAIN ROUTINE
#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename',
                        help = 'filename for data file')
    parser.add_argument('-c','--contour_plot',action='store_true', 
                        default=True,
                        help = 'create contour plots')
    parser.add_argument('-m','--multiline_plot',action='store_true', 
                        default=False,
                        help = 'create multiline plots')
    parser.add_argument('-i','--histogram',action='store_true',
                        default=False,
                        help = 'create histograms')
    args=parser.parse_args()
    
    if os.path.isfile(args.filename):
        e = fileio.Experiment(args.filename)
        print "  Reading '{}'".format(args.filename)
        print "    {} proposals x {} generations".format(len(e.data[0]), e.gen)
    else:
        print '  File not found.'
    
    
    if args.contour_plot:
        title='N={:,} fitness={:.1f} homogamy={:.1f}'\
          ''.format(int(e.constant_pop_size),
                    float(e.aa_fitness),
                    float(e.aa_homogamy))
        X = e.select('gen',0)
        for use in ['print']:
            for var in ['a','aa','F']:
                filename = os.path.splitext(args.filename)[0] + \
                           '.contour_plot.{var}.for_{use}.pdf'.format(**locals())
                Ya = e.select(var)
                contour_plot(filename, X, Ya,
                             title=title,
                             xlabel='Generations',
                             ylabel=var,
                             use=use)
                print "  Writing '{}'".format(filename) 
    if args.multiline_plot:
        title='N={:,} fitness={:.1f} homogamy={:.1f}'\
          ''.format(int(e.constant_pop_size),
                    float(e.aa_fitness),
                    float(e.aa_homogamy))
        X = e.select('gen',0)
        for use in ['print']:
            for var in ['a','aa','F']:
                filename = os.path.splitext(args.filename)[0] + \
                           '.multiline_plot.{var}.for_{use}.pdf'.format(**locals())
                Ya = e.select(var)
                multiline_plot(filename, X, Ya,
                               title=title,
                               xlabel='Generations',
                               ylabel=var,
                               use=use)
                print "  Writing '{}'".format(filename)
    if args.histogram:
        title='N={:,} fitness={:.1f} homogamy={:.1f}'\
          ''.format(int(e.constant_pop_size),
                    float(e.aa_fitness),
                    float(e.aa_homogamy))
        for use in ['print']:
            for var in ['a','aa','F']:
                filename = os.path.splitext(args.filename)[0] + \
                           '.histogram.{var}.for_{use}.pdf'.format(**locals())
                Ya = e.select(var)
                hist_with_hpd(filename, Ya[-1],
                               title=title,
                               xlabel=var,
                               use=use)
    print "Done."