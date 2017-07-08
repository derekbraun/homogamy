#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
    Samir Jain, Eric Epstein, Trevor Klemp, Maggie Gray, Selman Jawed, Derek 
    Braun* (*derek.braun@gallaudet.edu)
    
    Produces graphs for the data files created by simulator.py.
    Last updated: 23-Jun-2017 by Derek Braun
'''


import os
import argparse
import numpy
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt, lines
import fileio

#see: http://matplotlib.sourceforge.net/users/customizing.html
BLUE = '#00457c'
BUFF = '#e8d4a2'


def hist_with_hpd(filename, dist, title=None, xlabel=None, hpd=0.95, 
                  num_of_bins=41, crop=0.998, format='{:,.0f}', rc=None):
    '''
        Plots a histogram with the two-tailed HPD region (default 95%) in a 
        different color. This coloring sets this routine apart from 
        pyplot.hist()
        
        Accepts:
            filename    a filename for the histogram to be written
            dist        a list of values
            title       an optional title
            xlabel      an optional xlabel
            
            hpd         optional hpd (default is 95%)
            num_of_bins number of histogram bins, normally set to 41
            crop        where you want the histogram to stop. best to leave alone
            rc              rcparams file to load to set graph style.
            
        Saves a .pdf plot to filename. Doesn't return anything.
    '''
    # test num_of_bins and raise an error if it won't work    
    if num_of_bins == 'auto':
        num_of_bins = min([len(dist)/100,100])
    if num_of_bins % 2 == 0:                    # force an odd number,
        num_of_bins += 1                        # because numpy.histogram 
    if len(dist) < num_of_bins:                 # is off with even numbers of bins
        raise('Cannot produce histogram. There must be more data than bins.')
    hist, bin_edges = numpy.histogram(dist, bins=num_of_bins)
    # find exact HPD values
    dist.sort()
    hpd_min_index = int((1-hpd)/2.*len(dist))
    hpd_min = dist[hpd_min_index]
    mean = numpy.mean(dist)
    for i, datum in enumerate(dist):
        if datum >= mean:
            mean_index = i
            break
    hpd_max_index = int((1+hpd)/2.*len(dist))
    hpd_max = dist[hpd_max_index]
    crop_min_index = int((1-crop)/2.*len(dist))
    crop_max_index = int((1+crop)/2.*len(dist))
    crop_min = dist[crop_min_index]
    crop_max = dist[crop_max_index]
    
    # organize bins based on HPD (for coloring) and
    # find bin indexes for cropping (crop_min_bin_index and crop_max_bin_index) and HPD labeling 
    # (hpd_min_bin_index and hpd_max_bin_index)
    crop_min_bin_index = False
    crop_max_bin_index = False
    hpd_min_bin_index = False
    hpd_max_bin_index = False
    mean_bin_index = False
    sig_bars = []
    insig_bars = []
    n = 0
    for i, bin in enumerate(hist):
        n += bin
        if not crop_min_bin_index and n > crop_min_index:
            crop_min_bin_index = i
        if hpd_min_index < n < hpd_max_index:
            if not hpd_min_bin_index:
                hpd_min_bin_index = i
            hpd_max_bin_index = i+1
            sig_bars.append(bin)
            insig_bars.append(0)
        else:
            sig_bars.append(0)
            insig_bars.append(bin)
        if not mean_bin_index and n >= mean_index:
            mean_bin_index = i            
        if not crop_max_bin_index and n >= crop_max_index:
            crop_max_bin_index = i
            break
    with matplotlib.rc_context(fname=rc):
        plt.clf()
        plt.figure()
        plt.subplot(111)
        plt.bar(numpy.arange(len(sig_bars)), sig_bars, width=1,
                facecolor=BLUE)
        plt.bar(numpy.arange(len(insig_bars)), insig_bars, width=1, 
                facecolor=BUFF)
        if xlabel:
            plt.xlabel(xlabel)
    
        # The ticks are crop_min, median, and crop_max
        tick_pos = [crop_min_bin_index, mean_bin_index, crop_max_bin_index]
        tick_labels = [format.format(crop_min), 'mean = ' + format.format(mean), 
                       format.format(crop_max)]
        plt.xticks(tick_pos, tick_labels)
    
        # turn off axis box, turn off y-axis, crop axes, and set ticks
        plt.box(on=None)
        ax = plt.axes()
        ax.yaxis.set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        plt.xlim(crop_min_bin_index, crop_max_bin_index) 
    
        # print and label hpd lines
        plt.text(hpd_min_bin_index-1, max(hist), format.format(hpd_min),
                    ha='right')
        plt.axvline(hpd_min_bin_index, ls=':',lw=2)
        plt.text(hpd_max_bin_index+1, max(hist), format.format(hpd_max),
                    ha='left')
        plt.axvline(hpd_max_bin_index, ls=':',lw=2)
        if title:
            plt.title(title)
        plt.grid(False)
        plt.savefig(filename, transparent=True)
        plt.close()


def multiline_plot(filename, X, Ya, title=None, xlabel=None, ylabel=None, 
                   rc=None):
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
            rc              rcparams file to load to set graph style.
    '''
    with matplotlib.rc_context(fname=rc):
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
                 gradients=32, scaling=1.0, rc=None):
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
            gradients       the number of gradients. More gradients means
                            a smoother density plot at the cost of a larger 
                            vector file.
            scaling         scaling factor used for creating subplots,
                            where the text and lines need to be scaled down.
            rc              rcparams file to load to set graph style.
    '''
    #_set_plot_params(use, scaling)
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

    with matplotlib.rc_context(fname=rc):
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
        print "    {:,} simulations x {:,} generations".format(len(e.data[0]), int(e.gen))
    else:
        print '  File not found.'
    
    
    if args.contour_plot:
        title='N={:,} fitness={:.1f} homogamy={:.1f}'\
          ''.format(int(e.constant_pop_size),
                    float(e.aa_fitness),
                    float(e.aa_homogamy))
        X = e.select('gen',0)
        for rc in ['print.rc']:
            for var in ['a','aa','F']:
                filename = os.path.splitext(args.filename)[0] + \
                           '.contour_plot.{var}.for_{rc}.pdf'.format(**locals())
                Ya = e.select(var)
                contour_plot(filename, X, Ya,
                             title=title,
                             xlabel='Generations',
                             ylabel=var,
                             rc=rc)
                print "  Writing '{}'".format(filename) 
    if args.multiline_plot:
        title='N={:,} fitness={:.1f} homogamy={:.1f}'\
          ''.format(int(e.constant_pop_size),
                    float(e.aa_fitness),
                    float(e.aa_homogamy))
        X = e.select('gen',0)
        for rc in ['print.rc']:
            for var in ['a','aa','F']:
                filename = os.path.splitext(args.filename)[0] + \
                           '.multiline_plot.{var}.for_{rc}.pdf'.format(**locals())
                Ya = e.select(var)
                multiline_plot(filename, X, Ya,
                               title=title,
                               xlabel='Generations',
                               ylabel=var,
                               rc=rc)
                print "  Writing '{}'".format(filename)
    if args.histogram:
        title='N={:,} fitness={:.1f} homogamy={:.1f}'\
          ''.format(int(e.constant_pop_size),
                    float(e.aa_fitness),
                    float(e.aa_homogamy))
        for rc in ['print.rc']:
            for var in ['a','aa','F']:
                filename = os.path.splitext(args.filename)[0] + \
                           '.histogram.{var}.for_{rc}.pdf'.format(**locals())
                dist = e.select(var)[-1]
                hist_with_hpd(filename, dist,
                               title=title,
                               xlabel=var,
                               rc=rc)
                print "  Writing '{}'".format(filename)
    print "Done."