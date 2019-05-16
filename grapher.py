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
                  num_of_bins=41, crop=0.998, format='{:,.0f}',
                  rc=None, rcfname=None):
    '''
        Plots a histogram with the two-tailed HPD region (default 95%) in a
        different color. This coloring sets this routine apart from
        pyplot.hist()

        Accepts:
            filename        a filename for the histogram to be written
            dist            a list of values
            title           an optional title
            xlabel          an optional xlabel

            hpd             optional hpd (default is 95%)
            num_of_bins     number of histogram bins, normally set to 41
            crop            where you want the histogram to stop. best to leave alone
            rc              a rcparams dict. This takes precedence over
                            rcfname.
            rcfname         a rcparams file to load to set graph style.

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
    hpd_min_index = int((1.-hpd)/2.*len(dist))
    hpd_min = dist[hpd_min_index]
    median = numpy.median(dist)
    for i, datum in enumerate(dist):
        if datum >= median:
            median_index = i
            break
    hpd_max_index = int((1.+hpd)/2.*len(dist))
    hpd_max = dist[hpd_max_index]
    crop_min_index = int((1.-crop)/2.*len(dist))
    crop_max_index = int((1.+crop)/2.*len(dist))
    crop_min = dist[crop_min_index]
    crop_max = dist[crop_max_index]

    # organize bins based on HPD (for coloring) and
    # find bin indexes for cropping (crop_min_bin_index and crop_max_bin_index) and HPD labeling
    # (hpd_min_bin_index and hpd_max_bin_index)
    crop_min_bin_index = False
    crop_max_bin_index = False
    hpd_min_bin_index = False
    hpd_max_bin_index = False
    median_bin_index = False
    sig_bars = []
    insig_bars = []
    n = 0
    for i, bin in enumerate(hist):
        n += bin
        if not crop_min_bin_index and n > crop_min_index:
            crop_min_bin_index = i
        if hpd_min_index <= n <= hpd_max_index:
            if hpd_min_bin_index is False:
                hpd_min_bin_index = i
            hpd_max_bin_index = i+1
            sig_bars.append(bin)
            insig_bars.append(0)
        else:
            sig_bars.append(0)
            insig_bars.append(bin)
        if median_bin_index is False and n >= median_index:
            median_bin_index = i
        if not crop_max_bin_index and n >= crop_max_index:
            crop_max_bin_index = i
            break
    with matplotlib.rc_context(rc=rc, fname=rcfname):
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
        tick_pos = [hpd_min_bin_index-0.5, median_bin_index, hpd_max_bin_index-0.5]
        tick_labels = ['{:.3f}'.format(hpd_min), '{:.3f}'.format(median),
                       '{:.3f}'.format(hpd_max)]
        plt.xticks(tick_pos, tick_labels)
        # print and label hpd lines
        plt.axvline(hpd_min_bin_index-0.5, ls=':')
        plt.axvline(median_bin_index, ls=':')
        plt.axvline(hpd_max_bin_index-0.5, ls=':')
        if title:
            plt.title(title)
        # format plot: turn off axis box, turn off y-axis, crop axes,
        # set tick visibility, set limits
        plt.box(on=None)
        ax = plt.axes()
        ax.yaxis.set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        #plt.xlim(crop_min_bin_index, crop_max_bin_index)
        plt.grid(False)
        plt.savefig(filename, transparent=True)
        plt.close()


def contour_plot(filename, X, Ya, title=None, xlabel=None, ylabel=None,
                 ylim=None, text_format='{:.3%}', rc=None, rcfname=None):
    '''
        Produces a contour plot that is like a continuous boxplot.
        The 25% and 75% quartiles, median, 2% and 98% credible intervals are all
        represented with lines. There is shading between the 25% and 75%
        quartiles just like with a boxplot.

        Accepts:
            filename        the file to be written
            X               an array of x values
            Ya              an array of an array of y values
            title           the plot title
            xlabel          x axis label string
            ylabel          y axis label string
            text_format     format string for label text
            ylim            y limit
            rc              a rcparams dict. This takes precedence over
                            rcfname.
            rcfname         a rcparams file to load to set graph style.
    '''
    # calculate gradient
    Y_98 = []
    Y_75 = []
    Y_medians = []
    Y_25 = []
    Y_2 = []
    for Y in Ya:
        Y = list(Y)
        Y = map(float, Y)
        Y.sort()
        Y_98.append(Y[int(0.98*len(Y))])
        Y_75.append(Y[int(0.75*len(Y))])
        Y_medians.append(numpy.median(Y))
        Y_25.append(Y[int(0.25*len(Y))])
        Y_2.append(Y[int(0.02*len(Y))])

    with matplotlib.rc_context(rc=rc, fname=rcfname):
        plt.clf()
        fig = plt.figure()
        plt.locator_params(axis='x', nbins=4)
        if title is not None:
            plt.title(title)
        ax = fig.add_subplot(111)
        ax.set_xlim(0,max(X))                            # necessary to have the graph start at 0,0
        ax.set_ylim(0,max(Y) if ylim is None else ylim)  # necessary to have the graph start at 0,0
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            if len(ylabel) < 4:
                ax.set_ylabel(ylabel, rotation=0)
            else:
                ax.set_ylabel(ylabel)
        #shade in between the 25% and 75% quartiles, just like a boxplot
        ax.fill_between(X, Y_75, Y_25, color=BUFF, alpha=0.5)
        ax.plot(X, Y_98, color=BLUE, lw=0.5)
        ax.plot(X, Y_75, color=BLUE)
        #ax.set_yticklabels([format.format(x) for x in ax.get_yticks().tolist()])
        #ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        ax.text(max(X)*1.02, Y_98[-1], text_format.format(Y_98[-1]),
                va='center', ha='left')
        ax.plot(X, Y_medians, color=BLUE, lw=1.5)
        ax.text(max(X)*1.02, Y_medians[-1], text_format.format(Y_medians[-1]),
                va='center', ha='left')
        ax.plot(X, Y_25, color=BLUE)
        ax.plot(X, Y_2, color=BLUE, lw=0.5)
        ax.text(max(X)*1.02, Y_2[-1], text_format.format(Y_2[-1]),
                va='center', ha='left')
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
    parser.add_argument('-c','--contour_plot',
                        action='store_true',
                        default=False,
                        help = 'create contour plot(s)')
    parser.add_argument('-y','--ylim',
                        action='store',
                        type=float,
                        default=None,
                        help = 'manual y limit for contour plot')
    parser.add_argument('-t','--text_format',
                        action='store',
                        default='{:.3%}',
                        help = 'text format for labels')
    parser.add_argument('-i','--histogram',
                        action='store_true',
                        default=False,
                        help = 'create histogram(s)')
    args=parser.parse_args()

    if os.path.isfile(args.filename):
        e = fileio.Experiment(args.filename)
        print "  Reading '{}'".format(args.filename)
        print "    {:,} simulations x {:,} generations".format(len(e.data[0]), int(e.generations))
    else:
        print '  File not found.'

    title='pop size={:,}    fitness={:.1f}    homogamy={:.1f}    n={:,} simulations'\
          ''.format(int(e.constant_pop_size),
                    float(e.aa_fitness),
                    float(e.aa_homogamy),
                    len(e.data[0]))
    rc = {'axes.titlesize': 10}

    if args.contour_plot:
        X = e.select('gen',0)
        for rcfname in ['print.rc']:
            for var in ['a','aa','F']:
                filename = os.path.splitext(args.filename)[0] + \
                           '.contour_plot.{var}.{rcfname}.pdf'.format(**locals())
                Ya = e.select(var)
                contour_plot(filename, X, Ya,
                             title=title,
                             xlabel='Generations',
                             ylabel=var,
                             ylim=args.ylim,
                             text_format=args.text_format,
                             rc=rc,
                             rcfname=rcfname)
                print "  Writing '{}'".format(filename)
                print "Done."
    elif args.histogram:
        for rcfname in ['print.rc']:
            for var in ['a','aa','F']:
                filename = os.path.splitext(args.filename)[0] + \
                           '.histogram.{var}.{rcfname}.pdf'.format(**locals())
                dist = e.select(var)[-1]
                hist_with_hpd(filename, dist,
                               title=title,
                               xlabel=var,
                               text_format=args.text_format,
                               rc=rc,
                               rcfname=rcfname)
                print "  Writing '{}'".format(filename)
                print "Done."
    else:
        print "Select either -c for contour plot or -h for histogram."
