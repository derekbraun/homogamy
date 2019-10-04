#!/usr/local/bin/python3 -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
    Samir Jain, Eric Epstein, Trevor Klemp, Maggie Gray, Selman Jawed, Derek
    Braun* (*derek.braun@gallaudet.edu)

    Produces graphs for the data files created by simulator.py.
    Last updated: 23-Jun-2017 by Derek Braun
'''

AXIS_LABELS = {'a'  : 'Allelic frequency',
               'aa' : 'Deaf individuals',
               'aa_homogamy' : 'Homogamy',
               'aa_fitness' : 'Fitness',
               'F' : 'F'}

import os
import sys
import argparse
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, lines
import fileio

#see: http://matplotlib.sourceforge.net/users/customizing.html
BLUE = '#00457c'
BUFF = '#e8d4a2'


def violin_plot(filename, violins, title=None, xlabel=None, categories=None,
                ylabel=None, ylim=None, y_format='{x:.3%}', yscale='linear',
                rc=None, rcfname=None):
    '''
        Produces a violin plot.

        Accepts:
            filename        the file to be written
            violins         a multi-dimensional array of Y values
            title           the plot title
            xlabel          x axis label string
            categories      x axis categories (strings)
            ylabel          y axis label string
            ylim            y limit
            y_format        y axis formatter string
            yscale          linear or log scale
            rc              a rcparams dict. This takes precedence over
                            rcfname.
            rcfname         a rcparams file to load to set graph style.
    '''
    with matplotlib.rc_context(rc=rc, fname=rcfname):
        plt.clf()
        fig = plt.figure()
        #if title is not None:
        #    plt.title(title)
        ax = fig.add_subplot(111)
        ax.set_yscale(yscale)
        ax.set_ylim(0.001 if yscale == 'log' else 0,
                    numpy.amax(violins)*1.1 if ylim is None else ylim)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            if yscale == 'log':
                ylabel = 'log {}'.format(ylabel)
            if len(ylabel) < 4:
                ax.set_ylabel(ylabel, rotation=0)
            else:
                ax.set_ylabel(ylabel)
        ax.violinplot(violins, showmeans = False, showmedians = True,
                      showextrema = False)
        Y_98, Y_75, Y_50, Y_25, Y_2 = numpy.percentile(violins, [98, 75, 50, 25, 2], axis=1)
        inds = numpy.arange(1, len(Y_50) + 1)
        ax.vlines(inds, Y_25, Y_75, linestyle='-', lw=5)
        ax.vlines(inds, Y_2, Y_98, linestyle='-', lw=1)
        ax.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter(y_format))

        if categories is not None:
            ax.get_xaxis().set_tick_params(direction='out')
            ax.xaxis.set_ticks_position('bottom')
            ax.set_xticks(numpy.arange(1, len(categories) + 1))
            ax.set_xticklabels(categories)
            ax.set_xlim(0.25, len(categories) + 0.75)
        plt.savefig(filename, transparent=True)
        plt.close()


def contour_plot(filename, X, Ya, title=None, xlabel=None, ylabel=None,
                 ylim=None, y_format='{x:.3%}', rc=None, rcfname=None):
    '''
        Produces a contour plot that is like a continuous boxplot.
        The 25% and 75% quartiles, median, 2% and 98% credible intervals are all
        represented with lines. There is shading between the 25% and 75%
        quartiles just like with a boxplot.

        Accepts:
            filename        the file to be written
            X               an array of x values
            Ya              a multi-dimensional array of Y values
            title           the plot title
            xlabel          x axis label string
            ylabel          y axis label string
            ylim            y limit
            y_format        y axis formatter string
            rc              a rcparams dict. This takes precedence over
                            rcfname.
            rcfname         a rcparams file to load to set graph style.
    '''

    Y_98, Y_75, Y_50, Y_25, Y_2 = numpy.percentile(Ya, [98, 75, 50, 25, 2], axis=1)

    with matplotlib.rc_context(rc=rc, fname=rcfname):
        plt.clf()
        fig = plt.figure()
        grid = plt.GridSpec(1,10, wspace=0)
        ax1 = fig.add_subplot(grid[0,:9])
        ax2 = fig.add_subplot(grid[0,9], sharey=ax1)
        ax1.locator_params(axis='x', nbins=4)
        #if title is not None:
        #    plt.title(title, ha='right')
        ax1.set_xlim(0,max(X))                            # necessary to have the graph start at 0,0
        ax1.set_ylim(0,numpy.amax(Ya) if ylim is None else ylim)  # necessary to have the graph start at 0,0
        if xlabel is not None:
            ax1.set_xlabel(xlabel)
        if ylabel is not None:
            if len(ylabel) < 4:
                ax1.set_ylabel(ylabel, rotation=0)
            else:
                ax1.set_ylabel(ylabel)
        #shade in between the 25% and 75% quartiles, just like a boxplot
        ax1.fill_between(X, Y_75, Y_25, color=BUFF, alpha=0.5)
        ax1.plot(X, Y_98, color=BLUE, lw=0.5)
        ax1.plot(X, Y_75, color=BLUE)
        ax1.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter(y_format))
        ax1.plot(X, Y_50, color=BLUE, lw=1.5)
        ax1.plot(X, Y_25, color=BLUE)
        ax1.plot(X, Y_2, color=BLUE, lw=0.5)

        ax2.violinplot(Ya[-1], showmeans = False, showmedians = True,
                      showextrema = False)
        ax2.vlines([1], Y_25[-1], Y_75[-1], linestyle='-', lw=4)
        ax2.vlines([1], Y_2[-1], Y_98[-1], linestyle='-', lw=1)
        ax2.text(1.2, Y_2[-1], y_format.format(x=Y_2[-1]),
                        va='center', ha='left')
        ax2.text(1.3, Y_50[-1], y_format.format(x=Y_50[-1]),
                va='center', ha='left')
        ax2.text(1.2, Y_98[-1], y_format.format(x=Y_98[-1]),
                va='center', ha='left')
        ax2.axis('off')
        plt.savefig(filename, transparent=True)
        plt.close()

#
#   MAIN ROUTINE
#
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filenames',
                        nargs = '+',
                        help = 'filename(s) for data file')
    parser.add_argument('-o', '--output_file',
                        action='store',
                        default = None,
                        help = 'the output filename ')
    parser.add_argument('-r', '--rcfname',
                        action='store',
                        default = 'print.rc',
                        help = 'the rcf file to use. Must exist in directory')
    parser.add_argument('-f', '--field',
                        action='store',
                        help = 'the variable to compare among populations; ' \
                               ' e.g. a, aa, or F.')
    parser.add_argument('-t','--title',
                        action='store',
                        default = None,
                        help = 'user-specified title for plot')
    parser.add_argument('--ylim',
                        action='store',
                        type = float,
                        default = None,
                        help = 'manual y limit')
    parser.add_argument('--y_format',
                        action='store',
                        default = '{x:.3%}',
                        help = 'manual y axis formatter string')
    parser.add_argument('--yscale',
                        action='store',
                        default = 'linear',
                        help = 'manual y axis formatter string')
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

    # if only one file, do a contour plot
    if len(args.filenames) == 1:
        print('Writing contour plot using {}'.format(args.rcfname))
        e = experiments[0]
        default_title='pop size={:,}   fitness={:.1f}   homogamy={:.1f}'\
              ''.format(int(e.constant_pop_size),
                        float(e.aa_fitness),
                        float(e.aa_homogamy))
        rc = {'axes.titlesize': 10}
        X = e.select('gen',0)
        if args.output_file:
            filename = args.output_file
        else:
            filename = os.path.splitext(args.filenames[0])[0] + \
                   '.contour_plot.{args.field}.{args.rcfname}.pdf'.format(**locals())
        Ya = e.select(args.field)
        contour_plot(filename, X, Ya,
                     title=args.title if args.title is not None else default_title,
                     xlabel='Generation',
                     ylabel=AXIS_LABELS[args.field],
                     ylim=args.ylim,
                     y_format=args.y_format,
                     rc=rc,
                     rcfname=args.rcfname)
        print('   {}'.format(filename))
        print("Done.\n")

    # if more than one file, do a violin plot
    else:
        print('Identifying independent variable(s)...')
        indep_vars = []
        for var in fileio.INDEP_VARS:
            init_value = getattr(experiments[0],var)
            for e in experiments:
                if init_value != getattr(e, var):
                    indep_vars.append(var)
                    break
        if len(indep_vars) == 0:
            print('No independent variable identified. Exiting...')
            exit()
        elif len(indep_vars) > 1:
            print('   {}'.format(indep_vars))
            print('Code cannot handle this. Exiting...')
            exit()
        else:
            print('   {}'.format(indep_vars[0]))
            categories = []
            for e in experiments:
                categories.append(getattr(e, indep_vars[0]))
            print('Writing violin plot using {}'.format(args.rcfname))
            default_title=''
            if args.output_file:
                filename = args.output_file
            else:
                filename = os.path.splitext(args.filenames[0])[0] + \
                    '.violin_plot.{args.field}.{args.rcfname}.pdf'.format(**locals())
            violins = []
            for e in experiments:
                violins.append(e.select_endpoint(args.field))
            violin_plot (filename, violins,
                         title=args.title if args.title is not None else default_title,
                         xlabel=AXIS_LABELS[indep_vars[0]],
                         categories=categories,
                         ylabel=AXIS_LABELS[args.field],
                         ylim=args.ylim,
                         y_format=args.y_format,
                         yscale=args.yscale,
                         rcfname=args.rcfname)
            print('   {}'.format(filename))
        print("Done.\n")
