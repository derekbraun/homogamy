#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-

'''
    Samir Jain, Eric Epstein, Derek Braun*
    *derek.braun@gallaudet.edu
    Summer of '14
    
    See README.md for to-do list and instructions for installing
    dependencies.
    
    We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/
'''

PARAMS = ['constant_pop_size','aa_fitness','aa_homogamy','a_freq',
          'experiment_date', 'source_file', 'data_file']
GALLAUDET_BLUE = '#003b65'
GALLAUDET_BUFF = '#e5d19e'
PEN_COLOR = 'black'

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
# The matplotlib.use command specifies the backend and *must* appear
# before importing pyplot. For choices of backends, see:
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
matplotlib.use('Agg') 
from matplotlib import pyplot as plt


class Experiment:
    '''
        This class is a convenient container for experimental parameters
        and data. It also handles saving data to tsv files and reading
        data from files.
        
        It operates on a little bit of black magic: it doesn't need to know
        the names of its own variables. That is, variables containing 
        experimental parameters (metadata) are passed to __init__ as kwargs
        and become variables in the class scope.
    '''
    def __init__(self, filename=None, **kwargs):
        '''
            If a filename is defined, that file will be opened and the file
            data will populate this class.
            
            Otherwise, variables containing experimental parameters (metadata)
            should be passed to __init__ as kwargs. They will become variables
            in the class scope. These parameters will be saved to tsv files
            by the write methods and may also be retrieved from the class scope.
        '''
        
        if filename is None:
            for key in kwargs.keys():
                setattr(self, key, kwargs[key])
            self.experiment_date = time.strftime('%Y %b %d')
            self.source_file = os.path.split(__file__)[-1].replace('.pyc','.py')
        elif filename is not None and len(kwargs) == 0:
            self._read()
        elif filename is not None and len(kwargs) > 0:
            raise BaseException, 'You passed both a filename and experimental '\
                                 'parameters to Experiment.__init__.'\
                                 'This is confusing and not allowed.'

    def write_headers(self, overwrite=False):
        if os.path.isfile(self.data_file) and not overwrite:
            return False
        else:
            h = []
            for param in PARAMS:
                h += ['# {}={}'.format(param, getattr(self, param))]
                h += self.column_headers
                print h
            f = open(self.data_file,'wb')
            o = csv.writer(f, dialect=csv.excel_tab)
            o.writerows([h])
            f.close()
            return True
    
    def write_data(self, data):
        '''
            Appends a block of data to an existing TSV file. Can and should be 
            called repeatedly as more data become available.
            
            Accepts: data, a table (a list of rows of columns).
            Returns: True if successful; False if the data file doesn't exist.
        '''
        if os.path.isfile(self.data_file):
            f = open(self.data_file,'ab')
            o = csv.writer(f, dialect=csv.excel_tab)
            o.writerows(table)
            f.close()
            return True
        else:
            return False
            
    def _read(self, filename):
        '''
            Internal method to extract data and experimental parameters from a
            TSV file.
        '''
        f = open(filename,'r')
        rows = csv.reader(f, dialect=csv.excel_tab)
        self.data = []
        for row in rows:
            if '#' in row[0] and '=' in row[0]:
                key = row[0].replace('#','').split('=')[0].strip()
                value = row[0].replace('#','').split('=')[1].strip()
                
                if not isattr(self, key):
                    setattr(self, key, value)
                    continue
            elif '#' in row[0]:
                self.column_headers = row
                continue
            else:
                self.data += [row]
        # This clever Python shorthand transposes a table. It will cause data
        # truncation if the table does not have consistent dimensions.
        self.data = zip(*self.data)
        
    def select(self, param):
        '''
            Selects the first instance of a data column with name param.
            Returns a list
        '''
        l = []
        for h, col in zip(self.column_headers, self.data):
            if param in h:
                l.append(int(col[0]))
                break
        return numpy.array(l)
        
    def hpd(self, param):
        '''
            Selects all instances of a data column with name param.
            Returns a list of tuples (2.5%, 50%, and 97.5%) for each column.
        '''
        lt = []
        for h, col in zip(self.column_headers, self.data):
            if param in h:
                col = map(float, col)
                col.sort()
                lt.append((col[int(0.975*len(col))], 
                           numpy.median(col), 
                           col[int(0.025*len(col))]))
        return lt


def contour_plot(ax, X, Yt, title=None, xlabel=None, ylabel=None,
                 titlesize=20, labelsize=20, ticklabelsize=16, scaling=1):
    '''
        Produces a type of line chart, where for each x, the median value is 
        shown as a line, and the area between the 5% and 95% CI are shaded.
        
        Accepts:
            ax              a matplotlib.pyplot axis instance
            X               an array of x values
            Yt              an array of y tuples (95%, median, 5%)
            title           optional axis title
            xlabel          optional x axis label string
            ylabel          optional y axis label string
            titlesize       optional font size (in points) for titles
            labelsize       optional font size (in points) for x and y labels
            ticklabelsize   optional font size (in points) for xtick and ytick
                            labels
            scaling         optional float which scales the graph, including 
                            line widths and fonts.
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


if __name__ == '__main__':
   print 'This module is not intended to be executed directly.'