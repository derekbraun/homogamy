#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-

'''
    The Experiment class standardizes file io.
'''


PARAMS = ['experiment_date',
          'source_code',
          'constant_pop_size',
          'gen',
          'a_freq', 
          'aa_fitness',
          'aa_homogamy']

import os
import time
import csv


def create_folder(path):
    '''
        Creates a folder, as needed.
        
        Accepts
            path        a folder
            
        Returns
            True        if a folder was created
            False       if a folder was not created
    '''
    if not os.path.isdir(path):
        os.makedirs(path)
        return True
    else:
        return False

class Experiment:
    '''
        This class is a container for experimental parameters and data. 
        It has methods for reading and saving data to tsv files, as well
        as for quick statistics and recovering data columns.
        
        This code is written to be very flexible. The names of the metadata
        headers are stored in the global variable PARAMS.
    '''
    def __init__(self, filename=None, **kwargs):
        '''
            If a filename is defined, that file will be opened and the file
            data will populate this class. Metadata in the file will 
            become variables in the class scope.
            
            If a filename is not defined, then metadata should be passed
            to __init__ as kwargs. They will become variables
            in the class scope.
        '''
        
        if filename is None:
            for key in kwargs.keys():
                setattr(self, key, kwargs[key])
            self.experiment_date = time.strftime('%Y %b %d')
            self.headers = None
        elif filename is not None and len(kwargs) == 0:
            self._read()
        elif filename is not None and len(kwargs) > 0:
            raise BaseException, 'You cannot pass both a filename and metatata'\
                                 ' to Experiment.__init__.'

    def write_metadata(self, overwrite=False):
        '''
            Writes metadata and headers to self.filename.
            
            Accepts:
                overwrite:      True or False
                
            Returns True if the file is written.
        '''
        if os.path.isfile(self.filename) and not overwrite:
            return False
        else:
            h = []
            for param in PARAMS:
                h += [['# {} = {}'.format(param, getattr(self, param))]]
            h += [self.headers]
            f = open(self.filename,'wb')
            o = csv.writer(f, dialect=csv.excel_tab)
            o.writerows(h)
            f.close()
            return True
    
    def write(self, rows):
        '''
            Appends a block of data to an existing .tsv file. Should be 
            called repeatedly as more data become available.
            
            Accepts:    rows    a list of rows of columns
            
            Returns:    True    if successful
                        False   if self.filename file doesn't exist
        '''
        if os.path.isfile(self.filename):
            f = open(self.filename,'ab')
            o = csv.writer(f, dialect=csv.excel_tab)
            o.writerows(rows)
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
                self.headers = row
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
        for h, col in zip(self.headers, self.data):
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
        for h, col in zip(self.headers, self.data):
            if param in h:
                col = map(float, col)
                col.sort()
                lt.append((col[int(0.975*len(col))], 
                           numpy.median(col), 
                           col[int(0.025*len(col))]))
        return lt