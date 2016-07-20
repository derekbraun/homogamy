#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-

'''
    See: README.md for To-Do List.
    Samir Jain, Eric Epstein, Derek Braun*
    *derek.braun@gallaudet.edu
    Summer of '14
    
    We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/
    
    Lite code intended only for debugging.
'''

import sys
import os
import time
import random
import argparse
import shutil
import csv
import simuOpt

    
import simuPOP as sim
import numpy as np
                                # Removal of Print
#A_FREQ = 0.01304                # Able to execute with low frequencies as of now.
A_FREQ = 0.5               

GEN = 70
EXPERIMENTS = [ # small pop, equal fitness, random mating
                {'constant_pop_size': 10000,      
                'aa_fitness'        : 1.0,
                'aa_homogamy'       : 0.0,
                'a_freq'            : A_FREQ,
                'gen'               : GEN},
                # small pop, equal fitness, assortative mating
               {'constant_pop_size' : 10000,     
                'aa_fitness'        : 1.0,
                'aa_homogamy'       : 0.9,
                'a_freq'            : A_FREQ,
                'gen'               : GEN},
                # small pop, aa homozygotes have 2x fitness, random mating
               {'constant_pop_size' : 10000,      
                'aa_fitness'        : 2.0,
                'aa_homogamy'       : 0.0,
                'a_freq'            : A_FREQ,
                'gen'               : GEN},
                # small pop, aa homozygotes have 2x fitness, assortative mating
               {'constant_pop_size' : 10000,      
                'aa_fitness'        : 2.0,
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
    sim.setRNG(random.seed(sim.getRNG().seed()))    
    pop = sim.Population(constant_pop_size, loci=[1], infoFields=['fitness'])
    pop.dvars().header = [] 
    pop.dvars().row = []
    pop.setVirtualSplitter(sim.GenotypeSplitter(loci=[0], alleles=[[0,0,0,1],[1,1]]))
    # Creates two virtual subpopulations needed in order to define a mating 
    # scheme: 
    #   Note: allele definitions are _unphased_ so (0,1) and (1,0) are equivalent
    #   alleles=[0,0,0,1] are individuals with 00 or 01 (AA/Aa)
    #   alleles=[1,1] are individuals with 11 (aa)
    pop.evolve(
        initOps= [sim.InitSex(),
                  # Assigns individuals randomly to be male or female.
                  # This can result in slightly more males or females, 
                  # which can cause errors if the wrong mating scheme is 
                  # selected.
                  sim.InitGenotype(freq=[1-a_freq, a_freq])
                  ],
        preOps = [sim.MapSelector(loci=[0], fitness={(0,0):1,
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
        matingScheme = sim.HeteroMating([
                            sim.HomoMating(chooser=sim.RandomParentsChooser(),
                                generator=sim.OffspringGenerator(
                                sim.MendelianGenoTransmitter()),
                                subPops=[(1,1)], weight=aa_homogamy),
                            sim.HomoMating(chooser=sim.RandomParentsChooser(),
                                generator=sim.OffspringGenerator(
                                sim.MendelianGenoTransmitter()),
                                weight=1-aa_homogamy)]),
                    
        
        postOps = [sim.Stat(popSize=True, alleleFreq=[0], subPops=[(0,0),(0,1)], 
                                genoFreq=[0], inbreeding=0), 
                  # Addition of genoFreq to establish a parameter for counting individuals
                  # with a specific genotype. 
                   sim.PyExec(r"header += ['gen','A', 'a',"\
                                   "'AA', 'Aa', 'aa',"\
                                   "'AA_size', 'Aa_size', 'aa_size',"\
                                   "'FAA','FAa','Faa']"),
                  # Addition of Aa_size and AA_size. AA/Aa_Size has a capitalization
                  # to prevent error-reading mistakes.
                   sim.PyExec(r"row += [gen, alleleFreq[0][0], alleleFreq[0][1],"\
                                   "genoFreq[0][(0,0)],"\
                                   "genoFreq[0][(0,1)]+genoFreq[0][(1,0)],"\
                                   "genoFreq[0][(1,1)],"\
                                   "genoNum[0][(0,0)],"\
                                   "genoNum[0][(0,1)]+genoNum[0][(1,0)],"\
                                   "genoNum[0][(1,1)],"\
                                   "(genoFreq[0][(0,0)]-alleleFreq[0][0]**2)/"\
                                   "(alleleFreq[0][0]-alleleFreq[0][0]**2),"\
                                   "1-(genoFreq[0][(0,1)]+genoFreq[0][(1,0)])/"\
                                   "(2*alleleFreq[0][0]*alleleFreq[0][1]),"\
                                   "(genoFreq[0][(1,1)]-alleleFreq[0][1]**2)/"\
                                   "(alleleFreq[0][1]-alleleFreq[0][1]**2)]")
                  # Addition of genoNum[x][(x,x)] to count the number of individuals
                  # with that specific genotype.
                  # You can add frequencies and sizes through addition. There are two
                  # sizes for both Aa carriers (0,1) and (1,0). Adding these can be done
                  # as seen above.      
                   ],
        gen = gen
    )
    return {'header':pop.dvars().header, 'row':pop.dvars().row}

           
if __name__ == '__main__':

    for experiment in EXPERIMENTS[1:2]:
        # This quick sample run obtains the headers for the data file.
        sample_run = simuAssortativeMatingWithFitness(**experiment)
        sample_run['header'][0] = '# ' + sample_run['header'][0]
        headers =  [['# experiment date = {timestamp}' \
                     ''.format(timestamp=time.strftime('%Y %b %d'))],
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
        print headers[0:6]
        print
        for h in sample_run['header'][0:12]:
            print "{h:>10}".format(h=h),
        print
        for h in sample_run['header'][0:12]:
            print " ---------",
        print
        for gen in range(0,len(sample_run['row'])/12):
            for datum in sample_run['row'][12*gen:12*(1+gen)]:
                if type(datum) is int or datum == int(datum):
                    print " {datum:>9}".format(datum=int(datum)),
                else:
                    print " {datum:>9.3f}".format(datum=datum),
            print
        print