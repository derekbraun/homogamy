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
    
    SimuPOP resources:
        Documentation: http://simupop.sourceforge.net/manual_svn/index.html
        User's guide:  http://simupop.sourceforge.net/manual_svn/build/userGuide.html
'''

# INSTALLATION:
# 
# 1.  Install this code from github.
# 2.  Install the major dependencies. Directions for OSX are on my blog:
#     https://derekbraun.wordpress.com/
# 3.  Install simuPOP
# 
# 
# TO-DO LIST
#         
# 1.  Define an Experiment class that manages experimental data and metadata.
#     Give the class methods for reading and writing data.
#     This should make code easier to follow.
#     
# 2.  simuPOP obviously does not work. Bo Peng suggested starting with:
#       matingScheme = HeteroMating([HomoMating(subPops=[(0,1)], weight=aa_homogamy),
#                                    HomoMating(weight=1-aa_homogamy)])
#
#       This is close to the predefined RandomMating scheme. The above 
#       suggestion did not work "out of the box" because additional
#       parameters need to be specified. From looking at the _reference manual_
#       This looked easy enough and would not take long, but I am currently 
#       in a time crunch for my presentation at NIDCD. 
# 
# 3.  Change plotting (and data collection) to show the freq of AA, Aa, aa,
#     A, and a; and calculate F. Plot F on one of the tables.
# 
# 5.  Generate and output a final table that gives final medians and HPDs
#     for AA, Aa, aa, A and a.


DEBUG_MODE = True

import sys
import os                      # Removed duplicate (import time).
import time
import random
import argparse
import shutil
import csv
import multiprocessing
import simuOpt

if DEBUG_MODE:
    PROPOSALS = 1
else:
    PROPOSALS = 10000
    simuOpt.setOptions(optimized=True, numThreads=0, quiet=True)
    
import GraphingMode as gm       # Additional imports for code stability.
import simuPOP as sim
import numpy as np
A_FREQ = 0.01304                # Able to execute with low frequencies as of now.
GEN = 40
EXPERIMENTS = [ # small pop, equal fitness, random mating
                {'constant_pop_size': 10000,      
                'aa_fitness'        : 1.0,
                'aa_homogamy'       : 0.0,
                'a_freq'            : A_FREQ,
                'gen'               : GEN}]
                # small pop, equal fitness, assortative mating
#                {'constant_pop_size' : 10000,     
#                 'aa_fitness'        : 1.0,
#                 'aa_homogamy'       : 0.9,
#                 'a_freq'            : A_FREQ,
#                 'gen'               : GEN},
#                 # small pop, aa homozygotes have 2x fitness, random mating
#                {'constant_pop_size' : 10000,      
#                 'aa_fitness'        : 1.2,
#                 'aa_homogamy'       : 0.0,
#                 'a_freq'            : A_FREQ,
#                 'gen'               : GEN},
#                 # small pop, aa homozygotes have 2x fitness, assortative mating
#                {'constant_pop_size' : 10000,      
#                 'aa_fitness'        : 1.2,
#                 'aa_homogamy'       : 0.9,
#                 'a_freq'            : A_FREQ,
#                 'gen'               : GEN}]

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
    pop = sim.Population(constant_pop_size, loci=1, infoFields=['fitness'])
    sim.initSex(pop)
    pop.dvars().header = []
    pop.dvars().row = []
    vsp = sim.GenotypeSplitter(loci=1, loci=1, alleles=[[0,0,0,1],[1,1]])
    pop.setVirtualSplitter(vsp)
    # Creates two virtual subpopulations needed in order to define a mating 
    # scheme: 
    #   Note: allele definitions are _unphased_ so (0,1) and (1,0) are equivalent
    #   alleles=[0,0,0,1] are individuals with 00 or 01 (AA/Aa)
    #   alleles=[1,1] are individuals with 11 (aa)
    pop.evolve(
        initOps= [sim.InitSex(subPops=[(0,0)], maleProp=0.5),
                  sim.InitSex(subPops=[(0,1)], maleProp=0.5),
                  # assigns individuals randomly to be male or female.
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
        matingScheme = sim.HeteroMating(matingSchemes=[sim.RandomMating(subPops=[(1,1)], weight=0-aa_homogamy),
                sim.RandomMating(weight=aa_homogamy)]),
                    
        
        postOps = [sim.Stat(popSize=True, alleleFreq=[0], subPops=[(0,0),(0,1)], 
                                genoFreq=[0]), 
                  # Addition of genoFreq to establish a parameter for counting individuals
                  # with a specific genotype. 
                   sim.PyExec(r"header += ['gen','AA/Aa_Size','aa_size','A_freq',"\
                                   "'a_freq', 'AA_Freq', 'aa_Freq', 'Aa_Freq', 'Aa_size',"\
                                   "'AA_size', 'FAA', 'FAa', 'Faa']"),
                  # Addition of Aa_size and AA_size. AA/Aa_Size has a capitalization
                  # to prevent error-reading mistakes.
                   sim.PyExec(r"row += [gen, subPopSize[0], subPopSize[1],"\
                                   "alleleFreq[0][0], alleleFreq[0][1],"\
                                   "genoFreq[0][(0,0)], genoFreq[0][(1,1)],"\
                                   "genoFreq[0][(0,1)]+genoFreq[0][(1,0)],"\
                                   "genoNum[0][(0,1)]+genoNum[0][(1,0)],"\
                                   "genoNum[0][(0,0)],"\
                                   "((genoFreq[0][(0,0)]-(alleleFreq[0][0]*"\
                                   "alleleFreq[0][0]))/(alleleFreq[0][0]-"\
                                   "(alleleFreq[0][0]*alleleFreq[0][0]))),"\
                                   "(1-((genoFreq[0][(0,1)]+genoFreq[0][(1,0)]))/"\
                                   "(2*alleleFreq[0][0]*alleleFreq[0][1])),"\
                                   "((genoFreq[0][(1,1)]-(alleleFreq[0][1]*"\
                                   "alleleFreq[0][1]))/(alleleFreq[0][1]-"\
                                   "(alleleFreq[0][1]*alleleFreq[0][1])))]")
                  # Addition of genoNum[x][(x,x)] to count the number of individuals
                  # with that specific genotype.
                  # You can add frequencies and sizes through addition. There are two
                  # sizes for both Aa carriers (0,1) and (1,0). Adding these can be done
                  # as seen above.                     
                   ],
        gen = gen
    )
    
    return {'header':pop.dvars().header, 'row':pop.dvars().row}


# MAIN FUNCTION
           
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
    #shutil.copyfile(source_fn, os.path.join(args.path,source_fn))
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
                    simuAssortativeMatingWithFitness with its parameters to the 
                    multiprocessing pool.
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
    multiplot_AA_Freqs = []
    multiplot_Aa_Freqs = []
    multiplot_aa_Freqs = []
    multiplot_A_freqs = []      # Addition of multiplot individual and frequency lists.
    multiplot_a_freqs = []
    multiplot_AA_Individ = []   
    multiplot_Aa_Individ = []   
    multiplot_aa_Individ = []   
    multiplot_params= []
    multiplot_titles = []
    n = 0                       # To ensure graphs are not overwritten.
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
        
        multiplot_params.append(params)
        # Select and type convert X, which should be the same for all data.
        X = []
        for h, col in zip(headers, data):
            if 'gen' in h:
                X.append(int(col[0]))
        # Select a_freq
        FAA = []
        FAa = []
        Faa = []
        A_freqs = []
        a_freqs = []
        AA_Freqs = []
        Aa_Freqs = []
        aa_Freqs = []     
        AA_Individ = []      # Addition of individual and frequency lists.
        Aa_Individ = []
        aa_Individ = [] 
         
        
        # Expansion of searching executives in order to find additional
        # data needed to create prove and test hypothesises, etc.
        for h, col in zip(headers, data):
            if 'A_freq' in h:
                A_freqs.append(col)
        for h, col in zip(headers, data):
            if 'a_freq' in h:
                a_freqs.append(col)
        for h, col in zip(headers, data):
            if 'AA_Freq' in h:
                AA_Freqs.append(col)
        for h, col in zip(headers, data):
            if 'Aa_Freq' in h:
                Aa_Freqs.append(col)
        for h, col in zip(headers, data):
            if 'aa_Freq' in h:
                aa_Freqs.append(col)
        for h, col in zip(headers, data):
            if 'AA_size' in h:
                AA_Individ.append(col)
        for h, col in zip(headers, data):
            if 'Aa_size' in h:
                Aa_Individ.append(col)
        for h, col in zip(headers, data):
            if 'aa_size' in h:
                aa_Individ.append(col)
        for h, col in zip(headers, data):
            if 'FAA' in h:
                FAA.append(col)   
        for h, col in zip(headers, data):
            if 'FAa' in h:
                FAa.append(col)
        for h, col in zip(headers, data):
            if 'Faa' in h:
                Faa.append(col)     
   
        multiplot_AA_Individ.append(AA_Individ)     # Additional Appending of individuals
        multiplot_Aa_Individ.append(Aa_Individ)     # and frequencies.
        multiplot_aa_Individ.append(aa_Individ)
        multiplot_A_freqs.append(A_freqs)
        multiplot_a_freqs.append(a_freqs)
        multiplot_AA_Freqs.append(AA_Freqs)
        multiplot_Aa_Freqs.append(Aa_Freqs)
        multiplot_aa_Freqs.append(aa_Freqs)

        title='N={:,}   '\
              'fitness={:.1f}    '\
              'homogamy={:.1f}'\
              ''.format(int(params['constant_pop_size']),
                        float(params['aa_fitness']),
                        float(params['aa_homogamy']))
        multiplot_titles.append(title)
        
        # Create individual (RECESSIVE ALLELE) contour charts
        for use in ['print']:
            n+=1
            rn = 'RecessiveAlleleFrequency'+str(n)+'.{}.pdf'.format(use)
            filename = os.path.join(args.path, rn)
            print "Saving chart to '{}'.".format(filename)
            gm.write_density_plot(filename, X, a_freqs,
                               title=title,
                               xlabel='Generations',
                               ylabel='Recessive Allele Frequency',
                               use=use)
                               
        # Create individual (CARRIER ALLELE) contour charts
        for use in ['print']:
            dn = 'CarrierGeneFrequency'+str(n)+'.{}.pdf'.format(use)
            filename = os.path.join(args.path, dn)
            print "Saving chart to '{}'.".format(filename)
            gm.write_density_plot(filename, X, Aa_Freqs,
                               title=title,
                               xlabel='Generations',
                               ylabel='Carrier Gene Frequency',
                               use=use)
                               
        # Create individual (DOMINANT ALLELE) contour charts
        for use in ['print']:
            dn = 'DominantAlleleFrequency'+str(n)+'.{}.pdf'.format(use)
            filename = os.path.join(args.path, dn)
            print "Saving chart to '{}'.".format(filename)
            gm.write_density_plot(filename, X, A_freqs,
                               title=title,
                               xlabel='Generations',
                               ylabel='Dominant Allele Frequency',
                               use=use)
                               
                               
        # Create individual (RECESSIVE AND CARRIER INDIVIDUAL COMPARISON) contour charts
        for use in ['print']:
            fn = 'RVSCIndividualFrequency'+str(n)+'.{}.pdf'.format(use)
            filename = os.path.join(args.path, fn)
            print "Saving chart to '{}'.".format(filename)
            Ylst = [aa_Individ, Aa_Individ]
            gm.write_simple_plot(filename, X, Ylst,
                                    title=title,
                                    xlabel='Generations',
                                    ylabel='R VS C Number of Individuals',
                                    use=use)
                                  
                                    
        # Create individual (DOMINANT INDIVIDUAL) contour charts
        for use in ['print']:
            fn = 'DominantIndividualFrequency'+str(n)+'.{}.pdf'.format(use)
            filename = os.path.join(args.path, fn)
            print "Saving chart to '{}'.".format(filename)
            Ylst = [AA_Individ]
            gm.write_simple_plot(filename, X, Ylst,
                                    title=title,
                                    xlabel='Generations',
                                    ylabel='Dominant Number of Individuals',
                                    use=use)
                                
       
    # Create summary (RECESSIVE FREQUENCY) contour charts
    for use in ['print']:
        bn = 'summaryRecessiveGeneFrequency.{}.pdf'.format(use)
        filename = os.path.join(args.path, bn)
        print "Saving summary chart to '{}'.".format(filename)
        gm.write_summary_density_plot(filename, 
                                   Xarr=[X for i in range(len(multiplot_a_freqs))],
                                   Yarr=multiplot_a_freqs,
                                   nrows=2,
                                   ncols=2,
                                   title='Comparison Between Recessive Gene Frequencies',
                                   multiplot_titles=multiplot_titles,
                                   xlabel='Generations',
                                   ylabel='Recessive Gene Frequency',
                                   use = use)
                                   
                                   
    # Create summary (DOMINANT FREQUENCY) contour charts
    for use in ['print']:
        bn1 = 'summaryDominantGeneFrequency.{}.pdf'.format(use)
        filename = os.path.join(args.path, bn1)
        print "Saving summary chart to '{}'.".format(filename)
        gm.write_summary_density_plot(filename, 
                                   Xarr=[X for i in range(len(multiplot_A_freqs))],
                                   Yarr=multiplot_A_freqs,
                                   nrows=2,
                                   ncols=2,
                                   title='Comparison Between Dominant Gene Frequencies',
                                   multiplot_titles=multiplot_titles,
                                   xlabel='Generations',
                                   ylabel='Dominant Gene Frequency',
                                   use = use)
        
        
    # Create summary (DOMINANT INDIVIDUAL) contour charts
    for use in ['print']:
        bn1 = 'summaryDominantIndividuals.{}.pdf'.format(use)
        filename = os.path.join(args.path, bn1)
        print "Saving summary chart to '{}'.".format(filename)
        gm.write_summary_simple_plot(filename, 
                                   Xarr=[X for i in range(len(multiplot_AA_Individ))],
                                   Yarr=multiplot_AA_Individ,
                                   nrows=2,
                                   ncols=2,
                                   title='Comparison Between Dominant Individual Frequencies',
                                   multiplot_titles=multiplot_titles,
                                   xlabel='Generations',
                                   ylabel='Dominant Number of Individuals',
                                   use = use)
                                
    # Create summary (RECESSIVE AND CARRIER INDIVIDUAL COMPARISON) contour charts
    for use in ['print']:
        bn1 = 'summaryRecessiveAndCarrier.{}.pdf'.format(use)
        filename = os.path.join(args.path, bn1)
        print "Saving summary chart to '{}'.".format(filename)
        Ya = [multiplot_Aa_Individ, multiplot_aa_Individ]
        Ya = Ya[0]+Ya[1]
        gm.write_summary_simple_plot(filename, 
                                   Xarr=[X for i in range(len(multiplot_AA_Individ))],
                                   Yarr=Ya,
                                   nrows=2,
                                   ncols=2,
                                   title='Comparison Between R VS C Individual Frequencies',
                                   multiplot_titles=multiplot_titles,
                                   xlabel='Generations',
                                   ylabel='R VS C Number of Individuals',
                                   use = use)

    # Create summary (DOMINANT/CARRIER/RECESSIVE FREQUENCY COMPARISON) contour charts
    for use in ['print']:
        bn1 = 'summaryDominant_Carrier_Recessive_Frequency_Comparison.{}.pdf'.format(use)
        filename = os.path.join(args.path, bn1)
        print "Saving summary chart to '{}'.".format(filename)
        Ya = [multiplot_AA_Freqs, multiplot_Aa_Freqs, multiplot_aa_Freqs]
        Ya = Ya[0]+Ya[1]+Ya[2]
        gm.write_summary_simple_plot(filename, 
                                   Xarr=[X for i in range(len(multiplot_AA_Individ))],
                                   Yarr=Ya,
                                   nrows=2,
                                   ncols=2,
                                   title='Comparison Between R VS C Genotype Frequencies',
                                   multiplot_titles=multiplot_titles,
                                   xlabel='Generations',
                                   ylabel='R VS C Geneotype Frequency',
                                   use = use)                              

                                   