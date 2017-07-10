#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
    Samir Jain, Eric Epstein, Trevor Klemp, Maggie Gray, Selman Jawed, Derek 
    Braun* (*derek.braun@gallaudet.edu)
    
    Simulation module which uses simuPOP. Simulation parameters are set via 
    setting global variables, and, for changing assortative mating, sometimes
    by changing the code.
    
    Last updated: 23-Jun-2017 by Derek Braun
'''


DEBUG_MODE = False

# Simulation Parameters
SIMULATIONS = 100
aa_HOMOGAMY = 0.0               # This variable MUST be a global b/c this is 
                                # the only way to get it into the deafChooser
                                # generator function
                                
import fileio
experiment = fileio.Experiment( constant_pop_size   = 10000,
                                aa_fitness          = 20.0,
                                aa_homogamy         = aa_HOMOGAMY,
                                a                   = 0.01304,
                                gen                 = 100)

import os
import time
import random
import argparse
import multiprocessing
import simuOpt
if DEBUG_MODE:
    PROPOSALS = 1
else:
    simuOpt.setOptions(optimized=True, numThreads=0, quiet=True)
import simuPOP as sim

def deafChooser(pop, subPop):
    '''
        Generator function which chooses parents. 
        I don't know if it's possible to pass additional parameters to 
        this generator but I would like to pass the percentage of
        deaf-deaf marriages (aa_homogamy).
    
        Upon initialization, this chooser pairs up (marries) couples,
        resulting in a monogamous mating scheme. Deaf are paired up first,
        so as to achieve the desired percentage of deaf-deaf marriages.
    
        Each time this generator is called, it returns a random couple.
        The couples do not change within a generation (there is no divorce or
        death). Implemented this way, roughly 80% of couples will have
        children, and roughly 20% will have more than one child.
    
        A monogamous mating scheme isn't entirely representative of human 
        behavior, but it's much closer to reality than an entirely random
        mating scheme where nearly every child will have different parents and
        there are almost no full siblings.
    '''
    all_males = []
    all_females = []
    hearing_males = []
    hearing_females = []
    deaf_males = []
    deaf_females = []
    couples = []
    
    # bin individuals
    for i in range(pop.subPopSize(subPop)):
        person = pop.individual(i)
        if person.sex() == sim.MALE:
            all_males.append(i)
            if list(person.genotype()) == [1, 1]:
                deaf_males.append(i)
            else:
                hearing_males.append(i)
        elif person.sex() == sim.FEMALE:
            all_females.append(i)
            if list(person.genotype()) == [1, 1]:
                deaf_females.append(i)
            else:
                hearing_females.append(i)
        else:
            print "simuPOP gender dysphoria error. Send scathing email to Bo Peng"


    # pair off deaf individuals first, to achieve the desired percentage
    # of deaf-deaf marriage
    random.shuffle(deaf_females)
    random.shuffle(deaf_males)
    random.shuffle(hearing_females)
    random.shuffle(hearing_males)
    while len(deaf_females) > 0 and len(deaf_males) > 0 and len(hearing_males) > 0:
        woman = deaf_females.pop()
        if random.random() < aa_HOMOGAMY:
            man = deaf_males.pop()
            couples += [(man, woman)]
        else:
            man = hearing_males.pop()
            couples += [(man, woman)]
    
    # move remaining deaf people into hearing (now general) bins, and reshuffle
    hearing_males += deaf_males        
    hearing_females += deaf_females
    random.shuffle(hearing_males)
    random.shuffle(hearing_females)
    
    while len(hearing_females) > 0 and len(hearing_males) > 0:
        woman = hearing_females.pop()
        man = hearing_males.pop()
        couples += [(man, woman)]
        
    # marry off remaining un-married, mostly hearing people 
    # (let's just say that this represents second marriages and out-of-wedlock
    # children. More importantly, why are we doing this? 
    # This is needed to nullify the slight advantage that we
    # just gave deaf people by making sure that every single deaf person was
    # married, which we didn't do for hearing people. Without this extra code,
    # in small populations, deaf people will have a slight fitness advantage.
    
    while len(hearing_females) > 0:
        women = hearing_females.pop()
        man = random.choice(all_males)
        couples += [(man, woman)]
    while len(hearing_males) > 0:
        man = hearing_males.pop()
        woman = random.choice(all_females)
        couples += [(man, woman)]

    # This is what's called whenever the generator function is called,
    # after the first time.
    # Hey, you! You're having a kid today!
    while True:
        yield random.choice(couples)

        
def simuAssortativeMatingWithFitness(constant_pop_size, gen, a, 
                                    aa_fitness, aa_homogamy):
    '''
        Accepts:
        constant_pop_size   population size, which remains constant throughout
        gen                 number of generations
        a                   starting frequency of the a allele
        aa_fitness          _relative_ fitness of deaf (aa) individuals
        aa_homogamy         the percent of assortative mating between
                            deaf individuals
        Returns a dict containing the results from the simulation:
        gen                 generation number
        AA/Aa_size          size of the AA/aa population
        aa_size             size of the aa population
        A                   frequency of the A allele
        a                   frequency of the a allele
        
        Adopted from: http://simupop.sourceforge.net/Cookbook/AssortativeMating
    '''             
    sim.setRNG(random.seed(sim.getRNG().seed()))    
    pop = sim.Population(constant_pop_size, loci=[1], infoFields=['fitness'])
    pop.dvars().headers = [] 
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
                  sim.InitGenotype(freq=[1-a, a])
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
        matingScheme = sim.HomoMating(
                        chooser = sim.PyParentsChooser(deafChooser),
                        generator = sim.OffspringGenerator(sim.MendelianGenoTransmitter())),
                    
        postOps = [sim.Stat(alleleFreq=[0], genoFreq=[0]), 
                   sim.PyExec(r"headers += ['gen','A', 'a',"\
                               "'AA', 'Aa', 'aa',"\
                               "'AA_size', 'Aa_size', 'aa_size',"\
                               "'F']"),
                   sim.PyExec(r"row += [gen, "\
                               "alleleFreq[0][0],"                             # A          \
                               "alleleFreq[0][1],"                             # a          \
                               "genoFreq[0][(0,0)],"                           # AA         \
                               "genoFreq[0][(0,1)]+genoFreq[0][(1,0)],"        # Aa         \
                               "genoFreq[0][(1,1)],"                           # aa         \
                               "genoNum[0][(0,0)],"                            # AA_size    \
                               "genoNum[0][(0,1)]+genoNum[0][(1,0)],"          # Aa_size    \
                               "genoNum[0][(1,1)],"                            # aa_size    \
                               "1.0-((genoFreq[0][(0,1)]+genoFreq[0][(1,0)])/" # F          \
                               "(2.0*alleleFreq[0][0]*alleleFreq[0][1]))"\
                               "if alleleFreq[0][0] != 0.0 and alleleFreq[0][1]"\
                               "!= 0.0 else 0.0]")
                   ],
        gen = gen
    )
    return {'headers':pop.dvars().headers, 'row':pop.dvars().row}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path',
                        help = 'results folder path. If the folder does not exist, '\
                               'it will be created',
                        nargs = '?',  # makes this argument optional
                        default = os.getcwd())
    parser.add_argument('-o','--overwrite',action='store_true',
                        help = 'overwrite tsv file')
    parser.add_argument('-s','--sample_run',action='store_true',
                        help = 'sample run only; display results')
    args=parser.parse_args()
    
    # This quick sample run obtains the headers for the data file.
    sample_run = simuAssortativeMatingWithFitness(experiment.constant_pop_size, 
                                                  experiment.gen,
                                                  experiment.a, 
                                                  experiment.aa_fitness,   
                                                  experiment.aa_homogamy)    
    if args.sample_run:    
        print experiment.metadata()
        for h in sample_run['headers'][0:10]:
            print "{h:>8}".format(h=h),
        print
        for h in sample_run['headers'][0:10]:
            print " -------",
        print
        for gen in range(0,len(sample_run['row'])/10):
            for datum in sample_run['row'][10*gen:10*(1+gen)]:
                if type(datum) is int or datum == int(datum):
                    print " {datum:>7,}".format(datum=int(datum)),
                else:
                    print " {datum:>7.5f}".format(datum=datum),
            print
    else:
        if fileio.create_folder(args.path):
            print "  Created folder '{}'".format(args.path)
        else:
            print "  Using folder '{}'".format(args.path)
        experiment.headers = sample_run['headers']
        experiment.source_code = os.path.split(__file__)[-1].replace('.pyc','.py')
        experiment.filename = os.path.join(args.path, 
                                           'pop{experiment.constant_pop_size}'\
                                           '_fitness{experiment.aa_fitness}'\
                                           '_homogamy{experiment.aa_homogamy:.2}.tsv'\
                                           ''.format(**locals()))
        if experiment.write_metadata():
            print "  Created '{}'".format(experiment.filename)
        elif args.overwrite:
            experiment.write_metadata(overwrite=True)
            print "  Overwrote '{}'".format(experiment.filename)
        else:
            print "  File '{}' exists.".format(experiment.filename)
            print "  Use --overwrite to overwrite this file."
            exit()
        
        print "  Running simulations..."
    
        def worker():
                '''
                    The worker function exists as a convenient way of passing
                    simuAssortativeMatingWithFitness with its parameters to the 
                    multiprocessing pool.
                '''
                return simuAssortativeMatingWithFitness(experiment.constant_pop_size, 
                                                        experiment.gen,
                                                        experiment.a, 
                                                        experiment.aa_fitness,   
                                                        experiment.aa_homogamy)['row']
        
        simulations = 0
        mp_chunk_size = cpu_count = multiprocessing.cpu_count()
        pool = multiprocessing.Pool()
        while simulations < SIMULATIONS:
            start_time = time.time()
            p = [pool.apply_async(worker) for i in range(mp_chunk_size)]
            table = [item.get() for item in p]
            experiment.write(table)
            simulations += mp_chunk_size
            rate = mp_chunk_size*60./(time.time()-start_time)
            print '  {simulations:,} simulations completed ' \
                  '({cpu_count} CPUs; {rate:,.1f} simulations/min)'\
                  ''.format(simulations=simulations,
                            rate=rate,
                            cpu_count=cpu_count)
            # mp_chunk_size is dynamically adjusted based on actual
            # execution speed such that file writes occur once per minute.
            mp_chunk_size = int(rate - rate%cpu_count)
            if simulations + mp_chunk_size > SIMULATIONS:
                mp_chunk_size = SIMULATIONS-simulations
    print '  Done.'