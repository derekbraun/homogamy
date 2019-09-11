#!/usr/local/bin/python3 -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
Samir Jain, Eric Epstein, Maggie Gray, Derek Braun*
(*derek.braun@gallaudet.edu)

Simulation module which uses simuPOP.
Simulation parameters are set via command-line arguments.

Last updated: 2-May-2019 by Derek Braun
'''

# Default Simulation Parameters from Nance and Kearsey (2004)
#   First 5 generations: fitness 0
#   Homogamy increases to 90% by gen 5
#   20 more generations (400 years) representing 1800-2200
#   Population size fixed at 200,000

# defaults; some of these can be set otherwise by passing arguments
CONSTANT_POP_SIZE   = 200      # Nance and Kearsey: 200k
a_FREQ              = 0.01304  # Nance and Kearsey: 0.01304
aa_FITNESS          = 1.0      # Nance and Kearsey: 1.0
aa_HOMOGAMY         = 0.9      # Nance and Kearsey: 0.9
DEAF_FREQ           = 0.0008   # Nance and Kearsey: two other genes each at 0.003^2 freq (tiny)
GENERATIONS         = 20       # Nance and Kearsey: 5 gen with 0 fitness + 20 gen with 1.0 fitness
SIMULATIONS         = 5000


import os
import time
import random
import argparse
import subprocess
import multiprocessing
import simuOpt
simuOpt.setOptions(optimized=True, numThreads=0, quiet=True)
import simuPOP as sim
import fileio


def customChooser(pop, subPop):
    '''
        Generator function which chooses parents.

        Upon initialization, this chooser mates couples in a monogamous
        mating scheme.

        The algorithm goes through each eligible person listwise.
        If that person is deaf, that person marries deaf or hearing based on the
        probability set by aa_homogamy. Further, if either
        parent is deaf, their number of offspring is based on the
        fitness set by aa_fitness (non-integer fitnesses are handled using a
        randomizer for the fractional amount).

        Uses:  (there is no way to pass variables to customChooser)
        pop.dvars().constant_pop_size
        pop.dvars().a
        pop.dvars().aa_fitness
        pop.dvars().aa_homogamy
        pop.dvars().deaf

        Accepts:
        pop, subpop (this is standard/required by simuPOP)

        Yields:
        (parent1, parent2) tuple of two individuals (standard/required by simuPOP)
    '''

    def mate_with_fitness(parent1, parent2):
        '''
            Mates a couple and creates offspring. Creates a number of entries in the
            final list (representing the parents for each child to be born), based on
            reproductive fitness. Non-integer number of children are handled by
            using a randomizer for the fractional amount.

            Accepts:
            parent1, parent2             sim.individual objects

            Returns a list of (parent1, parent2) parental pairs
            reflecting fitness.
        '''

        if parent1.genotype() == [1,1] or parent2.genotype() == [1,1]:
            r = float(pop.dvars().aa_fitness)
            l = []
            while r >= 1:
                l += [(parent1, parent2)]
                r -= 1
            if random.random() < r:
                l += [(parent1, parent2)]
            return l
        else:
            return [(parent1, parent2)]

    deaf_parents = []
    hearing_parents = []
    couples = []

    # bin individuals
    for person in pop.individuals():
        if person.genotype() == [1,1]:
            deaf_parents.append(person)
        else:
            hearing_parents.append(person)

    # move some "hearing" individuals into the deaf bin - making them deaf -
    # to reflect non-Cx26 causes of congenital deafness. These individuals will
    # mate with other deaf but will not pass down Cx26
    for i in range(pop.dvars().adv_deaf_target):
        if len(hearing_parents) > 0:
            deaf_parents.append(hearing_parents.pop())
        else:
            break
    random.shuffle(deaf_parents)

    # calculate how many deaf-deaf marriages we need, then marry them off
    dp = float(len(deaf_parents))
    target = int(round(pop.dvars().aa_homogamy * len(deaf_parents)/2))
    for i in range(target):
        if len(deaf_parents) >= 2:
            couples += mate_with_fitness(deaf_parents.pop(), deaf_parents.pop())
        else:
            break
    if dp > 0:
        pop.dvars().homogamy = 2*target/dp
    else:
        pop.dvars().homogamy = -1

    # Merge remaining parents, so that their alleles are not lost.
    # Then, mate them off
    remaining_parents = hearing_parents + deaf_parents
    random.shuffle(remaining_parents)
    while len(remaining_parents) > 2:
        couples += mate_with_fitness(remaining_parents.pop(), remaining_parents.pop())

    # This is what's called whenever the generator function is called.
    while True:
        yield random.choice(couples)


def simuAssortativeMatingWithFitness(e):
    '''
        Accepts:
        e               an Experiment object.

        Returns a dict containing the results from each gen of the simulation:
        gen             generation number.
        A               frequency of the A allele.
        a               frequency of the a allele.
        AA              frequency of AA individuals.
        Aa              frequency of Aa individuals.
        aa              frequency of aa individuals.
        deaf            frequency of deaf individuals (incl adventitious).
        AA_size         size of the AA subpopulation.
        Aa_size         size of the Aa subpopulation.
        aa_size         size of the aa subpopulation.
        deaf_size       size of the deaf subpopulation (incl adventitious).
        homogamy        calculated actual homogamy.
        F               calculated inbreeding coefficient.

        Adopted from: http://simupop.sourceforge.net/Cookbook/AssortativeMating
    '''
    sim.setRNG(random.seed(sim.getRNG().seed()))
    pop = sim.Population(e.constant_pop_size*1000, loci=[1])
    # These variables need to be set in order to be available to customChooser().
    # There appears to be no way to directly pass variables to customChooser().
    pop.dvars().constant_pop_size   = e.constant_pop_size
    pop.dvars().a                   = e.a
    pop.dvars().aa_fitness          = e.aa_fitness
    pop.dvars().aa_homogamy         = e.aa_homogamy
    pop.dvars().deaf                = e.deaf
    pop.dvars().adv_deaf_target     = int(round((e.deaf - e.a**2) * e.constant_pop_size * 1000))

    # These will hold the final data
    pop.dvars().headers             = []
    pop.dvars().row                 = []
    pop.evolve(
        initOps= [sim.InitGenotype(freq=[1-e.a, e.a])],
        matingScheme = sim.HomoMating(
                    chooser = sim.PyParentsChooser(customChooser),
                    generator = sim.OffspringGenerator(sim.MendelianGenoTransmitter())),
        postOps = [sim.Stat(alleleFreq=[0], genoFreq=[0]),
                   sim.PyExec(r"headers += ['gen','A', 'a',"\
                               "'AA', 'Aa', 'aa', 'deaf', 'AA_size', 'Aa_size', " \
                               "'aa_size',  'deaf_size', 'homogamy', 'F'] \n" \
                               "F = 1.0-((genoFreq[0][(0,1)]+genoFreq[0][(1,0)])/" # F          \
                               "(2.0*alleleFreq[0][0]*alleleFreq[0][1])) "\
                               "if alleleFreq[0][0]*alleleFreq[0][1] > 0. "\
                               "else 0. \n" \
                               "deaf_size = min(genoNum[0][(1,1)] + adv_deaf_target, constant_pop_size*1000) \n"\
                               "row += [gen, "                           # generation \
                               "alleleFreq[0][0], "                      # A          \
                               "alleleFreq[0][1], "                      # a          \
                               "genoFreq[0][(0,0)],"                     # AA         \
                               "genoFreq[0][(0,1)]+genoFreq[0][(1,0)], " # Aa         \
                               "genoFreq[0][(1,1)], "                    # aa         \
                               "deaf_size/(constant_pop_size*1000.), "   # deaf       \
                               "genoNum[0][(0,0)], "                     # AA_size    \
                               "genoNum[0][(0,1)]+genoNum[0][(1,0)], "   # Aa_size    \
                               "genoNum[0][(1,1)], "                     # aa_size    \
                               "deaf_size, "                             # deaf_size  \
                               "homogamy, "                              # homogamy   \
                               "F if F>0. else 0.]")                     # F          \
                   ],
        gen = e.generations
    )
    return {'headers':pop.dvars().headers, 'row':pop.dvars().row}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path',
                        nargs = '?',  # makes this argument optional
                        default = os.getcwd(),
                        help = 'results folder path. If the folder does not '\
                               'exist, it will be created.')
    parser.add_argument('-w','--write',
                        action = 'store_true',
                        help = 'run {:,} simulations and write to file.' \
                               ''.format(SIMULATIONS))
    parser.add_argument('-p', '--pop_size',
                        action = 'store',
                        type = float,
                        default = CONSTANT_POP_SIZE,
                        help = 'constant population size, in thousands. ' \
                               '(default {}).'.format(CONSTANT_POP_SIZE))
    parser.add_argument('--homogamy',
                        action = 'store',
                        type = float,
                        default = aa_HOMOGAMY,
                        help = 'deaf-deaf assortative mating (homogamy) ' \
                               '(default {}).'.format(aa_HOMOGAMY))
    parser.add_argument('-f', '--fitness',
                        action = 'store',
                        type = float,
                        default = aa_FITNESS,
                        help = 'the relative reproductive fitness of deaf ' \
                               'individuals (default {}).'.format(aa_FITNESS))
    args=parser.parse_args()

    experiment = fileio.Experiment(constant_pop_size   = args.pop_size,
                                   a                   = a_FREQ,
                                   aa_fitness          = args.fitness,
                                   aa_homogamy         = args.homogamy,
                                   deaf                = DEAF_FREQ,
                                   generations         = GENERATIONS)

    # This quick run obtains the headers for the data file.
    sample_run = simuAssortativeMatingWithFitness(experiment)
    experiment.cpu = subprocess.check_output(['/usr/sbin/sysctl', "-n", \
                                     "machdep.cpu.brand_string"]).decode().strip() + \
                                     " ({} threads)".format(multiprocessing.cpu_count())
    experiment.simuPOP_version = sim.__version__

    if not args.write:
        # just show the results from the quick sample run and exit
        print(experiment.metadata())
        numcols = sample_run['headers'][1:].index("gen") + 1
        for h in sample_run['headers'][0:numcols]:
            print("{h:>9}".format(h=h), end=' ')
        print()
        for h in sample_run['headers'][0:numcols]:
            print(" --------", end=' ')
        print()
        for gen in range(0,len(sample_run['row'])//numcols):
            for datum in sample_run['row'][numcols*gen:numcols*(1+gen)]:
                if type(datum) is int or datum == int(datum):
                    print(" {datum:>8,}".format(datum=int(datum)), end=' ')
                else:
                    print(" {datum:>8.6f}".format(datum=datum), end=' ')
            print()
        print('Done.')
        exit()
    else:
        if fileio.create_folder(args.path):
            print('Creating folder...\n   {}'.format(args.path))
        experiment.headers = sample_run['headers']
        experiment.filename = os.path.join(args.path,
                                           'pop{experiment.constant_pop_size}k'\
                                           '_fit{experiment.aa_fitness}'       \
                                           '_hom{experiment.aa_homogamy:.2}'   \
                                           '.tsv'.format(**locals()))
        experiment.write_metadata(overwrite=True)
        print('Creating file...\n   {}'.format(experiment.filename))
        print(experiment.metadata())
        print('Running {:,} simulations...'.format(SIMULATIONS))

        def _worker():
            '''
                The worker function exists as a convenient way of passing
                simuAssortativeMatingWithFitness with its parameters to the
                multiprocessing pool.
            '''
            return simuAssortativeMatingWithFitness(experiment)['row']

        def _format_time (time):
            h = time//3600
            m = (time - 3600*(time//3600))//60
            s = time%60
            if h:
                return '{:.0f}h {:.0f}m'.format(h, m)
            elif m:
                return '{:.0f}m {:.0f}s'.format(m, s)
            else:
                return '{:.1f}s'.format(s)


        mp_chunk_size = cpu_count = multiprocessing.cpu_count()
        pool = multiprocessing.Pool()
        sims = 0
        while sims < SIMULATIONS:
            start_time = time.time()
            p = [pool.apply_async(_worker) for i in range(mp_chunk_size)]
            table = [item.get() for item in p]
            experiment.write(table)
            sims += mp_chunk_size
            rate = mp_chunk_size/(time.time()-start_time)
            time_remaining = (SIMULATIONS-sims)/rate if rate > 0 else 0
            print('   {:,} simulations completed ' \
                  '({:,.0f}/min) '\
                  '{} remaining.'\
                  ''.format(sims, 60*rate, _format_time(time_remaining)))
            # mp_chunk_size is dynamically adjusted based on actual
            # execution speed such that file writes occur once per minute.
            mp_chunk_size = max(int(300*rate - 300*rate%cpu_count), cpu_count)
            if sims + mp_chunk_size > SIMULATIONS:
                mp_chunk_size = SIMULATIONS-sims
        print('Saving file...\n   {}'.format(experiment.filename))
        exit()
