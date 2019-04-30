#!/usr/local/bin/python -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
    Samir Jain, Eric Epstein, Maggie Gray, Derek Braun*
    (*derek.braun@gallaudet.edu)

    Simulation module which uses simuPOP.
    Simulation parameters are set via setting the experiment variables below.

    Last updated: 2-Apr-2019 by Derek Braun
'''


# Simulation Parameters
import fileio
experiment = fileio.Experiment( constant_pop_size   = 200,      # Nance and Kearsey: 200k
                                a                   = 0.01304,  # Nance and Kearsey: 0.01304
                                aa_fitness          = 1.5,      # Nance and Kearsey: 1.0
                                aa_homogamy         = 0.9,      # Nance and Kearsey: 0.9
                                deafness_freq       = 0.0008,   # Nance and Kearsey: look again
                                generations         = 20,       # Nance and Kearsey: 5 gen with 0 fitness + 20 gen with 1.0 fitness
                                simulations         = 5000)     # Nance and Kearsey: unspecified
import os
import time
import random
import argparse
import subprocess
import multiprocessing
import simuOpt
simuOpt.setOptions(optimized=True, numThreads=0, quiet=True)
import simuPOP as sim

def customChooser(pop, subPop):
    '''
        Generator function which chooses parents.

        Upon initialization, this chooser mates couples in a monogamous
        mating scheme.

        The algorithm goes through each eligible person listwise.
        If that person is deaf, that person marries deaf or not based on the
        probability set by aa_homogamy. Further, if either
        parent is deaf, their number of offspring is based on the
        fitness set by aa_fitness (non-integer fitnesses are handled using a
        randomizer for the fractional amount).
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
    adv_deaf = int(round((pop.dvars().deafness_freq - pop.dvars().a**2) * pop.dvars().constant_pop_size * 1000))
    while adv_deaf > 0 and len(hearing_parents) > 0:
        deaf_parents.append(hearing_parents.pop())
        adv_deaf -= 1
    random.shuffle(deaf_parents)
    pop.dvars().num_deaf = len(deaf_parents)

    # calculate how many deaf-deaf marriages we need, then marry them off
    homogamy_target = round(pop.dvars().aa_homogamy * len(deaf_parents)/2)
    ddm = 0.                            # deaf-deaf marriages
    dp = float(len(deaf_parents))       # deaf parents
    while len(deaf_parents) and homogamy_target > 0:
        couples += mate_with_fitness(deaf_parents.pop(), deaf_parents.pop())
        ddm += 1.
        homogamy_target -= 1
    pop.dvars().homogamy = 2.*ddm/dp

    # Move remaining deaf parents into the hearing bin, so that their alleles
    # are not lost. Then, mate off the rest
    hearing_parents += deaf_parents
    random.shuffle(hearing_parents)
    while len(hearing_parents):
        couples += mate_with_fitness(hearing_parents.pop(), hearing_parents.pop())

    # This is what's called whenever the generator function is called.
    while True:
        yield random.choice(couples)


def simuAssortativeMatingWithFitness(constant_pop_size, a, aa_fitness,
                                     aa_homogamy, deafness_freq, generations):
    '''
        Accepts:
        constant_pop_size   population size, which remains constant throughout.
        a                   starting frequency of the a allele.
        aa_fitness          _relative_ fitness of deaf (aa) individuals.
        aa_homogamy         the percent of assortative mating between
                            deaf individuals.
        deafness_freq       overall frequency of deaf individuals at the time of
                            reproductive age, including from causes other than
                            connexin 26.
        generations         number of generations.

        Returns a dict containing the results from each gen of the simulation:
        gen                 generation number.
        A                   frequency of the A allele.
        a                   frequency of the a allele.
        AA                  frequency of AA individuals.
        Aa                  frequency of Aa individuals.
        aa                  frequency of aa individuals.
        AA_size             size of the AA subpopulation.
        Aa_size             size of the Aa subpopulation.
        aa_size             size of the aa subpopulation.
        F                   calculated inbreeding coefficient.
        num_deaf            size of the deaf subpopulation (incl adventitious).
        homogamy            calculated actual homogamy.

        Adopted from: http://simupop.sourceforge.net/Cookbook/AssortativeMating
    '''
    sim.setRNG(random.seed(sim.getRNG().seed()))
    pop = sim.Population(constant_pop_size*1000, loci=[1])
    pop.dvars().constant_pop_size   = constant_pop_size
    pop.dvars().a                   = a
    pop.dvars().aa_fitness          = aa_fitness
    pop.dvars().aa_homogamy         = aa_homogamy
    pop.dvars().deafness_freq       = deafness_freq
    pop.dvars().generations         = generations
    pop.dvars().headers             = []
    pop.dvars().row                 = []
    pop.evolve(
        initOps= [sim.InitSex(),
                  # Assigns individuals randomly to be male or female.
                  # This can result in slightly more males or females,
                  # which can cause errors if the wrong mating scheme is
                  # selected.
                  sim.InitGenotype(freq=[1-a, a])
                  ],
        matingScheme = sim.HomoMating(
                        chooser = sim.PyParentsChooser(customChooser),
                        generator = sim.OffspringGenerator(sim.MendelianGenoTransmitter())),
        postOps = [sim.Stat(alleleFreq=[0], genoFreq=[0]),
                   sim.PyExec(r"headers += ['gen','A', 'a',"\
                               "'AA', 'Aa', 'aa',"\
                               "'AA_size', 'Aa_size', 'aa_size',"\
                               "'F', 'nm_deaf', 'homogmy']"),
                   sim.PyExec(r"F = 1.0-((genoFreq[0][(0,1)]+genoFreq[0][(1,0)])/" # F          \
                              "(2.0*alleleFreq[0][0]*alleleFreq[0][1]))"\
                              "if alleleFreq[0][0] != 0.0 and alleleFreq[0][1]"\
                              "!= 0.0 else 0.0"),
                   sim.PyExec(r"row += [gen, "\
                               "alleleFreq[0][0],"                             # A          \
                               "alleleFreq[0][1],"                             # a          \
                               "genoFreq[0][(0,0)],"                           # AA         \
                               "genoFreq[0][(0,1)]+genoFreq[0][(1,0)],"        # Aa         \
                               "genoFreq[0][(1,1)],"                           # aa         \
                               "genoNum[0][(0,0)],"                            # AA_size    \
                               "genoNum[0][(0,1)]+genoNum[0][(1,0)],"          # Aa_size    \
                               "genoNum[0][(1,1)],"                            # aa_size    \
                               "F if F>0.0 else 0.0,"                          # F          \
                               "num_deaf,"                                      # nm_deaf    \
                               "homogamy]")                                     # homogmy
                   ],
        gen = generations
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
                                                  experiment.a,
                                                  experiment.aa_fitness,
                                                  experiment.aa_homogamy,
                                                  experiment.deafness_freq,
                                                  experiment.generations)
    experiment.source_code = os.path.split(__file__)[-1].replace('.pyc','.py')
    experiment.cpu = subprocess.check_output(['/usr/sbin/sysctl', "-n", \
                                     "machdep.cpu.brand_string"]).strip() + \
                                     " ({} threads)".format(multiprocessing.cpu_count())


    if args.sample_run:
        print experiment.metadata()
        for h in sample_run['headers'][0:12]:
            print "{h:>8}".format(h=h),
        print
        for h in sample_run['headers'][0:12]:
            print " -------",
        print
        for gen in range(0,len(sample_run['row'])/12):
            for datum in sample_run['row'][12*gen:12*(1+gen)]:
                if type(datum) is int or datum == int(datum):
                    print " {datum:>7,}".format(datum=int(datum)),
                else:
                    print " {datum:>7.5f}".format(datum=datum),
            print
    else:
        if fileio.create_folder(args.path):
            print "Created folder '{}'".format(args.path)
        else:
            print "Using folder '{}'".format(args.path)
        experiment.headers = sample_run['headers']
        experiment.source_code = os.path.split(__file__)[-1].replace('.pyc','.py')
        experiment.filename = os.path.join(args.path,
                                           'pop{experiment.constant_pop_size}k'\
                                           '_fit{experiment.aa_fitness}'\
                                           '_hom{experiment.aa_homogamy:.2}'\
                                           '.tsv'.format(**locals()))
        if experiment.write_metadata():
            print "Created '{}'".format(experiment.filename)
        elif args.overwrite:
            experiment.write_metadata(overwrite=True)
            print "Overwrote '{}'".format(experiment.filename)
        else:
            print "'{}' exists.".format(experiment.filename)
            print "  Use --overwrite to overwrite this file."
            exit()
        print experiment.metadata()
        print "Running {:,} simulations...".format(experiment.simulations)

        def _worker():
                '''
                    The worker function exists as a convenient way of passing
                    simuAssortativeMatingWithFitness with its parameters to the
                    multiprocessing pool.
                '''
                return simuAssortativeMatingWithFitness(experiment.constant_pop_size,
                                                        experiment.a,
                                                        experiment.aa_fitness,
                                                        experiment.aa_homogamy,
                                                        experiment.deafness_freq,
                                                        experiment.generations)['row']

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
        simulations = 0
        while simulations < experiment.simulations:
            start_time = time.time()
            p = [pool.apply_async(_worker) for i in range(mp_chunk_size)]
            table = [item.get() for item in p]
            experiment.write(table)
            simulations += mp_chunk_size
            rate = mp_chunk_size/(time.time()-start_time)
            time_remaining = (experiment.simulations-simulations)/rate if rate > 0 else 0
            print '{:,} simulations completed ' \
                  '({:,.0f}/min) '\
                  '{} remaining.'\
                  ''.format(simulations, 60*rate, _format_time(time_remaining))
            # mp_chunk_size is dynamically adjusted based on actual
            # execution speed such that file writes occur once per minute.
            mp_chunk_size = max(int(300*rate - 300*rate%cpu_count), cpu_count)
            if simulations + mp_chunk_size > experiment.simulations:
                mp_chunk_size = experiment.simulations-simulations
    print 'Saved to {}.'.format(experiment.filename)
    print 'Done.'
