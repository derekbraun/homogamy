homogamy

Samir Jain, Eric Epstein, Maggie Gray, Derek Braun*
(*derek.braun@gallaudet.edu)

Contains routines to standardize file I/O and to select columns and
keywords from the resulting data files. Greatly simplifies coding
elsewhere.

To clone this project: git clone https://github.com/derekbraun/homogamy.git

NANCE AND Kearsey
First 5 generations: fitness 0
Homogamy increases to 90% by gen 5
20 more generations (400 years) representing 1800-2200
Population size fixed at 200,000


TO-DO LIST

1.  Fix variable names


1.  Final simulations parameters:
    x 0% homogamy, normal fitness, adventitious deafness
    x 90% homogamy, normal fitness, adventitious deafness
    x 0% homogamy, 2x fitness, adventitious deafness
    x 90% homogamy, 2x fitness, adventitious deafness
    0 % homogamy, normal fitness, no adventitious deafness
    90 % homogamy, normal fitness, no adventitious deafness
    0 % homogamy, normal fitness, adventitious deafness, small population (50k?)
    90% homogamy, normal fitness, adventitious deafness, small population (50k?)

5.  Fix grapher.py to make more beautiful graphs:
    a. incorporate scaling (e.g. per 100,000)
    b. more lines than just 5% and 95%
    c. fix label formatting..

6.  Maybe: Base population size on actual historical population growth
    Determine whether population model has any real effect on frequencies.
    (report this in the manuscript)




MODULES

fileio.py       Contains routines to standardize file I/O and to obtain data
                from files. Routines allow the user to select columns and
                keywords from the resulting data files.


simulator.py    Simulation module which uses simuPOP.
                Simulation parameters are set via setting global variables
                and, for changing assortative mating, sometimes by changing
                the code.

                simulator.py --help to get parameters
                simulator.py --test-run to perform a single run


grapher.py      Produces graphs for the data files created by Simulator.py.

                grapher.py --help to get parameters


stats.py        Performs statistical analyses comparing data files created by
                simulator.py.

                stats.py --help to get parameters.




USEFUL REFERENCES FOR simuPOP OBJECTS

This is a list of simuPOP objects and methods that are particularly important
for troubleshooting. SimuPOP's variable interface is murky at best and really
does not follow Python conventions. During simulations, simuPOP keeps variables
in a dict called Population.vars() or as attributes of an object
called Population.dvars()

pop.genotype(0)                             lists all genotypes.
                                            Doesn't work for vsps.
pop.dvars(sp).alleleFreq[locus][allele]     the allele frequency at
                                            locus [locus] and allele [allele]
                                            for subpopulation sp. Note that
                                            these are dicts, not lists.
pop.individuals()                           an iterable to iterate through each
                                            individual in a population. Useful
                                            for either getting or setting
                                            attributes for each individual.
                                            e.g.
                                            "for individual in pop.individuals"
pop.individual(n).sex()                     either displays 1 or 2 to
                                            represent M/F
pop.individual(n).genotype()                displays the genotype.
pop.individual(n).affected()                displays whether affected
pop.individual(n).setAffected(True)         sets the individual to Affected
                                            without changing genotype

Population.evolve() has an "operator" called Stat which calculates things and
stores them in the local namespace under an unspecified name (argh!!!).
But this name could probably be figured out by inspecting Population.vars().


TROUBLESHOOTING simuPOP

sim.dump(pop)                   Outputs all data about all individuals
                                in a semi-tabular format.
describeEvolProcess(param)      Accepts the same parameters as Simulator.evolve.
                                Outputs a description of how an evolutionary
                                process will be executed.
moduleInfo                      Outputs a dict with information about the
                                version and build of simuPOP. See:
                                http://simupop.sourceforge.net/manual_svn/build/refManual_ch2_sec6.html#function-describeevolprocess



RUNNING simuPOP from the Python interactive

import simuPOP as sim
pop = sim.Population(10000, loci=[1])
