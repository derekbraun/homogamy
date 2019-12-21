*Derek C. Braun, Brian H. Greenwald, Samir Jain, Eric Epstein, Brienna Herold, Maggie Gray
(*derek.braun@gallaudet.edu)

To clone this project: git clone https://github.com/derekbraun/homogamy.git

MODULES

fileio.py       Contains routines to standardize file I/O and to obtain data
                from files. Routines allow the user to select columns and
                keywords from the resulting data files.


simulator.py    Simulation module which uses simuPOP.
                Simulation parameters are set via setting variables in the
                experiment class.

                simulator.py --help to get parameters
                simulator.py to perform a single run


grapher.py      Produces graphs for the data files created by simulator.py.

                grapher.py --help to get parameters


stats.py        Performs statistical analyses comparing data files created by
                simulator.py.

                stats.py --help to get parameters.

my_math.py      Calculates ending values using math equations from Crow & Felsenstein (1968).
                Used to validate/compare with simulation results.


simulator.bash  Runs the simulations for this publication (takes forever!).

grapher.bash    Produces graphs for this publication.



USEFUL REFERENCES FOR simuPOP OBJECTS

This is a list of simuPOP objects and methods that are particularly important
for troubleshooting. During simulations, simuPOP keeps variables
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
stores them in the local namespace.


TROUBLESHOOTING simuPOP

sim.dump(pop)                   Outputs all data about all individuals
                                in a semi-tabular format.
describeEvolProcess(param)      Accepts the same parameters as Simulator.evolve.
                                Outputs a description of how an evolutionary
                                process will be executed.
moduleInfo                      Outputs a dict with information about the
                                version and build of simuPOP. See:
                                http://simupop.sourceforge.net/manual_svn/build/refManual_ch2_sec6.html#function-describeevolprocess



RUNNING a sample simuPOP from the Python interactive

import simuPOP as sim
pop = sim.Population(10000, loci=[1])
