homogamy

Samir Jain, Eric Epstein, Trevor Klemp, Maggie Gray, Selman Jawed, Derek 
Braun* (*derek.braun@gallaudet.edu)
    
Contains routines to standardize file I/O and to select columns and
keywords from the resulting data files. Greatly simplifies coding
elsewhere.

To clone this project: git clone https://github.com/derekbraun/homogamy.git


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



TO-DO LIST

1.              Fix deafChooser in simulator.py to incorporate fitness,
                defined as the relative number of children (major)
                
2.              Fix grapher.py to incorporate both rc and rcfiles (easy)

3.              Fix grapher.py to incorporate scaling (e.g. per 100,000)
                Also fix label formatting.

4.              Adjust grapher.py output (label sizes, axis lines, etc).


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
pop.individuals                             an iterable to iterate through each
                                            individual in a population. Useful
                                            for either getting or setting
                                            attributes for each individual.
                                            e.g. 
                                            "for individual in pop.individuals"
pop.individual(n).sex()                     either displays 1 or 2 to 
                                            represent M/F
pop.individual(n).genotype()                displays the genotype.


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