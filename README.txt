*Derek C. Braun, Brian H. Greenwald, Samir Jain, Eric Epstein, Brienna Herold, Maggie Gray
(*derek.braun@gallaudet.edu)

To clone this project: git clone https://github.com/derekbraun/homogamy.git

DATASET

*.tsv           Contain data from various experiments used in our publication.


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
                  stats.py --help to get parameters

my_math.py      Calculates ending values using math equations from Crow & Felsenstein (1968).
                Used to validate/compare with simulation results and also
                cited in our publication.

simulator.bash  Runs the simulations for this publication (warning: slow!).

grapher.bash    Produces the graphs used in our publication.
