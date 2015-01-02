INSTALLATION:

1.  Install this code from github.

2.  Install the major dependencies (python, numpy, matplotlib). 
    Directions for OSX are on my blog: https://derekbraun.wordpress.com/
    
3.  Install simupop


TO-DO LIST
        
1.  Define an Experiment class that manages experimental data and metadata.
    Give the class methods for reading and writing data.
    This should make code easier to follow.
    
2.  simuPOP obviously does not work. Why not? Look through Bo Peng's book
    before consulting with Bo himself.

3.  Change plotting (and data collection) to AA_freq instead of AA_size.
    (I will need to go through simuPOP examples to make this work)
    (I forgot why I even cared in the first place, because freqs can be
    calculated rapidly from sizes.)
    
4.  Calculate the inbreeding coefficient, F, after each generation. We
    might want to use the simupop object as an iterable to accomplish
    this instead of using multiple PyEval statements.
    There should be examples on how to do this in Bo Peng's book.

5.  Generate and output a final table that gives final medians and HPDs
    for AA_freq, Aa_freq, aa_freq, A_freq, a_freq, F.

6.  Plot the inbreeding coefficient, F, calculated at each time point,
    on top of each of the figures.