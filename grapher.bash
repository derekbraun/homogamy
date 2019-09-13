#!/bin/bash

# Figure 1: Effect over time of homogamy on the frequencies of genetically
# deaf individuals and a recessive deafness allele
./grapher.py -f "aa" pop200k_hom0.9_fit1.0.tsv -t ""
./grapher.py -f "a" pop200k_hom0.9_fit1.0.tsv -t ""

# Figure 2: Effect of homogamy on the frequencies of genetically deaf
# individuals and a recessive deafness allele

./grapher.py -f "aa" pop200k_hom0.?_fit1.0.tsv --ylim 0.026
./grapher.py -f "a" pop200k_hom0.?_fit1.0.tsv --ylim 0.026

# Figure 3: Synergy of homogamy and relative fitness on the frequencies of
# genetically deaf individuals and a recessive deafness allele

./grapher.py -f "aa" pop200k_hom0.9_fit?.?.tsv --y_format "{x:.1%}"
./grapher.py -f "a" pop200k_hom0.9_fit?.?.tsv --y_format "{x:.1%}"
./grapher.py -f "aa" pop200k_hom0.0_fit?.?.tsv --y_format "{x:.1%}"
./grapher.py -f "a" pop200k_hom0.0_fit?.?.tsv --y_format "{x:.1%}"

mv *.pdf ~/Desktop
