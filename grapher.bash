#!/bin/bash
#
# Figure 1: Effect over time of homogamy on the frequencies of genetically
# deaf individuals and a recessive deafness allele
./grapher.py -f "aa" pop200k_hom0.9_fit1.0.tsv -t "" -o "Fig_1a.png"
./grapher.py -f "a" pop200k_hom0.9_fit1.0.tsv -t "" --ylim 0.025 -o "Fig_1b.png"

# Figure 2: Effect of homogamy on the frequencies of genetically deaf
# individuals and a recessive deafness allele

./grapher.py -f "aa" pop200k_hom0.?_fit1.0.tsv --ylim 0.0005 --ylabel "Deaf individuals after 20 generations" -o "Fig_2a.png"
./grapher.py -f "a" pop200k_hom0.?_fit1.0.tsv --ylim 0.025 --ylabel "Allelic frequency after 20 generations" -o "Fig_2b.png"

# Figure 3: Synergy of homogamy and relative fitness on the frequencies of
# genetically deaf individuals and a recessive deafness allele

./grapher.py -f "a" pop200k_hom0.0_fit?.?.tsv --ylim 1.0 --ylabel "Allelic frequency after 20 generations" --yformat "{x:.1%}" --yscale "log" -o "Fig_3a.png"
./grapher.py -f "a" pop200k_hom0.9_fit?.?.tsv --ylim 1.0 --ylabel "Allelic frequency after 20 generations" --yformat "{x:.1%}" --yscale "log" -o "Fig_3b.png"

mv *.png ~/Desktop
