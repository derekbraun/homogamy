#!/bin/bash
#
EXT="eps"  # encapsulated postscript file

# Figure 1: Effect over time of homogamy on the frequencies of genetically
# deaf individuals and a recessive deafness allele

./grapher.py pop200k_hom0.9_fit1.0.tsv \
             --field "aa" \
             --title "" \
             --output_file "Fig_1a.$EXT"

./grapher.py pop200k_hom0.9_fit1.0.tsv \
             --field "a" \
             --title "" \
             --ylim 0.025 \
             --output_file "Fig_1b.$EXT"

# Figure 2: Effect of homogamy on the frequencies of genetically deaf
# individuals and a recessive deafness allele

./grapher.py pop200k_hom0.?_fit1.0.tsv \
             --field "aa" \
             --ylim 0.0005 \
             --ylabel "Deaf individuals after 20 generations" \
             --output_file "Fig_2a.$EXT"

./grapher.py pop200k_hom0.?_fit1.0.tsv \
             --field "a" \
             --ylim 0.025 \
             --ylabel "Allelic frequency after 20 generations" \
             --output_file "Fig_2b.$EXT"

# Figure 3: Synergy of homogamy and relative fitness on the frequencies of
# genetically deaf individuals and a recessive deafness allele

./grapher.py pop200k_hom0.0_fit?.?.tsv \
            --field "a" \
            --ylim 1.0 \
            --ylabel "Allelic frequency after 20 generations" \
            --yformat "{x:.1%}" \
            --yscale "log" \
            --output_file "Fig_3a.$EXT"

./grapher.py  pop200k_hom0.9_fit?.?.tsv \
            --field "a" \
            --ylim 1.0 \
            --ylabel "Allelic frequency after 20 generations" \
            --yformat "{x:.1%}" \
            --yscale "log" \
            --output_file "Fig_3b.$EXT"

mv *.$EXT ~/Desktop
