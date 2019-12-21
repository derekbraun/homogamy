#!/usr/local/bin/python3 -u
# -*- coding: utf-8 -*-
# We generally follow PEP 8: http://legacy.python.org/dev/peps/pep-0008/

'''
*Derek C. Braun, Brian H. Greenwald, Samir Jain, Eric Epstein, Brienna Herold, Maggie Gray
(*derek.braun@gallaudet.edu)

Calculates ending values using math equations from Crow & Felsenstein (1968).
Used to validate/compare with simulation results.
'''

# defaults; some of these can be set otherwise by passing arguments
a_FREQ              = 0.01304  # Nance and Kearsey: 0.01304
aa_HOMOGAMY         = 0.9      # Nance and Kearsey: 0.9
DEAF_FREQ           = 0.0008   # Nance and Kearsey: two other genes each at 0.003^2 freq (tiny)
GENERATIONS         = 20       # Nance and Kearsey: 5 gen with 0 fitness + 20 gen with 1.0 fitness

import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--homogamy',
                        action = 'store',
                        type = float,
                        default = aa_HOMOGAMY,
                        help = 'deaf-deaf assortative mating (homogamy) ' \
                               '(default {}).'.format(aa_HOMOGAMY))
    args=parser.parse_args()

    headers = ['gen', 'aa']
    rows = []

    q = a_FREQ
    p = 1-q
    Rt = q**2
    print('{:^7}   {:^7}   {:^7}'.format('gen', 'r', 'aa'))
    print('{:^7}   {:^7}   {:^7}'.format('-'*7, '-'*7, '-'*7))
    for gen in range(GENERATIONS):
        r = args.homogamy * Rt/(DEAF_FREQ-q**2+Rt)  # recalculate r after each gen for more accuracy
        Rt = (1-r)*q**2 + r*(q**2+Rt*(p-q))/(1-Rt)  # equation 3 from Crow & Felsenstein (1968);
        print('{:^7}   {:^7.3%}   {:^7.4%}'.format(gen, r, Rt))
    print('Done.')
    exit()
