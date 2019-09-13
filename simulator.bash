#!/bin/bash

for homogamy in 0.0 0.3 0.6 0.9
   do
      for fitness in 0.0 1.0 #0.0 0.5 1.0 1.5 2.0
         do
            ./simulator.py --write --homogamy $homogamy --fitness $fitness
         done
   done
