import os
import readline
import string
import sys
import mat
import random
from os import system
import re

## Set the key variables:
path_to_data = "Data/PGV.Bi/3D7.HB3/";
recomb_rate = 0.0001;
geno_rate   = 0.001;

## Compile the C++ code:
line = "g++ $(pkg-config --cflags gsl) -c hmm.cpp";
os.system(line);

line = "g++ hmm.o $(pkg-config --libs gsl) -o hmm"
os.system(line);

## Run the HMM to produce Viterbi (color) paths:
line = "./Code/hmm Data/PGV.Bi/3D7.HB3/ 0.00001 0.0001";
os.system(line);

## R script to create color figures of cross results

#>	parse out figure path

fig_path = "Figures/PGV.Bi/3D7.HB3";

#>	make figure path, if needed 

command = "mkdir "+str(fig_path);
print(command);
os.system(command);

#>	execute script

command = "cat test.r | R --no-save --args test test1 test2^C"
print(command);

## R script to create pipeline data files





