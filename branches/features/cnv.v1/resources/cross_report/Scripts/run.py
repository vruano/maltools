import os
import readline
import string
import sys
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
#os.system(line);

## R script sequence to create color figures of cross results

#>	parse out figure path

splits = path_to_data.split("/");

fig_sub_path = "Figures/"+splits[1];
fig_path     = "Figures/"+splits[1]+ "/" + splits[2];

crosses = splits[2].split(".");
cross_1      = crosses[0];
cross_2      = crosses[1];

#>	make figure path, if needed 

command = "mkdir "+str(fig_sub_path);
print(command);
os.system(command);

command = "mkdir "+str(fig_path);
print(command);
os.system(command);

#>	execute script

command = "cat Scripts/calc.errors.r | R --no-save --args " + str(path_to_data) + " " + str(fig_path) + str("/") + " " + str(cross_1) + " " +str(cross_2);
print(command);
os.system(command)

## R script to create pipeline data files

#>	parse out output path

splits = path_to_data.split("/");

out_sub_path = "Output/"+splits[1];
out_path     = "Output/"+splits[1]+ "/" + splits[2];

#>	make figure path, if needed 

command = "mkdir "+str(out_sub_path);
print(command);
os.system(command);

command = "mkdir "+str(out_path);
print(command);
os.system(command);

#>	execute script

command = "cat Scripts/plot.cnc.uniq.r | R --no-save --args " + str(path_to_data) + " " + str(out_path) + str("/");
print(command)
os.system(command);





