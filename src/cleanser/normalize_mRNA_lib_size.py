#!/usr/bin/env python
import sys
import statistics
import os
import ProgramName
import gzip
import numpy as np
#=========================================================================
# main()
#=========================================================================

if len(sys.argv)!=2:
    exit(ProgramName.get()+"<mRNA_lib.txt>\n")
(lib_size)=sys.argv[1]

norm_dict = {}
count = 0
tot_size = 0

with open(lib_size, "r") as lib_file:
    for line in lib_file:
        (cell_num, size) = line.strip().split("\t")
        count += 1
        tot_size += int(size)
        norm_dict[cell_num] = size

avg_size = int(tot_size)/int(count)

for key_cell in norm_dict:
    lib_size = int(norm_dict[key_cell])
    norm_size = lib_size/avg_size
    size = str(norm_size)
    print(key_cell, size, sep = "\t")


