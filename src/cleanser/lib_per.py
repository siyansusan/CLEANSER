#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author:Susan Liu
#=========================================================================
import sys
import statistics
import os
import ProgramName
import gzip

#=========================================================================
# main()
#=========================================================================

if len(sys.argv)!=3:
    exit(ProgramName.get()+" <filename.mtx.gz> <col>\n")
mtxFile,col=sys.argv[1:]

col=int(col)
d={}

with gzip.open(mtxFile,"r") as matr:
    matr.readline()
    matr.readline()
    matr.readline()
    for line in matr:
        line = line.decode("utf8")
        (guide,cell,lib)=line.strip().split()
        if col==1:
            sort_by=cell
        if col==0:
            sort_by=guide
        if sort_by in d:
            d[sort_by]+=int(lib)
        else:
            d[sort_by]=int(lib)

#lib_size=list(d.values())
#med_lib_size=statistics.median(lib_size)

#for key in d:
#    d[key]=med_lib_size/d[key]

for key in d:
    print(str(key)+"\t"+str(d[key]))


#TO FIND TOP
#sorted_d=sorted(d,key=d.get)
#print(sorted_d)
#print(d["44118"])
