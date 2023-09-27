#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
from MatrixMarket import MatrixMarket
from SlurmWriter import SlurmWriter
from Rex import Rex
rex=Rex()

JOB_NAME="GUIDES"
MAX_PARALLEL=300
NUM_WARMUP=300
NUM_KEEP=1000

def getGuideID(group):
    if(len(group)==0): raise Exception("no lines in group")
    line=group[0]
    guideID=line[0]
    return guideID

def writeCounts(group,countsFile):
    with open(countsFile,"wt") as OUT:
        for line in group:
            cellID=line[1]
            count=line[2]
            print(cellID,count,sep="\t",file=OUT)

def writeSlurm(countsFile,slurmFile,numZeros,dataDir,fileIndex):
    samplesFile=dataDir+"/samples/"+str(fileIndex)+".txt"
    posteriorsFile=dataDir+"/posteriors/"+str(fileIndex)+".txt"
    cmd="module load python/3.7.4-gcb01 ; cd /gpfs/fs1/data/gersbachlab/susan/modeling/guide-mixture/mia/dcas9_barnyard/dcas9_dc_72hr/; python3 guide-mixture.py guide-mixture "+countsFile+" "+samplesFile+" "+posteriorsFile+" "+str(NUM_WARMUP)+" "+str(NUM_KEEP)+" "+str(numZeros)+" /gpfs/fs1/data/gersbachlab/susan/modeling/guide-mixture/mia/dcas9_barnyard/dcas9_dc_72hr/normalized_cell_lib_guide.txt"
    slurm.addCommand(cmd)

def getTotalCells(M):
    header=M.getHeader()
    if(len(header)!=3): raise Exception("bad header")
    totalCells=header[1]
    return totalCells

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <guides.mtg.gz> <out-slurm-dir> <out-data-dir>\n")
(guideFile,slurmDir,dataDir)=sys.argv[1:]
if(rex.find("(.*)/$",dataDir)): dataDir=rex[1]
if(rex.find("(.*)/$",slurmDir)): slurmDir=rex[1]

slurm=SlurmWriter()
M=MatrixMarket(guideFile)
totalCells=None
while(True):
    group=M.nextGroup(0)
    if(group is None): break
    guideID=getGuideID(group)
    if(totalCells is None): totalCells=getTotalCells(M)
    countsFile=dataDir+"/counts/"+str(guideID)+".txt"
    writeCounts(group,countsFile)
    slurmFile=slurmDir+"/slurm-"+str(guideID)+".slurm"
    numZeros=totalCells-len(group)
    writeSlurm(countsFile,slurmFile,numZeros,dataDir,guideID)
    
slurm.mem(5000)
slurm.setQueue("all")
slurm.writeArrayScript(slurmDir,JOB_NAME,MAX_PARALLEL);



