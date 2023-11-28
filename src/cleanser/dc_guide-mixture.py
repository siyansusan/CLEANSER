#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2018 William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import math
import ProgramName
from Rex import Rex
rex=Rex()
import TempFilename
from StanParser import StanParser
from Stan import Stan
from DataFrame import DataFrame
from SummaryStats import SummaryStats
import getopt
import random

DEBUG=True
STDERR=TempFilename.generate(".stderr")
INPUT_FILE=TempFilename.generate(".staninputs")
INIT_FILE=TempFilename.generate(".staninit")
OUTPUT_TEMP=TempFilename.generate(".stanoutputs")

#def writeInitializationFile(filename):
#    OUT=open(filename,"wt")
#    r=random.uniform(0.01,0.99)
#    mu=random.gauss(0,1)
#    sigma=math.exp(random.gauss(0,2))
#    print("r <-",r,file=OUT)
#    print("mu <-",mu,file=OUT)
#    print("sigma <-",sigma,file=OUT)
#    OUT.close()

def writeInputsFile(stan,X,filename,numZeros,L):
    OUT=open(filename,"wt")
    print("N <-",str(len(X)),file=OUT)
    print("numZeros <-",numZeros,file=OUT)
    stan.writeOneDimArray("X",X,len(X),OUT)
    stan.writeOneDimArray("L",L,len(L),OUT)
    #print(str(len(X)))
    #print(str(len(L)))
    OUT.close()

def runSTAN(model,X,numWarmup,numSamples,infile,outfile,numZeros,L):
    # Create STAN object
    stan=Stan(model)

    # Write input and initialization files
    writeInputsFile(stan,X,INPUT_FILE,numZeros,L)
    global INIT_FILE
    #writeInitializationFile(INIT_FILE)

    # Run STAN model
    #INIT_FILE=None; ### DEBUGGING
    if(DEBUG):
        print(stan.getCmd(numWarmup,numSamples,INPUT_FILE,OUTPUT_TEMP,STDERR,INIT_FILE))
    stan.run(numWarmup,numSamples,INPUT_FILE,OUTPUT_TEMP,STDERR,INIT_FILE)

    # Parse MCMC output
    parser=StanParser(OUTPUT_TEMP)
    r=parser.getVariable("r")
    n_nbMean = parser.getVariable("n_nbMean")
    n_nbDisp = parser.getVariable("n_nbDisp")
    #LAMBDA = parser.getVariable("lambda")
    nbMean=parser.getVariable("nbMean")
    nbDisp = parser.getVariable("nbDisp")

    #sLAMBDA = parser.getVariable("slambda")
    #snbMean = parser.getVariable("snbMean")
    (med_n_nbMean,a,b) = parser.getMedianAndCI(0.95,"n_nbMean")
    (med_nbMean,a,b)=parser.getMedianAndCI(0.95,"nbMean")
    (assignments) =getAssignments(parser)
    return (r,nbMean, nbDisp, n_nbMean, n_nbDisp, med_n_nbMean, med_nbMean, assignments)

#pass lambda and nbmean to compare, lambda should be smaller than nbmean 
def getAssignments(parser):
    assignments=[]
    likeli_n_negbin=[]
    likeli_negbin=[]
    n_negbin_sum=[]
    negbin_sum=[]
    names=parser.getVarNames()
    for name in names:
        #print(name)
    #do something similar
        if(len(name)>=4 and name[:4]=="PZi."):
            samples=parser.getVariable(name)
            median=SummaryStats.median(samples)
            assignments.append(median)
        #if len(name)>=4 and name[:10]=="likeli_neg":
        #    samples=parser.getVariable(name)
        #    median=SummaryStats.median(samples)
        #    likeli_negbin.append(median)
        #if(len(name)>=4 and name[:10]=="likeli_poi"):
        #    samples=parser.getVariable(name)
        #    median=SummaryStats.median(samples)
        #    likeli_poisson.append(median)
        #if (name[:9]=="lambdasum"):
        #    samples=parser.getVariable(name)
        #    median=SummaryStats.median(samples)
        #    poisson_sum.append(median)
        #if (name[:5]=="nbsum"):
        #    samples=parser.getVariable(name)
        #    median=SummaryStats.median(samples)
        #    negbin_sum.append(median)


    return (assignments)

def writeAssignments(cellIDs,assignments,assignFile, n_mu, mu):
    with open(assignFile,"wt") as ASSIGN:
        N=len(assignments)
        z = 0
        if float(n_mu) >= float(mu):
            z = 1
        for i in range(N):
            #print(str(i))
            cellID=cellIDs[i]
            x=assignments[i]
            #poisson=likeli_poisson[i]
            #negbin= likeli_negbin[i]
            #lam = LAMBDA[i]
            #Mu = mu[i]
            #if float(lam) >= float(Mu):
            #    x = 0
            #if z == 1:
            #    x = 0
            print(cellID,x, sep="\t",file=ASSIGN)

def writeSamples(r,mu,disp,n_mu, n_disp, outfile):
    N=len(r)
    with open(outfile,"wt") as SAMPLES:
        print("r\tmu\tDisp\tn_mu\tn_disp",file=SAMPLES)
        for i in range(N):
            print(r[i],mu[i], disp[i], n_mu[i], n_disp[i], sep="\t",file=SAMPLES)

#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:")
if(len(args)!=8):
    #print(str(len(args)))
    exit(ProgramName.get()+" [-s stanfile] <model> <input.txt> <r-and-mu.txt> <Zi.txt> <#warmup> <#keep> <#zeros> <cell_lib>\n   -s = save raw STAN file\n")
(model,inFile,outfile,assignFile,numWarmup,numSamples,numZeros, lib)=args
stanFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
numZeros=int(numZeros)

cell_dict = {}

# Read inputs
#df=DataFrame.readTable(inFile,header=False,rowNames=False)
#df.toInt()
#cellID=df.getColI(0).getRaw()
#X=df.getColI(1).getRaw() # read counts
cellID = []
X = []
L = []

with open(lib, "r") as lib_file:
    for line in lib_file:
        (cell_ID, lib_size) = line.strip().split("\t")
        cell_dict[cell_ID] = lib_size

with open(inFile, "r") as count_file:
    for line in count_file:
        (cell_ID, guide_count) = line.strip().split("\t")
        L.append(cell_dict[cell_ID])
        X.append(guide_count)
        cellID.append(cell_ID)
#print(L[:5])
#print(X[:5])
#print(cellID[:5])
# Run STAN
(r,mu,disp,n_mu, n_disp, med_n_mu, med_mu,  assignments)=runSTAN(model,X,numWarmup,numSamples,INPUT_FILE,INIT_FILE,numZeros,L)

#print(str(len(LAMBDA)), str(len(assignments)))

print("r=",SummaryStats.median(r),
      "\tmu=",SummaryStats.median(mu),"\tdisp=",SummaryStats.median(disp),"\tn_mu=",SummaryStats.median(n_mu),"\tn_disp=",SummaryStats.median(n_disp),sep="")
#print("poisson_sum=", poisson_sum, "\tnb_sum=", nb_sum, sep="")

# Write samples into output file
writeSamples(r,mu,disp, n_mu, n_disp, outfile)
#print(str(len(LAMBDA)))
#print(str(len(mu)))
writeAssignments(cellID,assignments,assignFile, med_n_mu, med_mu)

# Clean up
if(not DEBUG):
    os.remove(STDERR)
    os.remove(INPUT_FILE)
if(stanFile is None): 
    if(not DEBUG): os.remove(OUTPUT_TEMP)
else: os.system("mv "+OUTPUT_TEMP+" "+stanFile)

