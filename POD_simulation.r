#############################################################################
# This R script uses a set of parameters (stored as a table
# in the file "trueparam.txt") and outputs simulated sequences (1 per file).
# Here, a 2-population expansion
# model is used as an example.
#
# Joane Elleouet joane.elleouet@alumni.ubc.ca
# December 2017
#############################################################################

#!/usr/bin/env Rscript

options(scipen=999) #disable scientific notation (python script ms2fasta does not like it)

args=commandArgs(trailingOnly=TRUE)

N0=10000 # this is the pop size of a hypothetical pop used to scale T

#dataset
mod=as.numeric(args[1])
sampsize=2*as.numeric(args[2])
nbloci=as.numeric(args[3])
loclength=as.numeric(args[4])
dir=args[5]
iter=as.numeric(args[6])

#load parameter file
#####################
truepar=read.table("/data/user/trueparam.txt",header=T)

# define parameters
####################

n1 = truepar$n1[iter]
n2 = truepar$n2[iter]
if (mod==2 | mod==4){n02 = truepar$n02[iter]} else {n02 = 2}
N1 = n1/N0
N2 = n2/N0
N02 = n02/N0
if (mod==3 | mod==4){m21 = truepar$m21[iter]} else {m21 = 0}
logTEXP=truepar$logTEXP[iter]
rec = 10^-8
mu = 9*10^-9

#COMPLEX PARAMETERS
TEXP = round(exp(logTEXP))
TCOAL = TEXP/(4*N0)
gr2 = -log(N02/N2)/TCOAL
thetamu = 4*N0*loclength*mu
R = 4*N0*(loclength-1)*rec
M21 = 4*N0*(loclength)*m21

#run scrm for each sequence
###########################
for (j in 1:nbloci) {
        system(paste("scrm",sampsize,1,"-t",thetamu,"-r",R,loclength,"-I",2,sampsize/2,sampsize/2,"-n",1,N1,"-n",2,N2,"-g",2,gr2,"-m",2,1,M21,"-ej",TCOAL,2,1,"-eg",TCOAL,2,0,"-SC abs", paste("> ",dir,"scrm",iter,".out_",j,sep=""),sep=" "))
}
