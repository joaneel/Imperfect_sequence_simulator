##############################################################
# This shell masterscript creates datasets following a chosen 
# demographic model under the coalescent. it then simulates
# the sequencing process to create "imperfect" datasets.
# Standard summary statistics are then calculated.
#
# Elleouet joane.elleouet@alumni.ubc.ca
# Mon Dec 11 2017
##############################################################

#!/bin/bash

INDIR=/home/user/

#command arguments
MOD=$1 #model number
NIND=$2 #number of diploid genotypes
LOCNUM=$3 #number of sequences
LOCLENGTH=$4 #length of sequences
INDP=$5 #mean individual sequencing depth
SEQERR=$6 #pre-nucleotide sequencing error rate
NPODS=$7 #number of datasets to create

OUTDIR=$INDIR\POD_$MOD\_$NIND\_$LOCNUM\_$LOCLENGTH\_$INDP\_$SEQERR\/ #output directory

#create output directory
if [ ! -d "$OUTDIR" ]
then
mkdir $OUTDIR
fi

date > $OUTDIR\runtime_2.txt

#initialize nosegsite template: for compatibility with ms2fasta python script
NCHROM=$(($NIND*2))
echo "positions: 1" > nosite_template.txt
printf '0\n%.0s' `seq 1 $NCHROM` >> nosite_template.txt

for i in `seq 1 $NPODS`
do
        #create coalescent simulation
        #############################
        #initialize output file
        echo '#ms file' > $OUTDIR\POD$i\_allmarkers.ms
        #simulate sequences with scrm: creates as many files as number of markers
        Rscript --vanilla /data/user/scripts/ABC_POD_simulation.r $MOD $NIND $LOCNUM $LOCLENGTH $OUTDIR $i
        
        for j in `seq 1 $LOCNUM`
        do
                # modify scrm output to make invariant sequences compatible with ms2fasta
                ######################################################################
                sed -i -e '/^segsites: 0/r nosite_template.txt' $OUTDIR\scrm$i\.out_$j
                sed -i s/'segsites: 0'/'segsites: 1'/g $OUTDIR\scrm$i\.out_$j

                #Python script 1 to transform ms file created by scrm into a fasta file
                ########################################################################
                /data/user/programs/ms2fasta $OUTDIR\scrm$i\.out_$j > $OUTDIR\ms2fasta$i\_$j.out
                #rm $OUTDIR\scrm.out_$i

                #R script 2 to sequence and call genotypes
                ##########################################
                Rscript --vanilla /data/user/scripts/ABC_process_fasta_realistic_POD.r $NIND 1 $LOCLENGTH $INDP $SEQERR $OUTDIR $j $i
                
                #transform fastas of individual markers into ms output
                ######################################################

                sed -i s/a/A/g $OUTDIR\POD$i\_all_$j.fasta
                sed -i s/c/C/g $OUTDIR\POD$i\_all_$j.fasta
                sed -i s/g/G/g $OUTDIR\POD$i\_all_$j.fasta
                sed -i s/t/T/g $OUTDIR\POD$i\_all_$j.fasta
                sed -i s/n/N/g $OUTDIR\POD$i\_all_$j.fasta

                /data/user/programs/fastaconvtr-master/bin/fastaconvtr -F fasta -f ms -p 1 -i $OUTDIR\POD$i\_all_$j.fasta -o $OUTDIR\POD$i\_all_$j.ms  #ingores Ns (discards site)

                #concatenate  ms files
                grep -v '^#' $OUTDIR\POD$i\_all_$j.ms >> $OUTDIR\POD$i\_allmarkers.ms
                rm $OUTDIR\scrm$i.out_$j
                rm $OUTDIR\ms2fasta$i\_$j.out
                rm $OUTDIR\POD$i\_all_$j.fasta
                rm $OUTDIR\POD$i\_all_$j.ms_*
                rm $OUTDIR\POD$i\_all_$j.ms.log
                rm $OUTDIR\POD$i\_all_$j.ms
        done #end of j
        
                #calculate summary statistics
        /data/user/programs/msABC20120315/msABC $(($NIND*2)) $LOCNUM -I 2 $NIND $NIND --obs $OUTDIR\POD$i\_allmarkers.ms --options $INDIR\ms_ssdefs.txt > $OUTDIR\POD$i\_allmarkers.ms.mean

        #R script 3: summarize summary statistics across loci
        ######################################################
        Rscript --vanilla /data/user/scripts/ABC_summarize_POD.r $OUTDIR $i
done #end of i

#Rscript 4: put together params and stats for PODs
####################################################
Rscript --vanilla /data/user/scripts/ABC_create_PODs_results_table.r $OUTDIR $NPODS $MOD $LOCNUM

date >> $OUTDIR\runtime_2.txt

          
