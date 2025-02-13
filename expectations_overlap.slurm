#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=expct_over # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. Increase memory for this sofware because does not handle multiple cpus
#SBATCH --output=expct_over.txt  #output


WD=/mnt/SCRATCH/cedi/phDSalmon/evolutionary_analysis/pauline_data #shortcut for the working directory

echo "Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job."

NbWindTaj=1053
NbWindFST=1053
NbSNPxpehh=1982
SizeWindTaj=20000
SizeWindFST=20000
SizeWindxpehh=1
SizeGenome=2499322922


### make 1000 random expectations
echo -e "Run\toverlap_FST_TAJ\toverlap_XPEHH_TAJ\toverlap_FST_XPEHH">>Overlap_expected.txt

for ((i=1; i<=1000; i++))
do
    ### make expectations for Tajima

    rm RandWindTaj.bed

    for ((j=1; j<="$NbWindTaj"; j++))
        do
            while true; do
                StartRanWind=$((1 + $(od -An -N4 -tu4 < /dev/urandom) % SizeGenome))

                EndRanWind=$((StartRanWind + SizeWindTaj))

                if (( EndRanWind <= SizeGenome )); then
                    break  
                fi
            done

            echo -e "1\t$StartRanWind\t$EndRanWind">>RandWindTaj.bed

        done


    ### make expectations for FST

    rm RandWindFST.bed

    for ((f=1; f<="$NbWindFST"; f++))
        do
            while true; do
                StartRanWind=$((1 + $(od -An -N4 -tu4 < /dev/urandom) % SizeGenome))

                EndRanWind=$((StartRanWind + SizeWindFST))

                if (( EndRanWind <= SizeGenome )); then
                    break  
                fi
            done

            echo -e "1\t$StartRanWind\t$EndRanWind">>RandWindFST.bed

        done


    ### make expectations for XPEHH

    rm RandWindxpehh.bed

    for ((f=1; f<="$NbSNPxpehh"; f++))
        do

            StartRanWind=$((1 + $(od -An -N4 -tu4 < /dev/urandom) % SizeGenome))


            echo -e "1\t$StartRanWind\t$StartRanWind">>RandWindxpehh.bed

        done

    ## make overlap

    singularity exec containers/bedtools_2_30_hc088bd4_0.sif bedtools intersect -a RandWindTaj.bed -b RandWindFST.bed -wa -wb > temp_overlap_FST_TAJ.txt

    Nboverlap_FST_TAJ=$(cat temp_overlap_FST_TAJ.txt | wc -l)

    singularity exec containers/bedtools_2_30_hc088bd4_0.sif bedtools intersect -a RandWindTaj.bed -b RandWindxpehh.bed -wa -wb > temp_overlap_XPEHH_TAJ.txt

    Nboverlap_XPEHH_TAJ=$(cat temp_overlap_XPEHH_TAJ.txt | wc -l)

    singularity exec containers/bedtools_2_30_hc088bd4_0.sif bedtools intersect -a RandWindxpehh.bed -b RandWindFST.bed -wa -wb > temp_overlap_FST_XPEHH.txt

    Nboverlap_FST_XPEHH=$(cat temp_overlap_FST_XPEHH.txt | wc -l)

    echo -e "$i\t$Nboverlap_FST_TAJ\t$Nboverlap_XPEHH_TAJ\t$Nboverlap_FST_XPEHH">>Overlap_expected.txt

done