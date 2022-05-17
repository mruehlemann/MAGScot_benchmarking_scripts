#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=100gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="mb_01_catalog"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

####################################
####    CATALOG CREATION    ########
####################################

cd $datafolder/$dataset
for i in $(seq 0 9); do curl -O https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_sample_${i}_contigs.tar.gz; done
for i in $(seq 0 9); do curl -O https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_sample_${i}_reads.tar.gz; done

for i in $(seq 1 9); do tar -xvzf marmgCAMI2_sample_${i}_contigs.tar.gz; done
for i in $(seq 1 9); do tar -xvzf marmgCAMI2_sample_${i}_reads.tar.gz; done
