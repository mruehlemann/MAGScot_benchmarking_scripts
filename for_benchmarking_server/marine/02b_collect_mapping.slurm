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

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/biscoreto_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

##################################
####    start collection   ########
##################################

cd $workfolder/$dataset/
Rscript /work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/biscoreto_benchmark/scripts/BiScoReTo_benchmark_scripts/02b_collect_mapping.R
sed -i "s/[.]var/-var/g" ${dataset}.all_depth.tsv
