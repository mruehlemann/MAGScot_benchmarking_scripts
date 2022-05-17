#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=100gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="mb_01_catalog"

##########################
####    SETUP     #######
###########################

echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/$dataset/rawdata"


##################################
####    start collection   ########
##################################

cd $workfolder/subgroups/${this_local}
Rscript /work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/scripts/MAGScoT_benchmark_scripts/hmp2/02b_collect_mapping.R
sed -i "s/[.]var/-var/g" ${dataset}.all_depth.tsv
