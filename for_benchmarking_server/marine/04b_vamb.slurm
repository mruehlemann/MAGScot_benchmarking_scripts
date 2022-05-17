#!/bin/bash
#SBATCH --reservation=sukmb412_13
#SBATCH --gres=gpu:1
#SBACTH --nodelist=medg01
#SBATCH -c 4
#SBATCH --mem=100gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="04b_vamb"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/biscoreto_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

###########################
####    VAMB        #######
###########################

module load miniconda3
source activate vamb_gpu_env

vamb --cuda -p ${SLURM_CPUS_PER_TASK} --outdir $workfolder/$dataset/VAMB_out --fasta $workfolder/${dataset}/${dataset}.catalogue.fna.gz --jgi $workfolder/${dataset}/${dataset}.all_depth.tsv -o C
