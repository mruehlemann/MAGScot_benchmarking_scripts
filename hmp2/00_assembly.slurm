#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="00_assembly"

###########################
####    SETUP     #######
###########################

echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/$dataset/rawdata"



####################################
####    assembly    ########
####################################
source activate metagenome_env

samples=($(ls $datafolder/ | grep tar | cut -d '.' -f 1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}

cd $workfolder/$dataset/
mkdir -p $sample

cd $TMPDIR

megahit -1 $datafolder/${sample}_R1.fastq.gz -2 $datafolder/${sample}_R2.fastq.gz -o megahit -t ${SLURM_CPUS_PER_TASK}

cat megahit/final.contigs.fa | sed "s/k119/$sample\_contig/" > $workfolder/$dataset/$sample/$sample.megahit.fasta
