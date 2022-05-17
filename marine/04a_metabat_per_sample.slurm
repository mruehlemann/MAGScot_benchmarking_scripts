#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="04a_metabat"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

###########################
####    METABAT     #######
###########################
source activate metagenome_env

s=$SLURM_ARRAY_TASK_ID

cd $TMPDIR

zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta

metabat2 -i sample.${s}.fasta -a $workfolder/$dataset/sample_${s}/sample_${s}.depth.txt -o sample_${s}_metabat/sample_${s}.metabat2 -t ${SLURM_CPUS_PER_TASK} -m 2000
for fa in $(ls sample_${s}_metabat/sample_${s}.metabat2.*.fa); do mv $fa ${fa}sta; done

grep '>' sample_${s}_metabat/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/$dataset/sample_${s}/sample_${s}.metabat2.contigs_to_bin.tsv
