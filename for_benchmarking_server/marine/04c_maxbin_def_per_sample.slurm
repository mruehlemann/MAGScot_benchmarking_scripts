#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="04c_maxbin"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/biscoreto_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

#############################
####    MAXBIN DEFAULT #######
#############################
source activate metagenome_env
module load gcc

cd $TMPDIR

s=$SLURM_ARRAY_TASK_ID

zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta

awk '{if(NR>1) print $1"\t"$3}' $workfolder/$dataset/sample_${s}/sample_${s}.depth.txt | grep "^S${s}C" > sample_${s}.abu

mkdir -p sample_${s}_maxbin2_out
perl /work_ifs/sukmb276/software/MaxBin-2.2.4/run_MaxBin_prodigal.pl -contig sample.${s}.fasta -prodigal $workfolder/$dataset/sample_${s}/sample_${s}.prodigal.faa -out sample_${s}_maxbin2_out/sample_${s}.maxbin2 -abund sample_${s}.abu -thread ${SLURM_CPUS_PER_TASK} -min_contig_length 2000

grep '>' sample_${s}_maxbin2_out/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/$dataset/sample_${s}/sample_${s}.maxbin2.contigs_to_bin.tsv
rm -r binning/maxbin_multi_out/
