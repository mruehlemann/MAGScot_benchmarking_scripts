#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="04c_maxbin"

##########################
####    SETUP     #######
###########################

echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/$dataset/rawdata"


#############################
####    MAXBIN DEFAULT #######
#############################
source activate metagenome_env
module load gcc


samples=($(ls $datafolder/ | grep tar | cut -d '.' -f 1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}
s=$sample

cd $TMPDIR

cat  $workfolder/$dataset/$sample/$sample.megahit.fasta | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta

awk '{if(NR>1) print $1"\t"$3}' $workfolder/$dataset/${s}/${s}.depth.txt  > ${s}.abu

mkdir -p ${s}_maxbin2_out
perl /work_ifs/sukmb276/software/MaxBin-2.2.4/run_MaxBin_prodigal.pl -contig sample.${s}.fasta -prodigal $workfolder/$dataset/${s}/${s}.prodigal.faa -out ${s}_maxbin2_out/${s}.maxbin2 -abund ${s}.abu -thread ${SLURM_CPUS_PER_TASK} -min_contig_length 2000

grep '>' ${s}_maxbin2_out/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/$dataset/${s}/${s}.maxbin2.contigs_to_bin.tsv
