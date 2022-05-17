#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="mb_01_catalog"

###########################
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

module load R

grep '>' */*.prodigal.faa | cut -d '/' -f 1 | uniq -c | awk '{print "proteins\t"$2"\t"$1}' > $dataset.counts.tsv
grep '>' */*.fasta | grep -v megahit | cut -d '/' -f 1 | uniq -c | awk '{print "contigs\t"$2"\t"$1}' >> $dataset.counts.tsv
