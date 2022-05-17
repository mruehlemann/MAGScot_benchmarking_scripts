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

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/$dataset/rawdata"


source activate metagenome_env

cd $TMPDIR

for s in $(ls $datafolder/ | grep tar | cut -d '.' -f 1); do
#ava=$(grep $s $workfolder/samples_available.txt)
echo $s
cat $workfolder/$dataset/$s/$s.megahit.fasta | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta
done

cat sample.*.fasta > ${dataset}.catalogue.fna

minimap2 -I100G -d ${dataset}.catalogue.mmi ${dataset}.catalogue.fna # make index

gzip ${dataset}.catalogue.fna

mkdir -p $workfolder/${dataset}
cp ${dataset}.catalogue.mmi ${dataset}.catalogue.fna.gz  $workfolder/${dataset}
