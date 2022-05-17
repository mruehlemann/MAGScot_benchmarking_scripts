#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="04d_concoct"

##########################
####    SETUP     #######
###########################

echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/$dataset/rawdata"


#############################
####    concoct       #######
#############################

source activate binning_env
module load samtools/1.9
cd $TMPDIR


samples=($(ls $datafolder/ | grep tar | cut -d '.' -f 1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}
s=$sample

cat  $workfolder/$dataset/$sample/$sample.megahit.fasta | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta

cut_up_fasta.py sample.${s}.fasta -c 10000 -o 0 --merge_last -b $s.minimap.bed > ${s}.filtered.10k.fna
cp $workfolder/$dataset/${s}/$s.minimap.bam .
samtools index ${s}.minimap.bam
concoct_coverage_table.py $s.minimap.bed $s.minimap.bam > ${s}.coverage_table.tsv
concoct --composition_file ${s}.filtered.10k.fna --coverage_file ${s}.coverage_table.tsv -c 1000 -r 151 -t ${SLURM_CPUS_PER_TASK} -l 2000 -s 1234 -i 500 -b ${s}

merge_cutup_clustering.py ${s}_clustering_gt2000.csv > ${s}_clustering_merged.csv

awk -F',' -v sample=$s '{if(NR>1) print sample"_concoct_bin_"$2".fasta\t"$1}'  ${s}_clustering_merged.csv > $workfolder/$dataset/${s}/${s}.concoct.contigs_to_bin.tsv

#done
