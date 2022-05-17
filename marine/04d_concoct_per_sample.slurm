#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="04d_concoct"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

#############################
####    concoct       #######
#############################

source activate binning_env
module load samtools/1.9
cd $TMPDIR

s=$SLURM_ARRAY_TASK_ID

zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta

cut_up_fasta.py sample.${s}.fasta -c 10000 -o 0 --merge_last -b sample_$s.minimap.bed > sample_${s}.filtered.10k.fna
cp $workfolder/$dataset/sample_${s}/sample_$s.minimap.bam .
samtools index sample_${s}.minimap.bam
concoct_coverage_table.py sample_$s.minimap.bed sample_$s.minimap.bam > sample_${s}.coverage_table.tsv
concoct --composition_file sample_${s}.filtered.10k.fna --coverage_file sample_${s}.coverage_table.tsv -c 1000 -r 151 -t ${SLURM_CPUS_PER_TASK} -l 2000 -s 1234 -i 500 -b sample_${s}

merge_cutup_clustering.py sample_${s}_clustering_gt2000.csv > sample_${s}_clustering_merged.csv

awk -F',' -v sample=sample_$s '{if(NR>1) print sample"_concoct_bin_"$2".fasta\t"$1}'  sample_${s}_clustering_merged.csv > $workfolder/$dataset/sample_${s}/sample_${s}.concoct.contigs_to_bin.tsv

#done
