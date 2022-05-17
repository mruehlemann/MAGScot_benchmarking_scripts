#!/bin/bash
#SBATCH -c 12
#SBATCH --mem=120gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="02_mapping"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

###########################
####    MAPPING     #######
###########################
source activate metagenome_env

mkdir -p bams

cd $TMPDIR

s=$SLURM_ARRAY_TASK_ID
echo $s

mkdir -p $workfolder/$dataset/sample_${s}/
minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr  $workfolder/$dataset/${dataset}.catalogue.mmi $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/reads/anonymous_reads.fq.gz > sample_$s.catalogue.minimap.sam
samtools view -F 3584 -b --threads ${SLURM_CPUS_PER_TASK} sample_$s.catalogue.minimap.sam | samtools sort > sample_$s.minimap.catalogue.bam
#done

jgi_summarize_bam_contig_depths --outputDepth $workfolder/$dataset/sample_${s}/sample_${s}.catalogue.depth.txt sample_$s.minimap.catalogue.bam

zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta
minimap2 -I100G -d sample_${s}.mmi sample.${s}.fasta # make index
minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr sample_${s}.mmi $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/reads/anonymous_reads.fq.gz > sample_$s.minimap.sam
samtools view -F 3584 -b --threads ${SLURM_CPUS_PER_TASK} sample_$s.minimap.sam | samtools sort > sample_$s.minimap.bam
#done

jgi_summarize_bam_contig_depths --outputDepth $workfolder/$dataset/sample_${s}/sample_${s}.depth.txt sample_$s.minimap.bam
cp sample_$s.minimap.bam $workfolder/$dataset/sample_${s}/
