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

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/$dataset/rawdata"



####################################
####    assembly    ########
####################################
source activate metagenome_env

samples=($(ls $datafolder/ | grep tar | cut -d '.' -f 1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}

cd $workfolder/$dataset/

source activate metagenome_env

cd $TMPDIR

s=$sample
echo $s

minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr  $workfolder/$dataset/${dataset}.catalogue.mmi $datafolder/${sample}_R1.fastq.gz $datafolder/${sample}_R2.fastq.gz > sample_$s.catalogue.minimap.sam
samtools view -F 3584 -b --threads ${SLURM_CPUS_PER_TASK} sample_$s.catalogue.minimap.sam | samtools sort > sample_$s.minimap.catalogue.bam
#done

jgi_summarize_bam_contig_depths --outputDepth $workfolder/$dataset/${s}/${s}.catalogue.depth.txt sample_$s.minimap.catalogue.bam

cat  $workfolder/$dataset/$sample/$sample.megahit.fasta | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta
minimap2 -I100G -d sample_${s}.mmi sample.${s}.fasta # make index
minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr sample_${s}.mmi $datafolder/${sample}_R1.fastq.gz $datafolder/${sample}_R2.fastq.gz > sample_$s.minimap.sam
samtools view -F 3584 -b --threads ${SLURM_CPUS_PER_TASK} sample_$s.minimap.sam | samtools sort > $s.minimap.bam
#done

jgi_summarize_bam_contig_depths --outputDepth $workfolder/$dataset/${s}/${s}.depth.txt $s.minimap.bam
cp $s.minimap.bam $workfolder/$dataset/${s}/
