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

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/$dataset/rawdata"

cd $workfolder

##################################
####    start collection   ########
##################################

samples=($(ls $datafolder/ | grep tar | cut -d '.' -f 1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}
s=$sample


cd $workfolder/$dataset/${s}

cat  $workfolder/$dataset/$sample/$sample.megahit.fasta | cut -f 1 -d ' ' |  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print $1"\t"length($2)"\n"'} > ${s}.contig_lengths.tsv

cat ../metawrap_out/${s}/metawrap_50_10_bins.contigs | awk '{print $2"\t"$1"\tmetawrap"}' > ${s}.metawrap.contigs_to_bin.tsv

module load R

Rscript $workfolder/scripts/MAGScoT_benchmark_scripts/hmp2/06a_score_bins.R

source activate MAGScoT_env
Rscript /work_ifs/sukmb276/github/MAGScoT/MAGScoT.sR -i ${s}.final_ctb.tsv --hmm ${s}.hmm -o ${s}.final --score_only
conda deactivate

cd $TMPDIR
source activate binning_env
mkdir -p finalbins
cat  $workfolder/$dataset/$sample/$sample.megahit.fasta | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > ${s}.fasta

cat $workfolder/$dataset/${s}/${s}.final_ctb.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py ${s}.fasta /dev/stdin --output_path finalbins
conda deactivate

# score with checkm
conda activate checkm_env

checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f $workfolder/$dataset/${s}/${s}.final.checkm.out --tab_table finalbins checkm_out
rm -r checkm
conda deactivate
