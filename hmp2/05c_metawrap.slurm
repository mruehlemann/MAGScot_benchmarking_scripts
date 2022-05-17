#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/log/%A_%a.out
#SBATCH --job-name="mb_01_catalog"

##########################
####    SETUP     #######
###########################

echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark/$dataset/rawdata"


##################################
####    start collection   ########
##################################

samples=($(ls $datafolder/ | grep tar | cut -d '.' -f 1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}
s=$sample

cd $TMPDIR

cat  $workfolder/$dataset/$sample/$sample.megahit.fasta | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > ${s}.fasta

source activate binning_env
mkdir -p metabat2_bins
cat $workfolder/$dataset/${s}/${s}.metabat2.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py ${s}.fasta /dev/stdin --output_path metabat2_bins
mkdir -p maxbin2_bins
cat $workfolder/$dataset/${s}/${s}.maxbin2.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py ${s}.fasta /dev/stdin --output_path maxbin2_bins
mkdir -p concoct_bins
cat $workfolder/$dataset/${s}/${s}.concoct.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py ${s}.fasta /dev/stdin --output_path concoct_bins
mkdir -p vamb_bins
cat $workfolder/$dataset/${s}/${s}.vamb.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py ${s}.fasta /dev/stdin --output_path vamb_bins
conda deactivate

source activate metawrap_env
export PATH=$PATH:/work_ifs/sukmb276/github/metaWRAP/bin

echo $s metawrap_AB start `date +%Y%m%d%H%M%S` > $workfolder/$dataset/${s}/${s}.metawrap.log
/work_ifs/sukmb276/github/metaWRAP/bin/metawrap-modules/bin_refinement2.sh -o vamb_and_maxbin -A vamb_bins -B maxbin2_bins -t ${SLURM_CPUS_PER_TASK} -m 80 -c 50
echo $s metawrap_AB end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.metawrap.log

echo $s metawrap_ABCD start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.metawrap.log
/work_ifs/sukmb276/github/metaWRAP/bin/metawrap-modules/bin_refinement2.sh -o all -A vamb_and_maxbin/metawrap_50_10_bins -B concoct_bins -C metabat2_bins -t ${SLURM_CPUS_PER_TASK} -m 80 -c 50
echo $s metawrap_ABCD end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.metawrap.log

mkdir -p $workfolder/$dataset/metawrap_out/${s}
cp -r all $workfolder/$dataset/metawrap_out/${s}
cat all/metawrap_50_10_bins.contigs | awk '{print $2"\t"$1"\tmetawrap"}' >  $workfolder/$dataset/${s}/${s}.contig_to_bin.out
