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

workfolder="/home/mruehlemann/biscoreto_benchmark"
cd $workfolder

datafolder="/home/mruehlemann/transfer"

dataset="marine"

#############################
####    MAXBIN DEFAULT #######
#############################
conda activate benchmarking_env
conda activate --stack python27_stack_env

for s in $(seq 5 9); do
cd $workfolder
mkdir -p sample_$s
cd $workfolder/sample_${s}
#module load ruby

SLURM_CPUS_PER_TASK=8
mkdir -p metawrap
cd metawrap

conda activate --stack binning_env
mkdir -p metabat2_bins
cat $datafolder/sample_${s}.metabat2.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/sample_${s}.fasta /dev/stdin --output_path metabat2_bins
mkdir -p maxbin2_bins
cat $datafolder/sample_${s}.maxbin2.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/sample_${s}.fasta /dev/stdin --output_path maxbin2_bins
mkdir -p concoct_bins
cat $datafolder/sample_${s}.concoct.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/sample_${s}.fasta /dev/stdin --output_path concoct_bins
mkdir -p vamb_bins
cat $datafolder/sample_${s}.vamb.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/sample_${s}.fasta /dev/stdin --output_path vamb_bins
conda deactivate


conda activate --stack metawrap_env
mwfolder="/home/mruehlemann/software/metaWRAP"
export PATH=$PATH:$mwfolder/bin

cmdbench -j sample_${s}.metawrap_AB.json $mwfolder/bin/metawrap-modules/bin_refinement.sh -o vamb_and_maxbin -A vamb_bins -B maxbin2_bins -t ${SLURM_CPUS_PER_TASK} -m 240 -c 50

cmdbench -j sample_${s}.metawrap_CD.json $mwfolder/bin/metawrap-modules/bin_refinement.sh -o metabat_and_concoct -A metabat2_bins -B concoct_bins -t ${SLURM_CPUS_PER_TASK} -m 240 -c 50

cmdbench -j sample_${s}.metawrap_all.json $mwfolder/bin/metawrap-modules/bin_refinement.sh -o all -A vamb_and_maxbin/metawrap_50_10_bins -B metabat_and_concoct/metawrap_50_10_bins -t ${SLURM_CPUS_PER_TASK} -m 240 -c 50

mkdir -p $workfolder/$dataset/metawrap_out
cp -r all $workfolder/$dataset/metawrap_out/sample_${s}
conda deactivate


done
