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

workfolder="/home/mruehlemann/biscoreto_benchmark"

#cd $workfolder

datafolder="/home/mruehlemann/transfer"

dataset="marine"

#############################
####    MAXBIN DEFAULT #######
#############################
conda activate benchmarking_env

for s in $(seq 1 9); do
cd $workfolder
mkdir -p sample_$s
cd $workfolder/sample_${s}
conda activate --stack dastool_env
#module load ruby

SLURM_CPUS_PER_TASK=8

mkdir -p das_tool
cd das_tool

awk '{print $2"\t"$1}' $datafolder/sample_${s}.maxbin2.contigs_to_bin.tsv > maxbin2.out
awk '{print $2"\t"$1}' $datafolder/sample_${s}.metabat2.contigs_to_bin.tsv > metabat2.out
awk '{print $2"\t"$1}' $datafolder/sample_${s}.concoct.contigs_to_bin.tsv > concoct.out
awk '{print $2"\t"$1}' $datafolder/sample_${s}.vamb.contigs_to_bin.tsv > vamb.out

cmdbench -j sample_${s}.dastool.json DAS_Tool -i maxbin2.out,metabat2.out,concoct.out,vamb.out -c $datafolder/sample_${s}.fasta -o sample_${s}.dastool --labels=maxbin,metabat,concoct,vamb -t ${SLURM_CPUS_PER_TASK}

conda deactivate

conda activate --stack binning_env
mkdir -p dastool_cleanbins
cat sample_${s}.dastool_DASTool_contig2bin.tsv | awk -vs=$s 'BEGIN{i=0; bid=""}{if($2!=bid){i=i+1;bid=$2;nbid="sample_"s".dastoolbin_"i".fasta"}; print nbid"\t"$1"\tDASTool"}' > sample_${s}.dastool.contig_to_bin.out
cat sample_${s}.dastool.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/sample_${s}.fasta /dev/stdin --output_path dastool_cleanbins

conda deactivate

# 
# conda activate --stack checkm_env
#
# checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f ${s}.dastool.checkm.out --tab_table dastool_cleanbins checkm_out
# rm -r checkm
# conda deactivate

done
