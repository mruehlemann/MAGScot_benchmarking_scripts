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

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

#############################
####    MAXBIN DEFAULT #######
#############################
cd $TMPDIR

s=$SLURM_ARRAY_TASK_ID

source activate dastool_env
module load ruby

export PATH=$PATH:/work_ifs/sukmb276/software/DAS_Tool/

cd $workfolder/$dataset/sample_${s}

awk '{print $2"\t"$1}' sample_${s}.maxbin2.contigs_to_bin.tsv > $TMPDIR/maxbin2.out
awk '{print $2"\t"$1}' sample_${s}.metabat2.contigs_to_bin.tsv > $TMPDIR/metabat2.out
awk '{print $2"\t"$1}' sample_${s}.concoct.contigs_to_bin.tsv > $TMPDIR/concoct.out
awk '{print $2"\t"$1}' sample_${s}.vamb.contigs_to_bin.tsv > $TMPDIR/vamb.out

cd $TMPDIR
cp $workfolder/$dataset/sample_${s}/sample_${s}.prodigal.faa .
zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample_${s}.fasta

echo sample_$s dastool start `date +%Y%m%d%H%M%S` > $workfolder/$dataset/sample_${s}/sample_${s}.dastool.log

DAS_Tool -i maxbin2.out,metabat2.out,concoct.out,vamb.out -c sample_${s}.fasta -o sample_${s}.dastool -p sample_${s}.prodigal.faa --write_bin_evals --labels=maxbin,metabat,concoct,vamb -t ${SLURM_CPUS_PER_TASK}

echo sample_$s dastool end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.dastool.log


conda deactivate
source activate binning_env
mkdir -p $workfolder/$dataset/sample_${s}/dastool_cleanbins
cat sample_${s}.dastool_DASTool_contig2bin.tsv | awk -vs=$s 'BEGIN{i=0; bid=""}{if($2!=bid){i=i+1;bid=$2;nbid="sample_"s".dastoolbin_"i".fasta"}; print nbid"\t"$1"\tDASTool"}' > $workfolder/$dataset/sample_${s}/sample_${s}.dastool.contig_to_bin.out
cat $workfolder/$dataset/sample_${s}/sample_${s}.dastool.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py sample_${s}.fasta /dev/stdin --output_path $workfolder/$dataset/sample_${s}/dastool_cleanbins

conda deactivate
# source activate metagenome_env
#
# echo sample_$s dastool_checkm start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.dastool.log
#
# checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f $workfolder/samples/${s}/${s}.dastool.checkm.out --tab_table $workfolder/$dataset/sample_${s}/dastool_cleanbins checkm_out
#
# echo sample_$s dastool_checkm end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.dastool.log
