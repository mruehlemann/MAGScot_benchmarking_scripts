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
cd $workfolder

datafolder="/home/mruehlemann/transfer"

dataset="marine"

#############################
####    MAXBIN DEFAULT #######
#############################
conda activate benchmarking_env

for s in $(seq 0 8); do
cd $workfolder
mkdir -p sample_$s
cd $workfolder/sample_${s}
conda activate --stack biscoreto_env

biscofolder="/home/mruehlemann/software/BiScoReTo"
SLURM_CPUS_PER_TASK=8

mkdir -p bisocreto
cd bisocreto

# awk '{print $1"\t"$2"\tvamb"}'  $datafolder/sample_${s}.vamb.contigs_to_bin.tsv > sample_${s}.contigs_to_bin.tsv
# awk '{print $1"\t"$2"\tconcoct"}'  $datafolder/sample_${s}.concoct.contigs_to_bin.tsv >> sample_${s}.contigs_to_bin.tsv
# awk '{print $1"\t"$2"\tmetabat2"}'  $datafolder/sample_${s}.metabat2.contigs_to_bin.tsv >> sample_${s}.contigs_to_bin.tsv
# awk '{print $1"\t"$2"\tmaxbin2"}'  $datafolder/sample_${s}.maxbin2.contigs_to_bin.tsv >> sample_${s}.contigs_to_bin.tsv
#
#
# mkdir ${s}_tmp_workfolder
#
# cat $datafolder/sample_${s}.fasta | cmdbench -j sample_${s}.prodigal.json "parallel -j ${SLURM_CPUS_PER_TASK} --block 999k --recstart '>' --pipe prodigal -p meta -a ${s}_tmp_workfolder/${s}.{#}.faa -d ${s}_tmp_workfolder/${s}.{#}.ffn -o tmpfile 2> /dev/null"
#
# cat ${s}_tmp_workfolder/$s.*.faa > sample_${s}.prodigal.faa
# cat ${s}_tmp_workfolder/$s.*.ffn > sample_${s}.prodigal.ffn
# rm -r ${s}_tmp_workfolder/
#
# cmdbench -j sample_${s}.hmm_tigr.json hmmsearch -o sample_${s}.hmm.tigr.out --tblout sample_${s}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} $biscofolder/hmm/gtdbtk_rel207_tigrfam.hmm sample_${s}.prodigal.faa
# #pfam_scan.pl -cpu ${SLURM_CPUS_PER_TASK} -outfile ${sample}.hmm.pfam.hit.out -fasta $outdir/$sample/${sample}.prodigal.faa -dir /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/
# cmdbench -j sample_${s}.hmm_pfam.json hmmsearch -o sample_${s}.hmm.pfam.out --tblout sample_${s}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} $biscofolder/hmm/gtdbtk_rel207_Pfam-A.hmm sample_${s}.prodigal.faa
#
# cat sample_${s}.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > sample_${s}.tigr
# cat sample_${s}.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > sample_${s}.pfam
# cat sample_${s}.pfam sample_${s}.tigr > sample_${s}.hmm

cmdbench -j sample_${s}.biscoreto.json "Rscript $biscofolder/BiScoReTo.R -i sample_${s}.contigs_to_bin.tsv --hmm sample_${s}.hmm -o sample_${s}.biscoreto"
cmdbench -j sample_${s}.biscoreto_nomerge.json "Rscript $biscofolder/BiScoReTo.R -i sample_${s}.contigs_to_bin.tsv --hmm sample_${s}.hmm -o sample_${s}.biscoreto_nomerge --skip_merge_bins"

conda deactivate

# conda activate --stack binning_env
# mkdir -p biscoreto_cleanbins
# cat sample_${s}.biscoreto.refined.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/sample_${s}.fasta /dev/stdin --output_path biscoreto_cleanbins
# conda deactivate


# conda activate --stack checkm_env
#
# cmdbench -j sample_${s}.biscoreto_checkm.json checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f sample_${s}.biscoreto.checkm.out --tab_table biscoreto_cleanbins checkm_out
#
# rm -r checkm_out
# conda deactivate

done
