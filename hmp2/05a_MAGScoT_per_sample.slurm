#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=100gb
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

cd $workfolder/$dataset/${s}

grep "^${sample}_" $workfolder/$dataset/VAMB_out/clusters.tsv > ${s}.vamb.contigs_to_bin.tsv

awk '{print $1"\t"$2"\tvamb"}'  ${s}.vamb.contigs_to_bin.tsv > ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tconcoct"}'  ${s}.concoct.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tmetabat2"}'  ${s}.metabat2.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tmaxbin2"}'  ${s}.maxbin2.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv

source activate MAGScoT_env
cd $TMPDIR

# echo $s MAGScoT_hmm start `date +%Y%m%d%H%M%S` > $workfolder/$dataset/${s}/${s}.MAGScoT.log
#
# hmmsearch -o ${s}.hmm.tigr.out --tblout ${s}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/tigrfam/tigrfam.hmm $workfolder/$dataset/${s}/${s}.prodigal.faa
# #pfam_scan.pl -cpu ${SLURM_CPUS_PER_TASK} -outfile ${sample}.hmm.pfam.hit.out -fasta $outdir/$sample/${sample}.prodigal.faa -dir /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/
# hmmsearch -o ${s}.hmm.pfam.out --tblout ${s}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/Pfam-A.hmm $workfolder/$dataset/${s}/${s}.prodigal.faa
#
# cat ${s}.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > ${s}.tigr
# cat ${s}.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > ${s}.pfam
# cat ${s}.pfam ${s}.tigr > $workfolder/$dataset/${s}/${s}.hmm
#
# echo $s MAGScoT_hmm end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.MAGScoT.log

cd $workfolder/$dataset/${s}

echo $s MAGScoT start `date +%Y%m%d%H%M%S` > $workfolder/$dataset/${s}/${s}.MAGScoT.log
Rscript /work_ifs/sukmb276/github/MAGScoT/MAGScoT.R -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm -o ${s}.MAGScoT_sm50
echo $s MAGScoT end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.MAGScoT.log

echo $s MAGScoT_nomerge start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.MAGScoT.log
Rscript /work_ifs/sukmb276/github/MAGScoT/MAGScoT.R -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm -o ${s}.MAGScoT_nomerge --skip_merge_bins
echo $s MAGScoT_nomerge end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.MAGScoT.log

echo $s MAGScoT_ms80 start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.MAGScoT.log
Rscript /work_ifs/sukmb276/github/MAGScoT/MAGScoT.R -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm -o ${s}.MAGScoT_sm80 -s 0.8
echo $s MAGScoT_ms80 end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.MAGScoT.log

source deactivate

source activate binning_env
mkdir -p MAGScoT_cleanbins
cat  $workfolder/$dataset/$sample/$sample.megahit.fasta | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > ${s}.fasta

cat ${s}.MAGScoT.refined.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py ${s}.fasta /dev/stdin --output_path MAGScoT_cleanbins
source deactivate

#
# source activate metagenome_env
# cd $TMPDIR
# echo $s MAGScoT_checkm start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.MAGScoT.log
#
# checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f $workfolder/$dataset/${s}/${s}.MAGScoT.checkm.out --tab_table $workfolder/$dataset/${s}/MAGScoT_cleanbins checkm_out
#
# echo $s MAGScoT_checkm end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.MAGScoT.log
