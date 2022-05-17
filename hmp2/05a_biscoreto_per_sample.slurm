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

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/biscoreto_benchmark"
cd $workfolder

dataset="hmpgut"

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/biscoreto_benchmark/$dataset/rawdata"


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

source activate biscoreto_env
cd $TMPDIR

# echo $s biscoreto_hmm start `date +%Y%m%d%H%M%S` > $workfolder/$dataset/${s}/${s}.biscoreto.log
#
# hmmsearch -o ${s}.hmm.tigr.out --tblout ${s}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/tigrfam/tigrfam.hmm $workfolder/$dataset/${s}/${s}.prodigal.faa
# #pfam_scan.pl -cpu ${SLURM_CPUS_PER_TASK} -outfile ${sample}.hmm.pfam.hit.out -fasta $outdir/$sample/${sample}.prodigal.faa -dir /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/
# hmmsearch -o ${s}.hmm.pfam.out --tblout ${s}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/Pfam-A.hmm $workfolder/$dataset/${s}/${s}.prodigal.faa
#
# cat ${s}.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > ${s}.tigr
# cat ${s}.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > ${s}.pfam
# cat ${s}.pfam ${s}.tigr > $workfolder/$dataset/${s}/${s}.hmm
#
# echo $s biscoreto_hmm end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.biscoreto.log

cd $workfolder/$dataset/${s}

echo $s biscoreto start `date +%Y%m%d%H%M%S` > $workfolder/$dataset/${s}/${s}.biscoreto.log
Rscript /work_ifs/sukmb276/github/BiScoReTo/BiScoReTo.R -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm -o ${s}.biscoreto_sm50
echo $s biscoreto end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.biscoreto.log

echo $s biscoreto_nomerge start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.biscoreto.log
Rscript /work_ifs/sukmb276/github/BiScoReTo/BiScoReTo.R -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm -o ${s}.biscoreto_nomerge --skip_merge_bins
echo $s biscoreto_nomerge end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.biscoreto.log

echo $s biscoreto_ms80 start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.biscoreto.log
Rscript /work_ifs/sukmb276/github/BiScoReTo/BiScoReTo.R -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm -o ${s}.biscoreto_sm80 -s 0.8
echo $s biscoreto_ms80 end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.biscoreto.log

source deactivate

source activate binning_env
mkdir -p biscoreto_cleanbins
cat  $workfolder/$dataset/$sample/$sample.megahit.fasta | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > ${s}.fasta

cat ${s}.biscoreto.refined.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py ${s}.fasta /dev/stdin --output_path biscoreto_cleanbins
source deactivate

#
# source activate metagenome_env
# cd $TMPDIR
# echo $s biscoreto_checkm start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.biscoreto.log
#
# checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f $workfolder/$dataset/${s}/${s}.biscoreto.checkm.out --tab_table $workfolder/$dataset/${s}/biscoreto_cleanbins checkm_out
#
# echo $s biscoreto_checkm end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/${s}/${s}.biscoreto.log
