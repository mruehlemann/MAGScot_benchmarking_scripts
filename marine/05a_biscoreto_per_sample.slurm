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

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/biscoreto_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

#############################
####    MAXBIN DEFAULT #######
#############################

s=$SLURM_ARRAY_TASK_ID

cd $workfolder/$dataset/sample_${s}

grep "^S${s}C" $workfolder/$dataset/VAMB_out/clusters.tsv > sample_${s}.vamb.contigs_to_bin.tsv

awk '{print $1"\t"$2"\tvamb"}'  sample_${s}.vamb.contigs_to_bin.tsv > sample_${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tconcoct"}'  sample_${s}.concoct.contigs_to_bin.tsv >> sample_${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tmetabat2"}'  sample_${s}.metabat2.contigs_to_bin.tsv >> sample_${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tmaxbin2"}'  sample_${s}.maxbin2.contigs_to_bin.tsv >> sample_${s}.contigs_to_bin.tsv

source activate biscoreto_env
cd $TMPDIR

echo sample_$s biscoreto_hmm start `date +%Y%m%d%H%M%S` > $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log

hmmsearch -o sample_${s}.hmm.tigr.out --tblout sample_${s}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/tigrfam/tigrfam.hmm $workfolder/$dataset/sample_${s}/sample_${s}.prodigal.faa
#pfam_scan.pl -cpu ${SLURM_CPUS_PER_TASK} -outfile ${sample}.hmm.pfam.hit.out -fasta $outdir/$sample/${sample}.prodigal.faa -dir /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/
hmmsearch -o sample_${s}.hmm.pfam.out --tblout sample_${s}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/Pfam-A.hmm $workfolder/$dataset/sample_${s}/sample_${s}.prodigal.faa

cat sample_${s}.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > sample_${s}.tigr
cat sample_${s}.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > sample_${s}.pfam
cat sample_${s}.pfam sample_${s}.tigr > $workfolder/$dataset/sample_${s}/sample_${s}.hmm

echo sample_$s biscoreto_hmm end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log

cd $workfolder/$dataset/sample_${s}

echo sample_$s biscoreto start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log
Rscript /work_ifs/sukmb276/github/BiScoReTo/BiScoReTo.R -i sample_${s}.contigs_to_bin.tsv --hmm sample_${s}.hmm -o sample_${s}.biscoreto
echo sample_$s biscoreto end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log

echo sample_$s biscoreto_nomerge start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log
Rscript /work_ifs/sukmb276/github/BiScoReTo/BiScoReTo.R -i sample_${s}.contigs_to_bin.tsv --hmm sample_${s}.hmm -o sample_${s}.biscoreto_nomerge --skip_merge_bins
echo sample_$s biscoreto_nomerge end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log

echo sample_$s biscoreto_ms80 start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log
Rscript /work_ifs/sukmb276/github/BiScoReTo/BiScoReTo.R -i sample_${s}.contigs_to_bin.tsv --hmm sample_${s}.hmm -o sample_${s}.biscoreto_ms80 -s 0.8
echo sample_$s biscoreto_ms80 end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log

source deactivate

source activate binning_env
mkdir -p biscoreto_cleanbins
zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample_${s}.fasta

cat sample_${s}.biscoreto.refined.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py sample_${s}.fasta /dev/stdin --output_path biscoreto_cleanbins

# source activate metagenome_env
# cd $TMPDIR
# echo sample_$s biscoreto_checkm start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log
#
# checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.checkm.out --tab_table $workfolder/$dataset/sample_${s}/biscoreto_cleanbins checkm_out
#
# echo sample_$s biscoreto_checkm end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log


for
echo sample_$s biscoreto_ms80 start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log
Rscript /work_ifs/sukmb276/github/BiScoReTo/BiScoReTo.R -i sample_${s}.contigs_to_bin.tsv --hmm sample_${s}.hmm -o sample_${s}.biscoreto_ms80 -s 0.8
echo sample_$s biscoreto_ms80 end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.biscoreto.log
