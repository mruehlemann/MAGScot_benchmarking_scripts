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

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/MAGScoT_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

######################################################
######### CDS CALLING AND ANNOTATION     ############
######################################################
source activate metagenome_env

s=$SLURM_ARRAY_TASK_ID

cd $TMPDIR
mkdir ${s}_tmp_workfolder
zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample.${s}.fasta

echo $s preproc_prodigal start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.preproc.log
cat sample.${s}.fasta | ~/Isilon/software/bin/parallel -j ${SLURM_CPUS_PER_TASK} --block 999k --recstart '>' --pipe prodigal -p meta -a ${s}_tmp_workfolder/${s}.{#}.faa -d ${s}_tmp_workfolder/${s}.{#}.ffn -o tmpfile
cat ${s}_tmp_workfolder/$s.*.faa > $workfolder/${dataset}/sample_${s}/sample_${s}.prodigal.faa
cat ${s}_tmp_workfolder/$s.*.ffn > $workfolder/${dataset}/sample_${s}/sample_${s}.prodigal.ffn
echo $s preproc_prodigal end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.preproc.log

echo $s hmm start `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.MAGScoT_preproc.log
hmmsearch -o sample_${s}.hmm.tigr.out --tblout sample_${s}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/tigrfam/tigrfam.hmm $workfolder/$dataset/sample_${s}/sample_${s}.prodigal.faa
#pfam_scan.pl -cpu ${SLURM_CPUS_PER_TASK} -outfile ${sample}.hmm.pfam.hit.out -fasta $outdir/$sample/${sample}.prodigal.faa -dir /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/
hmmsearch -o sample_${s}.hmm.pfam.out --tblout sample_${s}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/Pfam-A.hmm $workfolder/$dataset/sample_${s}/sample_${s}.prodigal.faa
echo $s hmm end `date +%Y%m%d%H%M%S` >> $workfolder/$dataset/sample_${s}/sample_${s}.MAGScoT_preproc.log

cat sample_${s}.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > sample_${s}.tigr
cat sample_${s}.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > sample_${s}.pfam
cat sample_${s}.pfam sample_${s}.tigr > $workfolder/$dataset/sample_${s}/sample_${s}.hmm
