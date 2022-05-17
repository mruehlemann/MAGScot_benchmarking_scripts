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
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

###########################
####    METABAT     #######
###########################
s=$SLURM_ARRAY_TASK_ID

cd $workfolder/$dataset/sample_${s}

mkdir -p amber_files
awk '{if($2 > 0 && substr($1,1,3)=="Otu") print $1}'  $datafolder/$dataset/simulation_short_read/abundance0.tsv > amber_files/otus_in_sample_$s.txt
head -n 5 $datafolder/$dataset/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/binning_gs.tsv | sed "s/$dataset/$dataset\_sample_${s}/" > amber_files/sample_${s}_binning_gs.tsv
for o in $(cat amber_files/otus_in_sample_$s.txt | sed 's/[.]/\[.\]/'); do
  grep -w $o $datafolder/$dataset/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/binning_gs.tsv | awk -v o=$o 'BEGIN{totsize=0; size2k=0} {totsize=totsize+$4; if($4>=2000) size2k=size2k+$4} END{print o"\t"totsize"\t"size2k"\t"size2k/totsize"\t"}'
done | tee amber_files/sample_${s}.binstats

for o in $(awk '{if($4>.4) print $1}' amber_files/sample_${s}.binstats); do
  grep -w $o $datafolder/$dataset/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/binning_gs.tsv
done >> amber_files/sample_${s}_binning_gs.tsv

### Gold standard with contigs > 2kb
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.gsa2k.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.gsa2k.binning
awk '!/^@/ {if($4 >= 2000) print $1"\t"$2}' amber_files/sample_${s}_binning_gs.tsv >> amber_files/sample_${s}.gsa2k.binning

### concoct
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.concoct.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.concoct.binning
awk '{print $2"\t"$1}' sample_${s}.concoct.contigs_to_bin.tsv >> amber_files/sample_${s}.concoct.binning

### maxbin2
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.maxbin2.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.maxbin2.binning
awk '{print $2"\t"$1}' sample_${s}.maxbin2.contigs_to_bin.tsv >> amber_files/sample_${s}.maxbin2.binning

### metabat2
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.metabat2.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.metabat2.binning
awk '{print $2"\t"$1}' sample_${s}.metabat2.contigs_to_bin.tsv >> amber_files/sample_${s}.metabat2.binning

### vamb
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.vamb.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.vamb.binning
awk '{print $2"\t"$1}' sample_${s}.vamb.contigs_to_bin.tsv >> amber_files/sample_${s}.vamb.binning

### MAGScoT
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.MAGScoT.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.MAGScoT.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${s}.MAGScoT.refined.contig_to_bin.out >> amber_files/sample_${s}.MAGScoT.binning

### MAGScoT no merge
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.MAGScoT_nomerge.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.MAGScoT_nomerge.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${s}.MAGScoT_nomerge.refined.contig_to_bin.out >> amber_files/sample_${s}.MAGScoT_nomerge.binning

### MAGScoT min sharing 80%
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.MAGScoT_ms80.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.MAGScoT_ms80.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${s}.MAGScoT_ms80.refined.contig_to_bin.out >> amber_files/sample_${s}.MAGScoT_ms80.binning

### Dastool
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.dastool.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.dastool.binning
awk '{print $2"\t"$1}' sample_${s}.dastool.contig_to_bin.out >> amber_files/sample_${s}.dastool.binning

cat ../metawrap_out/sample_${s}/all/metawrap_50_10_bins.contigs | awk '{print $2"\t"$1"\tmetawrap"}' > sample_${s}.metawrap.contigs_to_bin.tsv
## metawrap
echo -e "@SampleID:marine_sample_${s}" > amber_files/sample_${s}.metawrap.binning
echo -e "@@SEQUENCEID\tBINID" >> amber_files/sample_${s}.metawrap.binning
awk '{print $2"\t"$1}' sample_${s}.metawrap.contigs_to_bin.tsv >> amber_files/sample_${s}.metawrap.binning

### gold standard binning files
cat amber_files/sample_${s}.gsa2k.binning | grep -v '@' | awk '{print $2"_gsa2k\t"$1"\tgsa2k"}' | sed 's/[.]/\_/g' > sample_${s}.gs.contig_to_bin.tsv

zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' |  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print $1"\t"length($2)"\n"'} > sample_${s}.contig_lengths.tsv

module load R

Rscript $workfolder/scripts/MAGScoT_benchmark_scripts/$dataset/06a_score_bins.R

source activate MAGScoT_env
Rscript /work_ifs/sukmb276/github/MAGScoT/MAGScoT.R -i sample_${s}.final_ctb.tsv --hmm sample_${s}.hmm -o sample_${s}.final --score_only
conda deactivate

cd $TMPDIR
source activate binning_env
mkdir -p finalbins
zcat $datafolder/${dataset}/simulation_short_read/2018.08.15_09.49.32_sample_${s}/contigs/anonymous_gsa.fasta.gz | cut -f 1 -d ' ' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > sample_${s}.fasta

cat $workfolder/$dataset/sample_${s}/sample_${s}.final_ctb.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py sample_${s}.fasta /dev/stdin --output_path finalbins
conda deactivate

# score with checkm
conda activate checkm_env

checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f $workfolder/$dataset/sample_${s}/sample_${s}.final.checkm.out --tab_table finalbins checkm_out
rm -r checkm
conda deactivate
