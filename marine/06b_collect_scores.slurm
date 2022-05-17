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

#############################
####    MAXBIN DEFAULT #######
#############################

cd $dataset

mkdir -p amber_input
cat sample_*/amber_files/sample_*_binning_gs.tsv > amber_input/all_binning_gs.tsv
cat sample_*/amber_files/sample_*.gsa2k.binning > amber_input/all.gsa2k.binning
cat sample_*/amber_files/sample_*.concoct.binning > amber_input/all.concoct.binning
cat sample_*/amber_files/sample_*.maxbin2.binning > amber_input/all.maxbin2.binning
cat sample_*/amber_files/sample_*.metabat2.binning > amber_input/all.metabat2.binning
cat sample_*/amber_files/sample_*.vamb.binning > amber_input/all.vamb.binning
cat sample_*/amber_files/sample_*.MAGScoT.binning > amber_input/all.MAGScoT_ms50.binning
cat sample_*/amber_files/sample_*.MAGScoT_nomerge.binning > amber_input/all.MAGScoT_nomerge.binning
cat sample_*/amber_files/sample_*.MAGScoT_ms80.binning > amber_input/all.MAGScoT_ms80.binning
cat sample_*/amber_files/sample_*.dastool.binning > amber_input/all.dastool.binning
cat sample_*/amber_files/sample_*.metawrap.binning > amber_input/all.metawrap.binning

source activate amber_env
amber.py -g amber_input/all_binning_gs.tsv amber_input/all.gsa2k.binning amber_input/all.MAGScoT_ms50.binning amber_input/all.MAGScoT_ms80.binning amber_input/all.MAGScoT_nomerge.binning amber_input/all.dastool.binning amber_input/all.metawrap.binning \
  amber_input/all.concoct.binning   amber_input/all.maxbin2.binning  amber_input/all.metabat2.binning  amber_input/all.vamb.binning \
  -l "gsa2k,bisocreto_ms50,bisocreto_ms80, MAGScoT_nomerge, dastool, metawrap, concoct, maxbin2, metabat2, vamb" \
  -o amber_final_out

conda deactivate

grep '>' */*.prodigal.faa | cut -d '/' -f 1 | uniq -c | awk '{print "proteins\t"$2"\t"$1}' > $dataset.counts.tsv
grep '>' */*.fasta | grep -v megahit | cut -d '/' -f 1 | uniq -c | awk '{print "contigs\t"$2"\t"$1}' >> $dataset.counts.tsv
