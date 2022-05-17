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

workfolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2/biscoreto_benchmark"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/CAMI_Challenge_2"

dataset="marine"

#############################
####    MAXBIN DEFAULT #######
#############################


cd $workfolder/$dataset

source activate biscoreto_env
### score dastool final bins
for i in $(seq 0 9); do cd sample_${i}; Rscript ~/Isilon/github/BiScoReTo/BiScoReTo.R -i sample_${i}.dastool.contig_to_bin.out --hmm sample_${i}.hmm -o sample_${i}.dastool --score_only; cd ..; done
## score metwawrap final bins
#for i in $(seq 1 9); do cd sample_${i}; Rscript ~/Isilon/github/BiScoReTo/BiScoReTo.R -i sample_${i}.dastool.contig_to_bin.out --hmm sample_${i}.hmm -o sample_${i}.dastool --score_only; cd ..; done

## score gold standard with BiScoReTo
for i in $(seq 0 9); do cd sample_${i};
awk '{if(NF!=4){next}; print $2"\t"$1"\tGSA"}' $datafolder/$dataset/simulation_short_read/2018.08.15_09.49.32_sample_${i}/contigs/binning_gs.tsv > sample_${i}.gsa.contig_to_bin.out
awk '{if(NF!=4){next}; if($4>=2000) print $2"\t"$1"\tGSA2k"}' $datafolder/$dataset/simulation_short_read/2018.08.15_09.49.32_sample_${i}/contigs/binning_gs.tsv > sample_${i}.gsa_2k.contig_to_bin.out
Rscript ~/Isilon/github/BiScoReTo/BiScoReTo.R -i sample_${i}.gsa.contig_to_bin.out --hmm sample_${i}.hmm -o sample_${i}.gsa --score_only;
Rscript ~/Isilon/github/BiScoReTo/BiScoReTo.R -i sample_${i}.gsa_2k.contig_to_bin.out --hmm sample_${i}.hmm -o sample_${i}.gsa_2k --score_only;

cd ..; done


source deactivate

source activate amber_env

for i in $(seq 0 9); do
cd sample_${i};
mkdir -p amber

awk '{if(substr($1,1,3) == "Otu" && $2>0) print}' $datafolder/$dataset//simulation_short_read/abundance${i}.tsv > amber/genomes_in_sample_${i}.txt

for b in $(cut -f 1 amber/genomes_in_sample_${i}.txt); do
b2=$(echo $b | sed 's/[.]/\[.\]/g')
grep -w $b2 $datafolder/$dataset/simulation_short_read/2018.08.15_09.49.32_sample_${i}/contigs/binning_gs.tsv | awk -v b=$b 'BEGIN{totsize=0;size2k=0}{totsize=totsize+$4;if($4>=2000) size2k=size2k+$4}END{print b"\t"totsize"\t"size2k"\t"size2k/totsize}'
done | tee amber/sample_${i}_gs_binstats

head -n 5 $datafolder/$dataset/simulation_short_read/2018.08.15_09.49.32_sample_${i}/contigs/binning_gs.tsv | sed "s/marine/marine_sample_${i}/" > amber/sample_${i}_binning_gs
awk '{if($4>=.5) print $1}' amber/sample_${i}_gs_binstats | sed 's/[.]/\[.\]/g'  | xargs -I {} grep -w {} $datafolder/$dataset/simulation_short_read/2018.08.15_09.49.32_sample_${i}/contigs/binning_gs.tsv >> amber/sample_${i}_binning_gs

echo "@SampleID:marine_sample_${i}" > amber/sample_${i}_gsa2k.binning
echo -e "@@SEQUENCEID\tBINID" >> amber/sample_${i}_gsa2k.binning
awk '{if(NR>6 && $4 >= 2000) print}' amber/sample_${i}_binning_gs >> amber/sample_${i}_gsa2k.binning




### FORMAT:
# @SampleID:marmgCAMI2_short_read_pooled_gold_standard_assembly
# @@SEQUENCEID    BINID
echo "@SampleID:marine_sample_${i}" > amber/sample_${i}_biscoreto.binning
echo -e "@@SEQUENCEID\tBINID" >> amber/sample_${i}_biscoreto.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${i}.biscoreto.refined.contig_to_bin.out >> amber/sample_${i}_biscoreto.binning

# @@SEQUENCEID    BINID
echo "@SampleID:marine_sample_${i}" > amber/sample_${i}_biscoreto_nomerge.binning
echo -e "@@SEQUENCEID\tBINID" >> amber/sample_${i}_biscoreto_nomerge.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${i}.biscoreto_nomerge.refined.contig_to_bin.out >> amber/sample_${i}_biscoreto.binning

echo "@SampleID:marine_sample_${i}" > amber/sample_${i}_dastool.binning
echo -e "@@SEQUENCEID\tBINID" >> amber/sample_${i}_dastool.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${i}.dastool.contig_to_bin.out >> amber/sample_${i}_dastool.binning

echo "@SampleID:marine_sample_${i}" > amber/sample_${i}_vamb.binning
echo -e "@@SEQUENCEID\tBINID" >> amber/sample_${i}_vamb.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${i}.vamb.contigs_to_bin.tsv >> amber/sample_${i}_vamb.binning

echo "@SampleID:marine_sample_${i}" > amber/sample_${i}_concoct.binning
echo -e "@@SEQUENCEID\tBINID" >> amber/sample_${i}_concoct.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${i}.concoct.contigs_to_bin.tsv >> amber/sample_${i}_concoct.binning

echo "@SampleID:marine_sample_${i}" > amber/sample_${i}_maxbin2.binning
echo -e "@@SEQUENCEID\tBINID" >> amber/sample_${i}_maxbin2.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${i}.maxbin2.contigs_to_bin.tsv >> amber/sample_${i}_maxbin2.binning

echo "@SampleID:marine_sample_${i}" > amber/sample_${i}_metabat2.binning
echo -e "@@SEQUENCEID\tBINID" >> amber/sample_${i}_metabat2.binning
awk '{if(NR>1) print $2"\t"$1}' sample_${i}.metabat2.contigs_to_bin.tsv >> amber/sample_${i}_metabat2.binning


#amber.py -g amber/sample_${i}_binning_gs amber/sample_${i}_biscoreto.binning amber/sample_${i}_dastool.binning amber/sample_${i}_vamb.binning amber/sample_${i}_concoct.binning amber/sample_${i}_maxbin2.binning amber/sample_${i}_metabat2.binning \
# -l "biscoreto,dastool,vamb,concoct,maxbin2,metabat2" -o results/


 cd $workfolder/$dataset

done

mkdir amber_all
cat sample_*/amber/sample_*_binning_gs > amber_all/all_binning_gs
cat sample_*/amber/sample_*_gsa2k.binning > amber_all/gsa2k.binning
cat sample_*/amber/sample_*_biscoreto.binning > amber_all/biscoreto.binning
cat sample_*/amber/sample_*_biscoreto_nomerge.binning > amber_all/biscoreto_nomerge.binning
cat sample_*/amber/sample_*_concoct.binning > amber_all/concoct.binning
cat sample_*/amber/sample_*_dastool.binning > amber_all/dastool.binning
cat sample_*/amber/sample_*_vamb.binning > amber_all/vamb.binning
cat sample_*/amber/sample_*_maxbin2.binning > amber_all/maxbin2.binning
cat sample_*/amber/sample_*_metabat2.binning > amber_all/metabat2.binning

amber.py -g amber_all/all_binning_gs amber_all/gsa2k.binning amber_all/biscoreto.binning amber_all/biscoreto_nomerge.binning amber_all/dastool.binning amber_all/vamb.binning amber_all/concoct.binning amber_all/maxbin2.binning amber_all/metabat2.binning \
 -l "gsa2k,biscoreto,biscoreto_nomerge,dastool,vamb,concoct,maxbin2,metabat2" -o amber_results/



 module load R

 Rscript
