# benchmarking
conda create -y -n benchmarking_env python pip
conda activate benchmarking_env
pip install cmdbench
conda deactivate

# BiScoReTo

conda create -n biscoreto_env
conda activate biscoreto_env

conda install -y -c conda-forge mamba
mamba install -y -c conda-forge -c bioconda -c r r-base r-optparse r-dplyr r-readr r-funr hmmer prodigal parallel

mkdir -p software
cd software

git clone https://github.com/mruehlemann/BiScoReTo

cd

conda deactivate

# checkm

conda create -n checkm_env
conda activate checkm_env

conda install -y -c conda-forge mamba
mamba install -y -c conda-forge -c bioconda python pysam checkm-genome

conda deactivate

# dastool

conda create -n dastool_env
conda activate dastool_env

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -y -c conda-forge mamba
mamba install -y -c bioconda das_tool

conda deactivate

# metawrap; only refiner needed

conda create -n metawrap_env
conda activate metawrap_env

conda install -y -c conda-forge mamba
mamba install -y -c conda-forge -c bioconda python pysam checkm-genome biopython checkm-genome pandas

cd software
git clone https://github.com/bxlab/metaWRAP


# metawrap; only refiner needed

conda create -n binning_env
conda activate binning_env

conda install -y -c bioconda concoct
