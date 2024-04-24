#!/bin/bash -l

#$ -S /bin/bash
#$ -N Parser
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q


module load miniconda/4.14.0
source activate TaxProfiler

set -x

TaxprofilerOutputFolder=$1
ParsingOutputFolder=$2

source /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Bash_Config_Taxprofiler.config

nextflow run /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParseTaxprofiler.nf --Taxprofiler_out ${TaxprofilerOutputFolder} --OutputDir ${ParsingOutputFolder} --TaxDump_gz ${taxdumpgz} --DepthTresh ${DepthTresh} --dbsheet ${DatabaseSheet} -c /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Mandalore_Taxprofiler.config

