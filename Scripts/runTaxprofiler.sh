#!/bin/bash -l

#$ -N Taxprofiler_TestClinical
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

module load miniconda/4.14.0
source activate TaxProfiler
module load singularity/v3.5.2

set -x

SampleSheet=$1
outDir=$2

source /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Bash_Config_Taxprofiler.config

nextflow run nf-core/taxprofiler --perform_shortread_qc --shortread_qc_tool 'fastp' --perform_shortread_hostremoval --shortread_hostremoval_index $hostindex --hostremoval_reference $hostfasta --save_preprocessed_reads --perform_shortread_complexityfilter --shortread_complexityfilter_tool 'bbduk' -profile singularity --input $SampleSheet --databases $DatabaseSheet --outdir $outDir -r 1.1.5 -c $MandaloreConfig --run_kraken2 --kraken2_save_readclassifications --kraken2_save_reads --run_krakenuniq --krakenuniq_save_reads --krakenuniq_save_readclassifications --run_diamond --diamond_output_format tsv --run_krona --save_hostremoval_unmapped 


