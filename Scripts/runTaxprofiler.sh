#!/bin/bash -l

#$ -N Taxprofiler_TestClinical
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q


module load miniconda/4.14.0
source activate TaxProfiler
module load singularity/v3.5.2

# We need to remove the human reads first as well, Keeping them when we are not mapping towards it just increases the alignment time for all steps! 

SampleSheet=$1
outDir=$2
#outDir=ClinicalSamplesTestTaxprofiler

DatabaseSheet=/medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/database_sheet.csv
taxdumpdir=/medstore/databases/taxprofiler/databases/taxpasta_taxonomy/taxdump
hostindex=/medstore/databases/taxprofiler/databases/GRCH38/
hostfasta=/medstore/databases/taxprofiler/databases/GRCH38/Homo_sapiens.GRCh38.dna.toplevel.canonical.fa
MandaloreConfig=/medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/database_sheet.csv


# Skip Kaiju, try extracting as much reads as possible

# Which tools can save reads?
# Centrifuge,
# Diamond, but then there will not be a report
# kraken2
# krakenuniq


echo "

nextflow run nf-core/taxprofiler --perform_shortread_qc --shortread_qc_tool 'fastp' --perform_shortread_hostremoval --shortread_hostremoval_index $hostindex --hostremoval_reference $hostfasta --save_preprocessed_reads --perform_shortread_complexityfilter --shortread_complexityfilter_tool 'bbduk' -profile singularity --input $SampleSheet --databases $DatabaseSheet --outdir $outDir -r 1.0.1 -c $MandaloreConfig --run_kraken2 --kraken2_save_readclassification --run_diamond --run_bracken --run_krona --run_krakenuniq --krakenuniq_save_readclassifications --run_metaphlan3 --save_hostremoval_unmapped -resume

"

nextflow run nf-core/taxprofiler --perform_shortread_qc --shortread_qc_tool 'fastp' --perform_shortread_hostremoval --shortread_hostremoval_index $hostindex --hostremoval_reference $hostfasta --save_preprocessed_reads --perform_shortread_complexityfilter --shortread_complexityfilter_tool 'bbduk' -profile singularity --input $SampleSheet --databases $DatabaseSheet --outdir $outDir -r 1.0.1 -c $MandaloreConfig --run_kraken2 --kraken2_save_readclassification --run_diamond --run_bracken --run_krona --run_krakenuniq --krakenuniq_save_readclassifications --run_metaphlan3 --save_hostremoval_unmapped -resume


echo "Fin"
