#!/bin/bash -l

#$ -S /bin/bash
#$ -N TaxprofilerRun
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

module load miniconda/4.14.0
source activate TaxProfiler

set -x

SampleSheet=/path/to/samplesheet.csv # Change here
outDir=TaxProfiler_out
Config=/medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Mandalore_Taxprofiler_mod.config
Metadata=/path/to/metadata.csv

#nextflow /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/main.nf -N sanna.abrahamsson@gu.se --OutputDir test --SampleSheet $1 -c /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Mandalore_Taxprofiler_mod.config --IgnoreReadExtraction_krakenuniq --Webstore --Webstore_mail -with-report test/summary_report.html -resume 

nextflow /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/main.nf -N sanna.abrahamsson@gu.se --OutputDir $outDir --SampleSheet $SampleSheet -c $Config --IgnoreReadExtraction_krakenuniq --Metadata $Metadata

