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

#nextflow /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/main.nf --OutputDir test --SampleSheet $1 -c /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Mandalore_Taxprofiler_mod.config --IgnoreReadExtraction_krakenuniq --Metadata /medstore/Development/Metagenomics/TaxProfiler/InsilicoTester_SE/metadata.csv -resume 

nextflow /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/main.nf -N sanna.abrahamsson@gu.se --OutputDir test --SampleSheet $1 -c /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Mandalore_Taxprofiler_mod.config --IgnoreReadExtraction_krakenuniq --Webstore --Webstore_mail -with-report test/summary_report.html -resume 


