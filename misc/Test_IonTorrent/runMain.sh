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


SampleSheet=/medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/misc/Test_IonTorrent/SampleSheet.csv
outDir=TaxProfiler_out
Metadata=/medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/misc/Test_IonTorrent/Metadata.csv
Config=/medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Mandalore_Taxprofiler_mod.config


nextflow run /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/main.nf --OutputDir $outDir --SampleSheet $SampleSheet -c $Config --IgnoreReadExtraction_krakenuniq --Metadata $Metadata



