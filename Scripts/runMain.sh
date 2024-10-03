#!/bin/bash -l
#$ -S /bin/bash
#$ -N TaxprofilerRun
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

module load micromamba/1.4.2
micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

# Setting up log

_timestamp_dir=`date '+%Y-%m-%d-%H-%M'`-`uuidgen -t`
_base_log_dir="/medstore/logs/pipeline_logfiles/TaxProfiler"
logpath=${_base_log_dir}/${_timestamp_dir}.log

# Setting up output names

Timestamp=$(date +%y%m%d-%H%M%S)

SampleSheet=$1 # Full path to sample sheet
Runname=$2 # Full name to run

outDir=${Runname}_TaxProfiler_out_${Timestamp}
Config=/medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Mandalore_Taxprofiler_mod_updatedDB.config

nextflow /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/main.nf --OutputDir $outDir --SampleSheet $SampleSheet -c $Config &>> ${logpath}

