# Taxprofiler-Metagenomics

Taxprofiler and Parser on medstore

# Background



# Instructions

TaxProfiler is already installed in the conda environment TaxProfiler

```

module load miniconda/4.14.0

source activate TaxProfiler

```


## Installation (Ignore if already installed and updated!)

Load Taxprofiler environment using miniconda

```

module load miniconda/4.14.0

source activate TaxProfiler

```

Download the newer version (tested on 1.1.5)


```

nf-core download taxprofiler -r 1.1.5

or 

nextflow pull nf-core/taxprofiler

```

## Workflow

1. Generate your $SampleSheet, comma separated containing the following metadata (SampleSheeet_Taxprofiler.csv)

Example below (illimuna):

```

sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
Sample1,RNA,ILLUMINA,/path/to/Sample1.R1.fastq.gz,/path/to/Sample1.R2.fastq.gz,
Sample2,RNA,ILLUMINA,/path/to/Sample2.R1.fastq.gz,/path/to/Sample2.R2.fastq.gz,
Sample3,DNA,ILLUMINA,/path/to/Sample3.R1.fastq.gz,/path/to/Sample3.R2.fastq.gz,
Sample4,DNA,ILLUMINA,/path/to/Sample4.R1.fastq.gz,/path/to/Sample4.R2.fastq.gz,


```

Example below (iontorrent):


```

sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
Sample1,RNA,ION_TORRENT,/path/to/Sample1.fastq.gz,,
Sample2,RNA,ION_TORRENT,/path/to/Sample2.fastq.gz,,
Sample3,DNA,ION_TORRENT,/path/to/Sample3.fastq.gz,,
Sample4,DNA,ION_TORRENT,/path/to/Sample4.fastq.gz,,

```

2. (*Optional*) Database is already set in the config file but if you want to change it generate your $DatabaseSheet, comma seperated containing the paths to the databases you want to include.

```
tool,db_name,db_params,db_path
kraken2,Kraken2_PlusPF_091023,--quick,/medstore/databases/taxprofiler/databases/kraken2/Kraken2_PlusPF_091023
diamond,Diamond_230321,,/medstore/databases/taxprofiler/databases/diamond/
krakenuniq,krakenuniq_MicrobialDB,,/medstore/databases/taxprofiler/databases/KrakenUniq_MicrobialDB
metaphlan,Metaphlan4,,/medstore/databases/taxprofiler/databases/Metaphlan4

```

3. (*Optional*) When comparisons are necessary

When there are clinical samples the hospital wants to run the patient against a control (this one usually contains a high amount of EBV). To be able to run the comparison module for the parser you need to add a metadata file. This file need to have have the following format.

```
Sample,Type,Group
Sample1,RNA,P
Sample2,RNA,C
Sample3,DNA,P
Sample4,DNA,C

```

* Sample name as the fastq file, remove the *fastq.gz extension.
* The type is either RNA or DNA.
* The group where P stands for Patient and C stands for Control.

The comparisons are made pairwise, All RNA P will be tested against all RNA C, one by one. No replicates with mean calculations. 

4. Submit with qsub

The main script runs taxprofiler (with the classifiers Diamond, Kraken2 and Krakenuniq), When this is done it runs the parser. If the user included the Metadata with the metadataflag you run the comparisons.

```

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

SampleSheet=/path/to/samplesheet.csv # Change here
outDir=TaxProfiler_out

Config=/medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/configs/Mandalore_Taxprofiler_mod.config

nextflow run /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/main.nf --OutputDir $outDir --SampleSheet $SampleSheet -c $Config --IgnoreReadExtraction_krakenuniq --Metadata /path/to/metadata.csv # change here


```

run it with following: 

```

qsub runMain.sh

```


## About the outputs

The counts are extracted at species level, the counts are normalized by using a pseudocount, dividing by the total sum in the samplem take the fraction times 100. This will result in percentage of the total. 

* TaxProfiler_out: Main folder output from taxprofiler
* Taxprofiler_Parsed: Intermediates like count files and reads from the parser.  
* excel sheet: Out_ParsedExcel.xlsx
* Comparisons: (if user supplied the Metadata flag) Comparisons.xlsx

5. Copy output to webstore

Copy module to webstore is still ongoing, for now copy the files manually to ```/webstore/clinical/development/Taxprofiler```.
For now only move the TaxProfiler_out folder containing the raw outputs from taxprofiler and the excel files Comparisons.xlsx and Out_ParsedExcel.xlsx to a new folder with the days date.


## Setting up the server config

The config is tailored for mandalore but can be tweaked. This config was created to make sure we avoid the memory issues as metagenomics analyses are very memory heavy.  
We set the ```queueSize=20``` this means we are not allowing more then 20 processes at the time. For the process_high_memory and the process_high we are using ```MaxForks=1```. This means that we are only allowing one run at the time for these. This affect the ```bowtie2``` (for removing the human reads), ```kraken2``` and ```krakenuniq```. We are only allowing one process, this is a bottleneck. One could generate one withName rule for each but then you need to define to which node you want to queue towards which is not optimal as we dont know which nodes are free beforehand. You only want 1 of these heavy processes in one node otherwise you risk a crash. We set a special process for ```Diamond```. Diamond is ```process_medium```, we overwrite the ```process_medium``` for ```Diamond``` with ```withName```. Here we have ```maxForks=10``` running it on ```kuat```. This is for now, but did not want to risk the diamond processes in a node were a heavy process is running. Diamond is slow but to many processes at the same time will clogg the memory.

# Updating the databases

Follow these instructions when it is time to update taxprofiler databases

We are running the old databases now so ignore this, but the newer described here will be implemented during validation

## Diamond

Diamond requires aa sequenses in fasta format, we will be using the non redundant database from refseq for Bacteria, Virus, archaea and fungi. We are can download what we need with the following script

current version: refseq_223 (15 March 2024)

```

module load miniconda/4.14.0

source activate TaxProfiler

python Download_fastas_Refseq_Diamond.py --out PathToNewDiamondDB 

```

The script townloads all nun redundant fastas (aa) for bacteria, virus, archaea and funig, it also downlaods the taxdumpfile ```taxdump.tar.gz``` and the protein2 accession file ```prot.accession2taxid.FULL.gz``` needed for the taxonomy. 

Extract the taxdump file

```

tar -xvf taxdump.tar.gz --directory taxdump

```

There is a lot of things that can go wrong when you download so many individual files, You can check that they are complete, If something is wrong with one of the files the database creation just interupts without an error message. If you run make database with a concat file you will get an inflation error if one of the sequence files are corrupt. 

Double check that all the faa files were downloaded correctly

```

for i in *faa.gz; do gzip -tv $i; done

```

Then we create the diamond database

```
#!/bin/bash -l

#$ -N DiamondDB
#$ -j y
#$ -cwd
#$ -pe mpi 6
#$ -q development.q

module load miniconda/4.14.0

source activate diamond

set -x

zcat fastas/*.faa.gz | diamond makedb --taxonmap prot.accession2taxid.FULL.gz --taxonnames taxdump/names.dmp --taxonnodes taxdump/nodes.dmp -d refseq_223_Bacteria_Virus_Archaea_Fungi_diamond -v --log --threads 6


```


```

## Kraken2

Create a new folder for your download in folder ```/medstore/databases/taxprofiler/databases/kraken2/```

You can download kraken2 updated databases from [here](https://benlangmead.github.io/aws-indexes/k2)

current version: PlusPF, k2_pluspf_20240112

**obs: Change file name to the new version you are downloading**

```

#!/bin/bash -l

#$ -N download
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q


set -x

wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240112.tar.gz

tar -zxvf k2_pluspf_20240112.tar.gz

```


## Kraken2Uniq

Create a new folder for your download in folder ```/medstore/databases/taxprofiler/databases/Krakenuniq```

You can download kraken2 updated databases from [here](https://benlangmead.github.io/aws-indexes/k2),

current version: MicrobialDB, kuniq_microbialdb_minus_kdb

**obs: Change file name to the new version you are downloading**

```

#!/bin/bash -l

#$ -N download
#$ -j y
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

set -x

wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/kuniq_microbialdb_minus_kdb.20230808.tgz

tar -zxvf kuniq_microbialdb_minus_kdb.20230808.tgz

wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/database.kdb

```


set the new paths to the database sheet

