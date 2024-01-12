# Taxprofiler-Metagenomics

Taxprofiler and Parser on medstore

# Background



# Instructions


## Installation


Load Taxprofiler environment using miniconda

```

module load miniconda/4.14.0

source activate TaxProfiler

```

Downlaod the newer version (tested on 1.1.2)


```

nf-core download taxprofiler -r 1.1.2

or 

nextflow pull nf-core/taxprofiler

```

## Workflow

1. Generate your $SampleSheet, comma separated containing the following metadata (SampleSheeet_Taxprofiler.csv)

Example below:

```

sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta
Sample1,Run1,ILLUMINA,/path/to/Sample1.R1.fastq.gz,/path/to/Sample1.R2.fastq.gz,

```


2. (*Optional*) Database is already set but if you want to change it generate your $DatabaseSheet, comma seperated containing the paths to the databases you want to include.

```
tool,db_name,db_params,db_path
kraken2,Kraken2_PlusPF_091023,--quick,/medstore/databases/taxprofiler/databases/kraken2/Kraken2_PlusPF_091023
diamond,Diamond_230321,,/medstore/databases/taxprofiler/databases/diamond/
krakenuniq,krakenuniq_MicrobialDB,,/medstore/databases/taxprofiler/databases/KrakenUniq_MicrobialDB
metaphlan,Metaphlan4,,/medstore/databases/taxprofiler/databases/Metaphlan4

```

3. Submit with qsub

(Make sure that the paths to the configs and database files are correct first) 


```

qsub /medstore/Development/Metagenomics/TaxProfiler/runTaxprofiler.sh <SampleSheet> <OutDir>

```

4. When Taxprofiler is done and you want to extract the reads and generate count tables. *Important* Need to include the entire path to the Taxprofiler output directory! 

```

qsub runParser_nextflow.sh <Full/Path/To/TaxprofilerOutPutDir> <ParserOutputDir>

```

## Setting up the server config

The config is tailored for mandalore but can be tweaked. This config was created to make sure we avoid the memory issues as metagenomics analyses are very memory heavy.  
We set the ```queueSize=20``` this means we are not allowing more then 20 processes at the time. For the process_high_memory and the process_high we are using ```MaxForks=1```. This means that we are only allowing one run at the time for these. This affect the ```bowtie2``` (for removing the human reads), ```kraken2``` and ```krakenuniq```. We are only allowing one process, this is a bottleneck. One could generate one withName rule for each but then you need to define to which node you want to queue towards which is not optimal as we dont know which nodes are free beforehand. You only want 1 of these heavy processes in one node otherwise you risk a crash. We set a special process for ```Diamond```. Diamond is ```process_medium```, we overwrite the ```process_medium``` for ```Diamond``` with ```withName```. Here we have ```maxForks=10``` running it on ```kuat```. This is for now, but did not want to risk the diamond processes in a node were a heavy process is running. Diamond is slow but to many processes at the same time will clogg the memory. 