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

Downlaod the newer version (tested on 1.0.1 for now)


```

nf-core download taxprofiler -r 1.0.1

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
kraken2,krak_230419,--quick,/medstore/databases/taxprofiler/databases/kraken2/kraken2_230419/AVB
bracken,bracken_150,;-r 150,/medstore/databases/taxprofiler/databases/kraken2/kraken2_230419/AVB
kaiju,Kaiju_140623,,/medstore/databases/taxprofiler/databases/Kaiju_140623/
diamond,Diamond_230321,,/medstore/databases/taxprofiler/databases/diamond/
krakenuniq,krakenuniq_MicrobialDB,,/medstore/databases/taxprofiler/databases/KrakenUniq_MicrobialDB
metaphlan3,Metaphlan3_manualDownload,,/medstore/databases/taxprofiler/databases/Metaphlan3_manualDownload

```

3. Submit with qsub

(Make sure that the paths to the configs and database files are correct first) 


```

qsub /medstore/Development/Metagenomics/TaxProfiler/runTaxprofiler.sh <SampleSheet> <OutDir>

```
