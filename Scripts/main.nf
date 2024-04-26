#!/usr/bin/env nextflow

//nextflow.enable.dsl = 2

// Help Section


def helpMessage(){
    log.info """
    Usage:
    A Parser for taxprofiler, generates tables, plots etc.

    The typical command for running the ASO pipeline is as follows:
        nextflow run ParseTaxprofiler.nf

    Mandatory arguments:
             --OutputDir	OutPut folder to store the results in
	     --SampleSheet	SampleSheet for TaxProfiler
             -c          	configuration file for nextflow (Works for Mandalores nodes)

    Optional arguments:
             --IgnoreReadExtraction_Diamond	Ignore extracting the reads from diamond classification
             --IgnoreReadExtraction_kraken2	Ignore extracting the reads from kraken2 classification
             --IgnoreReadExtraction_krakenuniq	Ignore extracting the reads from krakenuniq, this is pending to be fixed
             --help   This usage statement


    """
}

// show help message
params.help=false
params.IgnoreReadExtraction_Diamond=false
params.IgnoreReadExtraction_kraken2=false
params.IgnoreReadExtraction_krakenuniq=false


if (params.help){
   helpMessage()
   exit 0
}


//------------------------------Read-Parameters--------------------------------------------------

// Params from user
params.OutputDir=''
params.SampleSheet=''
SampleSheet_c=Channel.from(params.SampleSheet)
// Params from config
DatabaseSheet_c=Channel.from(params.DatabaseSheet)
MandaloreConfig_c=Channel.from(params.MandaloreConfig)
taxdumpgz_c=Channel.from(params.taxdumpgz)
hostindex_c=Channel.from(params.hostindex)
hostfasta_c=Channel.from(params.hostfasta)
DepthTresh_c=Channel.from(params.DepthTresh)

include {TAXPROFILER}	from	'./modules/Taxprofiler.nf'
include {PARSEKRAKEN2}	from	'./modules/ParseKraken2.nf'
include {PARSEDIAMOND}	from	'./modules/ParseDiamond.nf'
include {PARSEKRAKENUNIQ}	from	'./modules/ParseKrakenUniq.nf'
include {PARSEEXCEL}	from	'./modules/ParseExcel.nf'

//------------------------------WorkFlow-Parameters--------------------------------------------------



workflow {

	 // Workflow
	 TAXPROFILER(hostindex_c, hostfasta_c, SampleSheet_c, DatabaseSheet_c, MandaloreConfig_c)
	 // ParseKraken2
	 PARSEKRAKEN2(TAXPROFILER.out.Taxprofiler_out, DepthTresh_c, DatabaseSheet_c)
	 PARSEDIAMOND(TAXPROFILER.out.Taxprofiler_out, DepthTresh_c, DatabaseSheet_c, taxdumpgz_c)
	 PARSEKRAKENUNIQ(TAXPROFILER.out.Taxprofiler_out, DepthTresh_c, DatabaseSheet_c)
	 PARSEEXCEL(PARSEDIAMOND.out.Diamond_outs, PARSEKRAKEN2.out.Kraken2_outs, PARSEKRAKENUNIQ.out.KrakenUniq_outs)

}


