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
	     --Metadata				Path to metadata file required for comparison patient vs control. csv file with metadata, Sample (samplename), Type (DNA or RNA), Group (P (patient) or C (control)
	     --Webstore				Use flag if the results should be transferred to webstore
	     --Webstore_mail			Use together with webstore flag if a mail should be sent after completed transfer
             --help   This usage statement


    """
}

// show help message
params.help=false
params.IgnoreReadExtraction_Diamond=false
params.IgnoreReadExtraction_kraken2=false
params.IgnoreReadExtraction_krakenuniq=false
params.Metadata=false
params.Webstore=false
params.Webstore_mail=false
params.dummycomparison='did not run'

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

if (params.Metadata){
   Metadata_c=Channel.from(params.Metadata)
}

include {TAXPROFILER}	from	'./modules/Taxprofiler.nf'
include {PARSEKRAKEN2}	from	'./modules/ParseKraken2.nf'
include {PARSEDIAMOND}	from	'./modules/ParseDiamond.nf'
include {PARSEKRAKENUNIQ}	from	'./modules/ParseKrakenUniq.nf'
include {PARSEEXCEL}	from	'./modules/ParseExcel.nf'
include {COMPARISONS}	from	'./modules/Comparisons.nf'
include {COLLATE}	from	'./modules/CollateResults.nf'

//------------------------------WorkFlow-Parameters--------------------------------------------------


workflow {
	 // Workflow
	 TAXPROFILER(hostindex_c, hostfasta_c, SampleSheet_c, DatabaseSheet_c, MandaloreConfig_c)
	 // ParseKraken2
	 PARSEKRAKEN2(TAXPROFILER.out.Taxprofiler_out, DepthTresh_c, DatabaseSheet_c)
	 PARSEDIAMOND(TAXPROFILER.out.Taxprofiler_out, DepthTresh_c, DatabaseSheet_c, taxdumpgz_c)
	 PARSEKRAKENUNIQ(TAXPROFILER.out.Taxprofiler_out, DepthTresh_c, DatabaseSheet_c)
	 PARSEEXCEL(PARSEDIAMOND.out.Diamond_outs, PARSEKRAKEN2.out.Kraken2_outs, PARSEKRAKENUNIQ.out.KrakenUniq_outs)
	 // Optional to add comparison
	 if (params.Metadata){
	    COMPARISONS(PARSEDIAMOND.out.Diamond_outs, PARSEKRAKEN2.out.Kraken2_outs, PARSEKRAKENUNIQ.out.KrakenUniq_outs, Metadata_c)
	    }
	  // Copy to webstore without the comparison done
	 if (params.Webstore){
	    if (params.Metadata){
	    // Collate after the comparison is done
	    COLLATE(PARSEEXCEL.out.ExcelOut, COMPARISONS.out.comparisons)
	    }
	    else {
	    COLLATE(PARSEEXCEL.out.ExcelOut,params.dummycomparison)
	    }}}
