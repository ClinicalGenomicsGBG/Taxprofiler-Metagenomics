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
    	     --Taxprofiler_out	Output directory of Taxprofiler
	     --DepthTresh	Depth Treshold to count it as detected
	     --OutputDir	OutPut folder to store the results in 
	     -c			configuration file for nextflow (Works for Mandalores nodes)
    Optional arguments:
    	     --help   This usage statement

    """
}

// show help message
params.help=false

if (params.help){
   helpMessage()
   exit 0
}


process ParseMetaPhlan4 {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "process_single" 

	input:
		val Taxprofiler_out_p
		val DepthTresh_p
		
	output:
		path 'Metaphlan4', emit: Metaphlan4_outs

	script:
	"""

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseMetaphlan4.py --Taxprofiler_out ${Taxprofiler_out_p} --DepthTresh ${DepthTresh_p} --taxdumpfile /medstore/databases/taxprofiler/databases/taxpasta_taxonomy/taxdump

	"""		

}




process ParseKraken2 {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "process_single" 

	input:
		val Taxprofiler_out_p
		val DepthTresh_p
		
	output:
		path 'Kraken2', emit: Kraken2_outs

	script:
	"""

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseKraken2.py --Taxprofiler_out ${Taxprofiler_out_p} --DepthTresh ${DepthTresh_p} --taxdumpfile /medstore/databases/taxprofiler/databases/taxpasta_taxonomy/taxdump

	"""		

}


process ParseKrakenUniq {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "process_single" 

	input:
		val Taxprofiler_out_p
		val DepthTresh_p
		
	output:
		path 'KrakenUniq', emit: KrakenUniq_outs

	script:
	"""

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseKrakenUniq.py --Taxprofiler_out ${Taxprofiler_out_p} --DepthTresh ${DepthTresh_p} --taxdumpfile /medstore/databases/taxprofiler/databases/taxpasta_taxonomy/taxdump

	"""		

}



process ParseDiamond {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "process_single" 

	input:
		val Taxprofiler_out_p
		val DepthTresh_p
		
	output:
		path 'Diamond', emit: Diamond_outs

	script:
	"""

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseDiamond.py --Taxprofiler_out ${Taxprofiler_out_p} --DepthTresh ${DepthTresh_p} --taxdumpfile /medstore/databases/taxprofiler/databases/taxpasta_taxonomy/taxdump

	"""		

}





workflow {
	//--- Parameters From UserInput --- 
	params.Taxprofiler_out=''
	Taxprofiler_out_p=Channel.from(params.Taxprofiler_out)
	params.DepthTresh=''
	DepthTresh_p=Channel.from(params.DepthTresh)
	params.OutputDir=''

	//--- Workflow --- 
	ParseMetaPhlan4(Taxprofiler_out_p, DepthTresh_p)
	ParseKraken2(Taxprofiler_out_p, DepthTresh_p)
	ParseKrakenUniq(Taxprofiler_out_p, DepthTresh_p)
	ParseDiamond(Taxprofiler_out_p, DepthTresh_p)
}