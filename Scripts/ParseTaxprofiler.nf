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

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseKraken2.py --Taxprofiler_out ${Taxprofiler_out_p} --DepthTresh ${DepthTresh_p} 

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

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseKrakenUniq.py --Taxprofiler_out ${Taxprofiler_out_p} --DepthTresh ${DepthTresh_p} 

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

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseDiamond.py --Taxprofiler_out ${Taxprofiler_out_p} --DepthTresh ${DepthTresh_p} 

	"""		

}


// Ganon is tester now until updating of the pipeline, also think it is a bit memory heavy, need to think about how it runs! 
process ParseGanon {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "process_single" 

	input:
		val Taxprofiler_out_p
		val DepthTresh_p
		
	output:
		path 'Ganon', emit: Ganon_outs

	script:
	"""

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseGanon.py --Taxprofiler_out ${Taxprofiler_out_p} --DepthTresh ${DepthTresh_p} --taxdumpfile /medstore/databases/taxprofiler/databases/taxpasta_taxonomy/taxdump

	"""		

}


process ParseExcel {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "process_single" 

	input:
		val Diamond_outs
		val Kraken2_outs
		val KrakenUniq_outs
		
	output:
		path 'Diamond', emit: Diamond_outs

	script:
	"""

	/apps/bio/software/anaconda2/envs/TaxProfiler/bin/python /medstore/Development/Metagenomics/TaxProfiler/TestSamples_FromMicro/TaxprofilerRun/Scripts/Taxprofiler-Metagenomics/Scripts/ParserScripts/WriteExcel.py --Tools ${Diamond_outs} ${Kraken2_outs} ${Metaphlan4_outs} ${Ganon_outs}

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
	//ParseMetaPhlan4(Taxprofiler_out_p, DepthTresh_p)
	ParseKraken2(Taxprofiler_out_p, DepthTresh_p)
	// Turning of krakenuniq for now, it is to memory heavy 
	ParseKrakenUniq(Taxprofiler_out_p, DepthTresh_p) 
	ParseDiamond(Taxprofiler_out_p, DepthTresh_p)
	// Ganon is just for tester here as well, we might need to remove it also! 
	//ParseGanon(Taxprofiler_out_p, DepthTresh_p)
	//ParseExcel(ParseDiamond.out.Diamond_outs, ParseKraken2.out.Kraken2_outs)
}