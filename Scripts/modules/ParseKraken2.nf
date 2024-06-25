process PARSEKRAKEN2 {
	publishDir "${params.OutputDir}/Taxprofiler_Parsed", mode: 'copy', overwrite: true
       	label "single_thread"

	module 'miniconda/4.14.0'

	input:
		val Taxprofiler_out
		val DepthTresh_c
		val DatabaseSheet_c
		
	output:
		path 'Kraken2', emit: Kraken2_outs

	script:
	if (params.IgnoreReadExtraction_kraken2)

	"""
	
	source activate TaxProfiler

	python /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseKraken2.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c} --IgnoreReadExtraction

	"""		

	else

	"""
	
	source activate TaxProfiler
	
	python /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseKraken2.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c}

	"""
}

