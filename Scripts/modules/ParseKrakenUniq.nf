process PARSEKRAKENUNIQ {
	publishDir "${params.OutputDir}/Taxprofiler_Parsed",mode: 'copy', overwrite: true
       	label "single_thread"
	
	module 'miniconda/4.14.0'

	input:
		val Taxprofiler_out
		val DepthTresh_c
		val DatabaseSheet_c
		
	output:
		path 'KrakenUniq', emit: KrakenUniq_outs

	script:
	if (params.IgnoreReadExtraction_krakenuniq)
	"""

	source activate TaxProfiler

	python /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseKrakenUniq.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c} --IgnoreReadExtraction

	"""		

	else
	"""
	
	source activate TaxProfiler

	python /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseKrakenUniq.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c} 

	"""		


}

