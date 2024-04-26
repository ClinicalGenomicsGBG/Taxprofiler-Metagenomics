process PARSEDIAMOND {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "single_thread"

	module 'miniconda/4.14.0'

	input:
		val Taxprofiler_out
		val DepthTresh_c
                val DatabaseSheet_c
		val taxdumpgz_c

		
	output:
		path 'Diamond', emit: Diamond_outs

	script:
	if (params.IgnoreReadExtraction_Diamond)
	"""

	source activate TaxProfiler

	python /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseDiamond.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c} --taxdumpfile ${taxdumpgz_c} --IgnoreReadExtraction
	"""
	
	else 

	"""

	source activate TaxProfiler
	
	python /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParserScripts/ParseDiamond.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --taxdumpfile ${taxdumpgz_c} --Db_sheet ${DatabaseSheet_c}

	"""		

}

