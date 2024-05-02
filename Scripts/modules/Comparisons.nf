process COMPARISONS {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "single_thread"

	module 'miniconda/4.14.0'

	input:
		val Diamond_outs
		val Kraken2_outs
		val KrakenUniq_outs
		val Metadata_c
		
	output:
		path 'Comparisons.xlsx', emit: comparisons

	script:
	
	"""
	
	source activate TaxProfiler

	python /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParserScripts/runComparison.py  --Parsed_Taxprofiler_out "$launchDir/${params.OutputDir}/Taxprofiler_Parsed" --Metadata ${Metadata_c}

	"""		
}


