process COMPARISONS {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "single_thread"

	module 'micromamba/1.4.2'

	input:
		val Diamond_outs
		val Kraken2_outs
		val KrakenUniq_outs
		val Metadata_c
		
	output:
		path 'Comparisons.xlsx', emit: comparisons

	script:
	
	"""
	
	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	runComparison.py  --Tools ${Diamond_outs} ${Kraken2_outs} ${KrakenUniq_outs} --Metadata ${Metadata_c}

	"""		
}


