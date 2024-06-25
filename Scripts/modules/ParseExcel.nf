process PARSEEXCEL {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "single_thread" 

	module 'miniconda/4.14.0'

	input:
		val Diamond_outs
		val Kraken2_outs
		val KrakenUniq_outs
		
	output:
		path 'Out_ParsedExcel.xlsx', emit: ExcelOut

	script:
	"""
	
	source activate TaxProfiler

	python /medstore/Development/Metagenomics/TaxProfiler/Taxprofiler-Metagenomics/Scripts/ParserScripts/WriteExcel.py --Tools ${Diamond_outs} ${Kraken2_outs} ${KrakenUniq_outs}

	"""		

}

