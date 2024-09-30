process PARSEEXCEL {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "single_thread" 

	module 'micromamba/1.4.2'

	input:
		val Diamond_outs
		val Kraken2_outs
		val KrakenUniq_outs
		
	output:
		path 'Out_ParsedExcel.xlsx', emit: ExcelOut

	script:
	"""
	
	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	WriteExcel.py --Tools ${Diamond_outs} ${Kraken2_outs} ${KrakenUniq_outs}

	"""		

}

