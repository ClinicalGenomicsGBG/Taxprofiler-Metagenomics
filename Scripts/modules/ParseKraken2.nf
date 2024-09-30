process PARSEKRAKEN2 {
	publishDir "${params.OutputDir}/Taxprofiler_Parsed", mode: 'copy', overwrite: true
       	label "single_thread"

	module 'micromamba/1.4.2'

	input:
		val Taxprofiler_out
		val DepthTresh_c
		val DatabaseSheet_c
		
	output:
		path 'Kraken2', emit: Kraken2_outs
		path 'Kraken2/*_CountsForplotting.txt', emit: Kraken2_outs_countfiles
		path 'Kraken2/Extras/*_SpeciesDomainLinkage.txt', emit: Kraken2_outs_linkagefiles

	script:
	if (params.IgnoreReadExtraction_kraken2)

	"""
	
	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	ParseKraken2.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c} --IgnoreReadExtraction

	"""		

	else

	"""
	
	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8
	
	ParseKraken2.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c}

	"""
}

