process PARSEDIAMOND {
	publishDir "${params.OutputDir}/Taxprofiler_Parsed",mode: 'copy', overwrite: true
       	label "single_thread"

	module 'micromamba/1.4.2'

	input:
		val Taxprofiler_out
		val DepthTresh_c
                val DatabaseSheet_c
		val taxdumpgz_c

		
	output:
		path 'Diamond', emit: Diamond_outs
		path 'Diamond/*_CountsForplotting.txt', emit: Diamond_outs_countfiles
                path 'Diamond/Extras/*_SpeciesDomainLinkage.txt', emit: Diamond_outs_linkagefiles


	script:
	if (params.IgnoreReadExtraction_Diamond)
	"""

	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	ParseDiamond.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c} --taxdumpfile ${taxdumpgz_c} --IgnoreReadExtraction
	"""
	
	else 

	"""

	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8
	
	ParseDiamond.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --taxdumpfile ${taxdumpgz_c} --Db_sheet ${DatabaseSheet_c}

	"""		

}

