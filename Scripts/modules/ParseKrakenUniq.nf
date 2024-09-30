process PARSEKRAKENUNIQ {
	publishDir "${params.OutputDir}/Taxprofiler_Parsed",mode: 'copy', overwrite: true
       	label "single_thread"
	
	module 'micromamba/1.4.2'

	input:
		val Taxprofiler_out
		val DepthTresh_c
		val DatabaseSheet_c
		
	output:
		path 'KrakenUniq', emit: KrakenUniq_outs

	script:
	if (params.IgnoreReadExtraction_krakenuniq)
	"""

	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	which ParseKrakenUniq.py

	ParseKrakenUniq.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c} --IgnoreReadExtraction

	"""		

	else
	"""
	
	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	ParseKrakenUniq.py --Taxprofiler_out ${Taxprofiler_out} --DepthTresh ${DepthTresh_c} --Db_sheet ${DatabaseSheet_c} 

	"""		


}

