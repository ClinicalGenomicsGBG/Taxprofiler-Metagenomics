process COLLATE {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "single_thread"

	module 'micromamba/1.4.2'

	input:
		val ExcelOut
		val comparisons

	script:

	if (params.Webstore_mail)

	"""
	
	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	Collate.py --FolderToMove "$launchDir/${params.OutputDir}" --mail

	"""

	else

	"""

	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	Collate.py --FolderToMove "$launchDir/${params.OutputDir}"	
	

	"""

}


