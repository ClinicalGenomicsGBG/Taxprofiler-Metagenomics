process COLLATE {
	publishDir params.OutputDir,mode: 'copy', overwrite: true
       	label "single_thread"

	module 'miniconda/4.14.0'

	input:
		val ExcelOut
		val comparisons

	script:

	if (params.Webstore_mail)

	"""
	
	source activate TaxProfiler

	Collate.py --FolderToMove "$launchDir/${params.OutputDir}" --mail

	"""

	else

	"""

	Collate.py --FolderToMove "$launchDir/${params.OutputDir}"	
	

	"""

}


