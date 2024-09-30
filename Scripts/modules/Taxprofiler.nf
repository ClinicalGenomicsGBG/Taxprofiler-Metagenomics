process TAXPROFILER {
        publishDir params.OutputDir, mode: 'copy', overwrite: true
        label "single_thread"

	module 'micromamba/1.4.2'

        input:
	
		val hostindex_c
		val hostfasta_c
		val SampleSheet_c
		val DatabaseSheet_c
		val MandaloreConfig_c

        output:
                path 'TaxProfiler_out', emit: Taxprofiler_out

        script:
        """
	
	micromamba activate /medstore/projects/P23-015/Intermediate/MicroMambaEnvs/TaxProfiler_1.1.8

	nextflow run nf-core/taxprofiler --perform_shortread_qc --shortread_qc_tool 'fastp' --perform_shortread_hostremoval --shortread_hostremoval_index ${hostindex_c} --hostremoval_reference ${hostfasta_c} --save_preprocessed_reads --perform_shortread_complexityfilter --shortread_complexityfilter_tool 'bbduk' -profile singularity --input ${SampleSheet_c} --databases ${DatabaseSheet_c} --outdir TaxProfiler_out -r 1.1.8 -c ${MandaloreConfig_c} --run_kraken2 --kraken2_save_readclassifications --kraken2_save_reads --run_krakenuniq --krakenuniq_save_reads --krakenuniq_save_readclassifications --run_diamond --diamond_output_format tsv --run_krona --save_hostremoval_unmapped -resume
	
        """

}
