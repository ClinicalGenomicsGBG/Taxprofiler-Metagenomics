#!/usr/bin/py

import setuptools

setuptools.setup(
    name='Taxprofiler-medstore',  
    version='0.1',
    scripts=['Scripts/runParser_nextflow.sh','Scripts/ParseTaxprofiler.nf','Scripts/runTaxprofiler.sh', 'Scripts/ParserScripts/ParseDiamond.py',  'Scripts/ParserScripts/ParseGanon.py', 'Scripts/ParserScripts/ParseKraken2.py', 'Scripts/ParserScripts/ParseKrakenUniq.py', 'Scripts/ParserScripts/ParseMetaphlan4.py'],
    author="Sanna Abrahamsson",
    author_email="sanna.abrahamsson@gu.se",
    description="Installation for running Taxprofiler in Mandalore",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ClinicalGenomicsGBG/Taxprofiler-Metagenomics.git",
    packages=setuptools.find_packages(),
    package_data={'Taxprofiler-medstore': ['README.md', "database_sheet.csv","SampleSheeet_Taxprofiler.csv"]
                   },
    data_files= [('configs',['configs/Mandalore_Taxprofiler.config','configs/Bash_Config_Taxprofiler.config','configs/base.config'])],
    include_package_data=True,
    classifiers=[
         "Programming Language :: Python :: 3", 
         
     ],
 )
