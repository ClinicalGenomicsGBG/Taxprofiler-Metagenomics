#!/usr/bin/py

import argparse
import subprocess
import sys
import requests
import os

def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Downloads the protein fasta files from refseq togehter with current taxdump. It extracts data from virus, archaea, fungi and bacteria using non redundant database (wb)')
    parser.add_argument("--out", dest = 'OutFolder', required=True, help ="Path to Output Folder (required)")
    arguments = parser.parse_args(argv)
    return arguments

def download(OutFolder):
    if not os.path.exists(OutFolder):
        os.makedirs(OutFolder)
    if not os.path.exists(f'{OutFolder}/fastas'):
        os.makedirs(f'{OutFolder}/fastas')


    # because there might be errors when downloading and concatenating many files i take them one by one so if something goes wrong we can download that specific sample again! 
        
    # Bacteria and archaea, non redundant contains wp in filename
    Taxas=["bacteria", "archaea"]
    PathToRefseq="https://ftp.ncbi.nlm.nih.gov/refseq/release"
    for taxa in Taxas:
        print(f'Downloading {taxa}')
        r=requests.get(f'{PathToRefseq}/{taxa}')
        data=r.text
        for line in data.split("\n"):
            if ".protein.faa.gz" in line:
                if ".wp_protein." in line:
                    file_taxa=line.split("\"")[1]
                    fullpath=f'{PathToRefseq}/{taxa}/{file_taxa}'
                    with open(f'{OutFolder}/fastas/{file_taxa}', "wb") as o:
                        todownload=requests.get(fullpath)
                        o.write(todownload.content)
    
    # Virus and Fungi ands in a slighty different manner compared to bacteria, that is why we need two loops
    Taxas=["viral","fungi"]         
    for taxa in Taxas:
        print(f'Downloading {taxa}')
        r=requests.get(f'{PathToRefseq}/{taxa}')
        data=r.text
        for line in data.split("\n"):
            if ".protein.faa.gz" in line:
                file_taxa=line.split("\"")[1]
                fullpath=f'{PathToRefseq}/{taxa}/{file_taxa}'
                with open(f'{OutFolder}/fastas/{file_taxa}', "wb") as o:
                    todownload=requests.get(fullpath)
                    o.write(todownload.content)

                    
    # Download the taxdump
    print("Downloading taxdump")
    PathToTaxdump="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
    taxdumpfile="taxdump.tar.gz"
    fullpath=f'{PathToTaxdump}/{taxdumpfile}'
    with open(f'{OutFolder}/{taxdumpfile}', "wb") as o:
        todownload=requests.get(fullpath)
        o.write(todownload.content)

    # Download prot to accession
    print("Downloading accession")
    PathToaccession="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
    accessionfile="prot.accession2taxid.FULL.gz"
    fullpath=f'{PathToaccession}/{accessionfile}'
    with open(f'{OutFolder}/{accessionfile}', "wb") as o:
        todownload=requests.get(fullpath)
        o.write(todownload.content)

    

    
    
        
        
def main(OutFolder):
    download(OutFolder)
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.OutFolder)


    
