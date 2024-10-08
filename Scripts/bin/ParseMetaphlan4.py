#!/usr/bin/env python


import argparse
from Bio import SeqIO
from collections import Counter
import glob
import gzip
import logging
import os
import pandas as pd
from plotnine import ggplot, aes, geom_bar, theme, element_text, xlab, ggtitle
import psutil
import subprocess
import sys
import time
import xlsxwriter


def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes the output from taxprofiler and parses it')
    parser.add_argument("--Taxprofiler_out", dest = 'taxprofdict', required=True, help ="Full path to Output folder from taxprofiler (required)")
    parser.add_argument("--taxdumpfile", dest = 'taxdump', required=True, help ="Path to taxdump (required)")
    parser.add_argument("--DepthTresh", dest = 'dptresh',default=10, type=int, help ="Minimum depth required to be reported")
    arguments = parser.parse_args(argv)
    return arguments


def ParseMetaphlan4(taxprofdict, taxdump,dptresh):
    """
    """

    
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="metaphlan"
    for i in subfolders:
        subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
        for k in subfolders_2:
            if "Metaphlan4" in k:
                try:
                    os.mkdir("Metaphlan4")
                except FileExistsError:
                    pass
                # Taxpasta
                metaphlandb=k.split("/")[-1]
                outTaxpasta= "Metaphlan4/Metaphlan4_" + metaphlandb + ".tsv"
                command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p metaphlan -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*_profile.txt" %(outTaxpasta, taxdump, k)
                subprocess.call(command, shell=True)
                # Printing required files for plotting
                Sampleorder=[]
                with open(outTaxpasta, "r") as taxpasta:
                    header = taxpasta.readline().split("\t")
                    for i in header[4:]:                        
                        samplename=i.split(metaphlandb+".metaphlan_profile")[0]
                        Sampleorder.append([i,samplename,header.index(i)])
                for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                    samplename=i[1]
                    samplename="_".join(samplename.rsplit("_pe_", 1)) # Remove the PE that Metaphlan adds to the name
                    sampleindex=i[2]
                    outforplot="Metaphlan4/"+samplename+"CountsForplotting.txt"
                    with open(outTaxpasta,"r") as taxpasta, open(outforplot, "w") as o:
                        print("Taxonomy_nr\tTaxonomy_name\tCounts", file=o)
                        next(taxpasta)
                        next(taxpasta)
                        for l in taxpasta:
                            l=l.strip()
                            TAXID=int(l.split("\t")[0])
                            Counts=int(l.split("\t")[sampleindex])
                            rank=l.split("\t")[2]
                            if rank == "species":
                                if Counts >= dptresh:
                                    speciesid=l.split("\t")[0]
                                    speciesname=l.split("\t")[1]
                                    print(str(speciesid)+"\t"+str(speciesname)+"\t"+str(Counts), file=o)

                                    
def main(taxprofdict,taxdump,dptresh):
    
    ParseMetaphlan4(taxprofdict, taxdump, dptresh)

    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict,args.taxdump, args.dptresh)

    
