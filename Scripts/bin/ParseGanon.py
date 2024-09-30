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


def ParseGanon(taxprofdict, taxdump,dptresh):
    """
    """

    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="ganon"
    for i in subfolders:
        subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in ganon dir
        for k in subfolders_2:
            if "Manual" in k:
                try:
                    os.mkdir("Ganon")
                except FileExistsError:
                    pass
                # Taxpasta
                db=k.split("/")[-1]
                outTaxpasta= "Ganon/Ganon_" + db + ".tsv"
                command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p ganon -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*.tre" %(outTaxpasta, taxdump, k)
                #subprocess.call(command, shell=True)


                for tre in glob.glob(k+"/"+"*.tre"): # open the tre files generated from ganon report
                    samplename=tre.split("/")[-1].split("_report.tre")[0]
                    
                    outforplot="Ganon/"+samplename+"_CountsForplotting_2.txt"
                    with open(tre,"r") as inf, open(outforplot, "w") as o:
                        print("Taxonomy_nr\tTaxonomy_name\tCounts", file=o)
                        for l in inf:
                            l=l.strip()
                            rank=l.split("\t")[0]
                            TaxonomyID=l.split("\t")[1]
                            TaxonomyName=l.split("\t")[3]
                            UniqueReads=int(l.split("\t")[4])
                            if rank == "species":
                                if UniqueReads >= dptresh:
                                    print(str(TaxonomyID)+ "\t"+str(TaxonomyName)+"\t"+str(UniqueReads), file=o)
                
                """ 
                # Generate the count files from taxpasta, skip this if we have issues with the taxdump used 
                # Printing required files for plotting
                Sampleorder=[]
                with open(outTaxpasta, "r") as taxpasta:
                    header = taxpasta.readline().split("\t")
                    for i in header[4:]:                        
                        samplename=i.split("_report")[0]
                        Sampleorder.append([i,samplename,header.index(i)])
                for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                    samplename=i[1]
                    samplename="_".join(samplename.rsplit("_pe_", 1)) # Remove the PE that Ganon adds to the name
                    sampleindex=i[2]
                    outforplot="Ganon/"+samplename+"_CountsForplotting.txt"
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

                """
                                    
def main(taxprofdict,taxdump,dptresh):
    ParseGanon(taxprofdict, taxdump, dptresh)

    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict,args.taxdump, args.dptresh)

    
