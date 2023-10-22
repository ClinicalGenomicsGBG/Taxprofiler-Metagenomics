#!/usr/bin/py

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
import re
import subprocess
import sys
import time
import xlsxwriter


def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes the output from taxprofiler and parses it')
    parser.add_argument("--Taxprofiler_out", dest = 'taxprofdict', required=True, help ="Output folder from taxprofiler (required)")
    parser.add_argument("--taxdumpfile", dest = 'taxdump', required=True, help ="Path to taxdump (required)")
    parser.add_argument("--DepthTresh", dest = 'dptresh',default=10, type=int, help ="Minimum depth required to be reported")
    arguments = parser.parse_args(argv)
    return arguments




def ParseKraken2(taxprofdict, taxdump, dptresh):
    """
    
    """

    print("Parsing Kraken2...")
    
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="kraken2"
    Fastqfiles=[]
    Annotation={}
    
    for i in subfolders:
        if "bowtie2" in i:
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ]
            for k in subfolders_2:
                if "align" in k:
                    Fastqfiles=glob.glob(k+"/*.fastq.gz") # We have FastqFiles then we can extract
        if tool in i: # We can extract reads using krak
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
            for k in subfolders_2:
                if "krak_" in k:
                    krakdb=k.split("/")[-1]
                    try:
                        os.mkdir("Kraken2")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())



                    reports=glob.glob(k+"/*.report.txt")

                    for r in reports: # Looping through the reports! 
                        speciesStrainAnno={} # To keep the species annotation and info if there is a strain annotation, we use this when we extract the detected reads from classified report
                        print(r)
                        samplename=r.split(krakdb+".kraken2.kraken2.report")[0].split("/")[-1]
                        if "_pe_" in samplename:
                            samplename_base="_".join(samplename.rsplit("_pe_", 1)).rstrip("_") # Remove the PE that kraken2 adds to the name, if PE reads
                        if "_se_" in samplename:
                            samplename_base="_".join(samplename.rsplit("_se_", 1)).rstrip("_") # Remove the SE that kraken2 adds to the name, if SE reads
                        outforplot="Kraken2/"+samplename_base+"_CountsForplotting.txt"
                        print(outforplot)
                        with open(r, "r") as report, open(outforplot, "w") as o:
                            for l in report:
                                l=l.strip()
                                taxnr=l.split("\t")[4]
                                counts=int(l.split("\t")[1])
                                taxlevel=l.split("\t")[3]
                                taxname=l.split("\t")[5].lstrip()
                                if taxlevel=="S":
                                    speciesTresh=counts
                                    if counts >= dptresh: # If our taxlevel is species
                                        print(str(taxnr)+"\t"+taxname+"\t"+str(counts), file=o)
                                        speciesStrainAnno[taxname]=[taxnr]
                                        speciesLinkedTostrain=taxname # We save this as we can link the strains to this species
                                elif re.findall(r'(S(\w)+)',taxlevel): # If there is something after S, these are after the species string
                                    if speciesTresh >= dptresh: # We need to make sure that the species treshold is more than or equal to the depthtreshold, if the species is ok we append the strains if they are there!
                                        speciesStrainAnno[speciesLinkedTostrain].append(taxnr)
                                    
                                    
                         # Extract the reads, outputed from the bowtie directoty
                        if Fastqfiles:
                            detectedReads={}
                            print("Create the species subfolders for all detected reads and extracts the reads from kraken2 classified reads report")
                            classifiedreads=glob.glob(k+"/*.classifiedreads.txt")

                            for c in classifiedreads:
                                if samplename in c:
                                    with open(c, "r") as classifiedreadreport:
                                        for l in classifiedreadreport:
                                            l=l.strip()
                                            if l.split("\t")[0] == "C": # If classified extract read id
                                                readname=l.split("\t")[1]
                                                taxid=int(l.split("\t")[2])
                                                for key, values in speciesStrainAnno.items():
                                                    if str(taxid) in values: # The species identifier is always value[0], strain identifier will be the additional items in the value list

                                                        speciesIdentifierkey=key.replace(" ","").replace("(","").replace(")","").replace("/","")+"_"+str(values[0]) # Remove space, remove parantesis, remove / from the species names
                                                        outfoldersspecies="Kraken2/Classified_Reads/"+speciesIdentifierkey

                                                        if not speciesIdentifierkey in detectedReads:
                                                            detectedReads[speciesIdentifierkey]=[readname]
                                                        else:
                                                            detectedReads[speciesIdentifierkey].append(readname)
                                                            
                                                        try:
                                                            os.makedirs(outfoldersspecies)
                                                        except FileExistsError:
                                                            continue                                                        
                            for f in Fastqfiles:
                                if samplename_base in f:
                                    if f.endswith(".gz"): # If the files are gziped you need to use gzip open, save record to dict and get a basename from the fastq
                                        Records=SeqIO.to_dict(SeqIO.parse(gzip.open(f, "rt"),'fastq'))
                                        fname=f.split("/")[-1].replace(".unmapped","").split(".fastq.gz")[0]
                                    else:
                                        Records=SeqIO.to_dict(SeqIO.parse(f,'fastq'))
                                        fname=f.split("/")[-1].replace(".unmapped","").split(".fastq")[0]

                                    for key, values in detectedReads.items(): # Loop through detected reads, for each species for that sample we are extracting from the fastq file
                                        outfq="Kraken2/Classified_Reads/"+key+"/"+key+"_"+fname+".fastq" # Out fastq filename
                                        print(outfq)
                                        with open(outfq, "w") as o:
                                            Counter=0
                                            for v in values:
                                                try: # To allow for PE info within the read header
                                                    rec=Records[v].format("fastq").strip()
                                                    print(rec, file=o)
                                                    Counter+=1
                                                except KeyError:
                                                    continue

                                            if not len(values)==Counter: # Check so the amount of extracted reads is the same in fastq as the countfile
                                                print("Warning only %s reads extracted, should be %s" %(Counter,len(values)) )
                        else:
                            print("Warning no fastqfiles Available, no Read extraction")

                                                        
                            
def main(taxprofdict, taxdump, dptresh):
    ParseKraken2(taxprofdict, taxdump, dptresh)
    
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict,args.taxdump, args.dptresh)

