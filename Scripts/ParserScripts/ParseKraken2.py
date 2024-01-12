#!/usr/bin/py

import argparse

from Bio.SeqIO.QualityIO import FastqGeneralIterator
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
    parser.add_argument("--DepthTresh", dest = 'dptresh',default=10, type=int, help ="Minimum depth required to be reported (Default 10)")
    parser.add_argument("--IgnoreReadExtraction", dest = 'IgnoreReadExtraction',help ="If used the reads will not be extracted (Optional)", action='store_true')
    arguments = parser.parse_args(argv)
    return arguments



def ParseKraken2(taxprofdict, dptresh, IgnoreReadExtraction):
    """
    
    """

    print("Parsing Kraken2...")


    
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="kraken2"
    Fastqfiles=[]
    Annotation={}
    
    for i in subfolders:
        if tool in i: # We can extract reads using krak
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
            for k in subfolders_2:
                if "Kraken2_" in k:
                    krakdb=k.split("/")[-1]
                    try:
                        os.mkdir("Kraken2")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    reports=glob.glob(k+"/*.report.txt")
                    for r in reports: # Looping through the reports! 
                        speciesStrainAnno={} # To keep the species annotation and info if there is a strain annotation, we use this when we extract the detected reads from classified report
                        SpeciesCounts={}
                        samplename=r.split(".kraken2.kraken2.report.txt")[0].split("/")[-1]                    
                        outforplot="Kraken2/"+samplename+"_CountsForplotting.txt"
                        
                        with open(r, "r") as report, open(outforplot, "w") as o:
                            print("TaxID\tSpecies\tCounts", file=o)
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
                                        SpeciesCounts[taxname]=counts
                                        
                                elif re.findall(r'(S(\w)+)',taxlevel): # If there is something after S, these are after the  strain string
                                    if speciesTresh >= dptresh: # We need to make sure that the species treshold is more than or equal to the depthtreshold, if the species is ok we append the strains if they are there!
                                        speciesStrainAnno[speciesLinkedTostrain].append(taxnr)

                        if "_pe_" in r: # We have the PE flag added by taxprofiler, this was run in Paired END mode!
                            R1=k+"/"+samplename+".kraken2.classified_1.fastq.gz"
                            R2=k+"/"+samplename+".kraken2.classified_2.fastq.gz"
                            
                            fastqs=[R1, R2]
                            for f in fastqs:
                                if os.path.exists(f) and not IgnoreReadExtraction:
                                    SpeciesWithFastq={}
                                    with gzip.open(f, "rt") as handle:
                                        for header, seq, qual in FastqGeneralIterator(handle):
                                            fastqread="@"+header+"\n"+seq+"\n+\n"+qual
                                            taxid=header.split("|")[-1]
                                            for key, values in speciesStrainAnno.items():
                                                if taxid in values:
                                                    if not key in SpeciesWithFastq:
                                                        SpeciesWithFastq[key]=[(fastqread)]
                                                    else:
                                                        SpeciesWithFastq[key].append((fastqread))

                                    for key, values in SpeciesWithFastq.items():
                                        speciesIdentifierkey=key.replace(" ","_").replace("(","").replace(")","").replace("/","") # Remove space, remove parantesis, remove from the species names
                                        outfoldersspecies="Kraken2/Classified_Reads/"+speciesIdentifierkey
                                        try: # Generate the species folder if it is not there!
                                            os.makedirs(outfoldersspecies)
                                        except FileExistsError:
                                            pass
                                        sname=samplename.split(k)[0] # Remove the database part from the sample name!

                                        # R1 or R2?
                                        if "_1.fastq.gz" in f: 
                                            OutFastq="Kraken2/Classified_Reads/"+speciesIdentifierkey+"/"+sname +"_" +speciesIdentifierkey+"_1.fastq"
                                        elif  "_2.fastq.gz" in f:
                                            OutFastq="Kraken2/Classified_Reads/"+speciesIdentifierkey+"/"+sname +"_" +speciesIdentifierkey+"_2.fastq"

                                        else:
                                            print("Error: Classified fastq does not contain R1 or R2!")                                            
                                        with open(OutFastq, "w") as o:
                                            for read in values: # One read as item in the dictionary
                                                print(read, file=o)
                                                # Check that the amount of reported counts is the same as the amount of reported fastq reads!
                                        if not SpeciesCounts[key] == len(values):
                                            print("For sample %s, species %s, reads in fastq: %s, reported counts: %s" %(sname,key, len(values), SpeciesCounts[key]))

                                        
                        elif "_se_" in r: # We have the SE flag added by taxprofiler, this was run in Single END mode!
                            fastq=k+"/"+samplename+".kraken2.classified.fastq.gz" # We need to have it in gzip format!
                            if os.path.exists(fastq) and not IgnoreReadExtraction: # We have the classified reads fastq in kraken2 out, therefore we can extract the reads                            
                                SpeciesWithFastq={}
                                with gzip.open(fastq, "rt") as handle:
                                    for header, seq, qual in FastqGeneralIterator(handle):
                                        fastqread="@"+header+"\n"+seq+"\n+\n"+qual
                                        taxid=header.split("|")[-1]
                                        for key, values in speciesStrainAnno.items():
                                            if taxid in values:
                                                if not key in SpeciesWithFastq:
                                                    SpeciesWithFastq[key]=[(fastqread)]
                                                else:
                                                    SpeciesWithFastq[key].append((fastqread))
                                                    
                                for key, values in SpeciesWithFastq.items():
                                    speciesIdentifierkey=key.replace(" ","_").replace("(","").replace(")","").replace("/","") # Remove space, remove parantesis, remove from the species names
                                    outfoldersspecies="Kraken2/Classified_Reads/"+speciesIdentifierkey
                                    try: # Generate the species folder if it is not there!
                                        os.makedirs(outfoldersspecies)
                                    except FileExistsError:
                                        x="Folder already exists, continue"
                                    sname=samplename.split(k)[0] # Remove the database part from the sample name!
                                    OutFastq="Kraken2/Classified_Reads/"+speciesIdentifierkey+"/"+sname +"_" +speciesIdentifierkey+".fastq"
                                    with open(OutFastq, "w") as o:
                                        for read in values: # One read as item in the dictionary
                                            print(read, file=o)
                                    # Check that the amount of reported counts is the same as the amount of reported fastq reads!
                                    if not SpeciesCounts[key] == len(values):            
                                        print("For sample %s, species %s, reads in fastq: %s, reported counts: %s" %(sname,key, len(values), SpeciesCounts[key]))



                                    
def main(taxprofdict, dptresh, IgnoreReadExtraction):
    ParseKraken2(taxprofdict, dptresh, IgnoreReadExtraction)
    
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict, args.dptresh, args.IgnoreReadExtraction)

