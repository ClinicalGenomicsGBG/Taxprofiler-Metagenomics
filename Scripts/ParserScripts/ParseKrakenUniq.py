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
    parser.add_argument("--DepthTresh", dest = 'dptresh',default=10, type=int, help ="Minimum depth required to be reported")
    arguments = parser.parse_args(argv)
    return arguments



def ParseKrakenUniq(taxprofdict, dptresh):
    """
    
    """

    print("Parsing KrakenUniq...")
    
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="krakenuniq"
    Fastqfiles=[]
    Annotation={}
    
    for i in subfolders:
        if tool in i: # We can extract reads using krak
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
            for k in subfolders_2:
            
                if "krakenuniq_MicrobialDB" in k:
                    krakdb=k.split("/")[-1]
                    try:
                        os.mkdir("KrakenUniq")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    reports=glob.glob(k+"/*.report.txt")
                    for r in reports: # Looping through the reports! 
                        speciesStrainAnno={} # To keep the species annotation and info if there is a strain annotation, we use this when we extract the detected reads from classified report
                        SpeciesCounts={}
                        samplename=r.split(".krakenuniq.report.txt")[0].split("/")[-1]                    
                        outforplot="KrakenUniq/"+samplename+"_CountsForplotting.txt"

                        with open(r, "r") as report, open(outforplot, "w") as o:
                            print("TaxID\tSpecies\tCounts", file=o)
                            for l in report:
                                l=l.strip()
                                if not (l.startswith("#") or l.startswith("%")):
                                    taxnr=l.split("\t")[6]
                                    counts=int(l.split("\t")[1])
                                    taxlevel=l.split("\t")[7] # Rank
                                    taxname=l.split("\t")[-1].lstrip()
                                    if taxlevel=="species":
                                        speciesTresh=counts
                                        if counts >= dptresh: # If our taxlevel is species
                                            print(str(taxnr)+"\t"+taxname+"\t"+str(counts), file=o)
                                            speciesStrainAnno[taxname]=[taxnr]
                                            speciesLinkedTostrain=taxname # We save this as we can link the strains to this species
                                            SpeciesCounts[taxname]=counts
                                        
                                    elif taxlevel=="subspecies" or taxlevel=="strain" : # There might also be sub species and strains, 
                                        if speciesTresh >= dptresh: # We need to make sure that the species treshold is more than or equal to the depthtreshold, if the species is ok we append the strains if they are there!
                                            speciesStrainAnno[speciesLinkedTostrain].append(taxnr)



                        # To get to the read names we need to go for the classified report first

                        classified=k+"/"+samplename+".krakenuniq.classified.txt"

                        if os.path.exists(classified):
                            specieswithReadnames={}
                            with open(classified, "r") as inf:
                                for l in inf:
                                    l=l.strip()
                                    detected=l.split("\t")[0]
                                    readname=l.split("\t")[1]
                                    taxid=l.split("\t")[2]
                                    for key, values in speciesStrainAnno.items(): 
                                        if taxid in values:
                                            if not key in specieswithReadnames:
                                                specieswithReadnames[key]=[readname]
                                            else:
                                                specieswithReadnames[key].append(readname)
                            
                        # Only for SE as it is now, PE is missing untill update from taxprofiler! 
                        fastq=k+"/"+samplename+".classified.fastq.gz"

    
                        if os.path.exists(fastq):

                            outfolderClassifiedReads="KrakenUniq/Classified_Reads"
                            try:
                                os.makedirs(outfolderClassifiedReads)
                            except FileExistsError:
                                pass
                            
                            
                            Records=SeqIO.to_dict(SeqIO.parse(gzip.open(fastq, "rt"),'fastq'))



                            for key, values in specieswithReadnames.items():


                                if len(values)>=dptresh: 

                                    speciesIdentifierkey=key.replace(" ","_").replace("(","").replace(")","").replace("/","") # Remove space, remove parantesis, remove from the species names


                                    outfolderspecies=outfolderClassifiedReads+"/"+speciesIdentifierkey
                                    

                                    try:
                                        os.makedirs(outfolderspecies)
                                    except FileExistsError:
                                        pass
                                    outfastq=outfolderspecies+"/"+speciesIdentifierkey+"_"+samplename+".fastq"


                                    with open(outfastq,  "w") as o:
                                        for reads in values:
                                            rec=Records[reads].format("fastq").strip()
                                            print(rec, file=o)

                                    
def main(taxprofdict, dptresh):
    ParseKrakenUniq(taxprofdict, dptresh)
    
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict, args.dptresh)

