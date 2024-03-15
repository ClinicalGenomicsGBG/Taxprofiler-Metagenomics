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
    parser.add_argument("--Db_sheet", dest = 'dbsheet', help ="Path to Database sheet")
    parser.add_argument("--IgnoreReadExtraction", dest = 'IgnoreReadExtraction',help ="If used the reads will not be extracted (Optional)", action='store_true')
    arguments = parser.parse_args(argv)
    return arguments


def ParseKrakenUniq(taxprofdict, dptresh, dbsheet, IgnoreReadExtraction):
    """
    
    """
    
    # Extract name of the database from the database sheet

    with open(dbsheet, "r") as db:
        next(db)
        for l in db:
            l=l.strip()
            tool="krakenuniq"
            if tool in l.split(",")[0]:
                krakdb=l.split(",")[1]


    print("Parsing KrakenUniq...")
    
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    Fastqfiles=[]
    Annotation={}
    
    for i in subfolders:
        if tool in i: # We can extract reads using krak
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
            for k in subfolders_2:
                if krakdb in k:
                    try:
                        os.mkdir("KrakenUniq")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    reports=glob.glob(f'{k}/*.report.txt')
                    for r in reports: # Looping through the reports! 
                        speciesStrainAnno={} # To keep the species annotation and info if there is a strain annotation, we use this when we extract the detected reads from classified report
                        SpeciesCounts={}
                        samplename=r.split(".krakenuniq.report.txt")[0].split("/")[-1]                    
                        outforplot=f'KrakenUniq/{samplename}_CountsForplotting.txt'
                        speciesindex=0 # to check for the incrementation
                        
                        with open(r, "r") as report, open(outforplot, "w") as o:
                            print("TaxID\tSpecies\tCounts", file=o)
                            for l in report:
                                l=l.strip()
                                if not (l.startswith("#") or l.startswith("%")):
                                    taxnr=l.split("\t")[6]
                                    counts=int(l.split("\t")[1])
                                    taxcounts=int(l.split("\t")[2])
                                    taxlevel=l.split("\t")[7] # Rank
                                    taxname=l.split("\t")[-1].lstrip()
                                    taxaindex=l.split("\t")[-1].split(" ").count("")
                                    if taxlevel=="species":
                                        speciesindex=l.split("\t")[-1].split(" ").count("")
                                        speciesTresh=counts
                                        if counts >= dptresh: # If our taxlevel is species
                                            print(str(taxnr)+"\t"+taxname+"\t"+str(counts), file=o)
                                            speciesStrainAnno[taxname]=[taxnr]
                                            speciesLinkedTostrain=taxname # We save this as we can link the strains to this species
                                            SpeciesCounts[taxname]=counts
                                            
                                    else: # Check the incrementation for read extraction!
                                        if taxaindex < speciesindex: # If we are at a lower incrementation we always need to reset it! 
                                            speciesindex = 0
                                        if not taxcounts == 0 and taxaindex > speciesindex and speciesindex != 0: # We are searching for the increments here, incremented more than species and has a count to it
                                            if speciesTresh >= dptresh:
                                                speciesStrainAnno[speciesLinkedTostrain].append(taxnr)
                                                
                        # To get to the read names we need to go for the classified report first
                        
                        classified=f'{k}/{samplename}.krakenuniq.classified.txt'
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
                        
                        fastaSE=f'{k}/{samplename}.classified.fasta.gz'
                        fastaPE=f'{k}/{samplename}.merged.classified.fasta.gz'
                        
                        if (os.path.exists(fastaSE) or os.path.exists(fastaPE)) and not IgnoreReadExtraction:

                            if os.path.exists(fastaSE):
                                # We have SE reads
                                fasta=fastaSE
                            if os.path.exists(fastaPE):
                                # We have PE reads
                                fasta=fastaPE
                            
                            outfolderClassifiedReads="KrakenUniq/Classified_Reads"
                            try:
                                os.makedirs(outfolderClassifiedReads)
                            except FileExistsError:
                                pass
                            
                            Records=SeqIO.to_dict(SeqIO.parse(gzip.open(fasta, "rt"),'fasta'))
                            for key, values in specieswithReadnames.items():
                                
                                if len(values)>=dptresh: 
                                    speciesIdentifierkey=key.replace(" ","_").replace("(","").replace(")","").replace("/","") # Remove space, remove parantesis, remove from the species names
                                    outfolderspecies=outfolderClassifiedReads+"/"+speciesIdentifierkey
                            
                                    try:
                                        os.makedirs(outfolderspecies)
                                    except FileExistsError:
                                        pass
                                    outfasta='f{outfolderspecies}/{speciesIdentifierkey}_{samplename}.fasta'
                                    with open(outfasta,  "w") as o:
                                        for reads in values:
                                            rec=Records[reads].format("fasta").strip()
                                            print(rec, file=o)



                                            
def main(taxprofdict, dptresh, dbsheet, IgnoreReadExtraction):
    ParseKrakenUniq(taxprofdict, dptresh, dbsheet, IgnoreReadExtraction)
    
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict, args.dptresh, args.dbsheet, args.IgnoreReadExtraction)

