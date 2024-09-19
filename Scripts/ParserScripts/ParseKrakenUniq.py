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

                    try:
                        os.mkdir("KrakenUniq/Extras")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                        
                    reports=glob.glob(f'{k}/*.report.txt')
                    for r in reports: # Looping through the reports! 
                        speciesStrainAnno={} # To keep the species annotation and info if there is a strain annotation, we use this when we extract the detected reads from classified report
                        SpeciesCounts={}
                        samplename=r.split(".krakenuniq.report.txt")[0].split("/")[-1]                    
                        outforplot=f'KrakenUniq/{samplename}_CountsForplotting.txt'
                        outforplot_discaredSpecies=f'KrakenUniq/Extras/{samplename}_CountsForplotting_DiscaredSpecies.txt' # Species we remove when we are at lower treshold
                        outforplot_extras=f'KrakenUniq/Extras/{samplename}_CountsForplotting_Others.txt' # Groups removed above species
                        SpeciesDomainLinkage=f'KrakenUniq/Extras/{samplename}_SpeciesDomainLinkage.txt'

                        speciesindex=0 # to check for the incrementation

                        Anyspecieshits="no"
                        
                        with open(r, "r") as report, open(outforplot, "w") as o, open(outforplot_discaredSpecies, "w") as o_discardedSpecies, open(outforplot_extras, "w") as o_discarded_extras, open(SpeciesDomainLinkage, "w") as o_SpeciesDomainLinkage:
                            print("TaxID\tSpecies\tCounts", file=o)
                            print("TaxID\tSpecies\tCounts", file=o_discardedSpecies)
                            print("TaxID\tTaxname\tLineage\tCounts", file=o_discarded_extras) 
                            print("TaxID\tTaxname\tDomain", file=o_SpeciesDomainLinkage)
                            for l in report:
                                l=l.strip()
                                if not (l.startswith("#") or l.startswith("%")):
                                    taxnr=l.split("\t")[6]
                                    counts=int(l.split("\t")[1])
                                    taxcounts=int(l.split("\t")[2])
                                    taxlevel=l.split("\t")[7] # Rank
                                    taxname=l.split("\t")[-1].lstrip()
                                    taxaindex=l.split("\t")[-1].split(" ").count("")
                                    if taxlevel == "superkingdom" or taxname=="other entries": # In krakenUniq we have other entires as well, if not using this they end up behind archeae
                                        CurrentLineage=taxname
                                    if Anyspecieshits == "no":
                                        print(f'{taxnr}\t{taxname}\t{taxlevel}\t{counts}', file=o_discarded_extras) # We have not hit species yet, that it why we cannot use the first increment
                                    if taxlevel=="species":
                                        Anyspecieshits = "yes" # we hit our first species so we can use the increment now
                                        speciesindex=l.split("\t")[-1].split(" ").count("")
                                        speciesTresh=counts
                                        if counts >= dptresh: # If our taxlevel is species
                                            print(str(taxnr)+"\t"+taxname+"\t"+str(counts), file=o)
                                            print(f'{taxnr}\t{taxname}\t{CurrentLineage}', file=o_SpeciesDomainLinkage)
                                            speciesStrainAnno[taxname]=[taxnr]
                                            speciesLinkedTostrain=taxname # We save this as we can link the strains to this species
                                            SpeciesCounts[taxname]=counts
                                        else: # We are at species but we are not meeting the user defined treshold
                                            print(f'{taxnr}\t{taxname}\t{counts}', file=o_discardedSpecies)
                                    else: # Check the incrementation for read extraction! This is due to some bein subspecies, some are strain and some are no rank. They are all included in the total of 
                                        if taxaindex < speciesindex: # If we are at a lower incrementation we always need to reset it! 
                                            speciesindex = 0

                                            # We are futher up than the species index, then we can save it to discarded taxas above species
                                            print(f'{taxnr}\t{taxname}\t{taxlevel}\t{counts}', file=o_discarded_extras) # We would only hit this one if we hit species once before! 
                                            
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
                            
                        # Issue, PE outputs fasta, SE outputs fastq that is name fa.
                        
                        fastqSE=f'{k}/{samplename}.classified.fastq.gz'
                        fastaPE=f'{k}/{samplename}.merged.classified.fasta.gz'
                        
                        if (os.path.exists(fastqSE) or os.path.exists(fastaPE)) and not IgnoreReadExtraction:

                            if os.path.exists(fastqSE):
                                # We have SE reads
                                Records=SeqIO.to_dict(SeqIO.parse(gzip.open(fastqSE, "rt"),'fastq'))
                            if os.path.exists(fastaPE):
                                # We have PE reads
                                Records=SeqIO.to_dict(SeqIO.parse(gzip.open(fastaPE, "rt"),'fasta'))
                                
                            outfolderClassifiedReads="KrakenUniq/Classified_Reads"
                            try:
                                os.makedirs(outfolderClassifiedReads)
                            except FileExistsError:
                                pass
                            

                            for key, values in specieswithReadnames.items():
                                
                                if len(values)>=dptresh: 
                                    speciesIdentifierkey=key.replace(" ","_").replace("(","").replace(")","").replace("/","") # Remove space, remove parantesis, remove from the species names
                                    outfolderspecies=outfolderClassifiedReads+"/"+speciesIdentifierkey
                            
                                    try:
                                        os.makedirs(outfolderspecies)
                                    except FileExistsError:
                                        pass
                                    outfasta=f'{outfolderspecies}/{speciesIdentifierkey}_{samplename}.fasta'
                                    with open(outfasta,  "w") as o:
                                        for reads in values:
                                            rec=Records[reads].format("fasta").strip()
                                            print(rec, file=o)



                                            
def main(taxprofdict, dptresh, dbsheet, IgnoreReadExtraction):
    ParseKrakenUniq(taxprofdict, dptresh, dbsheet, IgnoreReadExtraction)
    
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict, args.dptresh, args.dbsheet, args.IgnoreReadExtraction)

