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

def ParseKrakenUniq(taxprofdict, taxdump, dptresh):
    """
    Parsing KrakenUniq
    """
    
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="krakenuniq"
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
                if "krakenuniq_" in k:
                    try:
                        os.mkdir("KrakenUniq")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    # Taxpasta
                    krakdb=k.split("/")[-1]
                    outTaxpasta= "KrakenUniq/"+ "krakenuniq_" + krakdb+".tsv"
                    command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p krakenuniq -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*.krakenuniq.report.txt" %(outTaxpasta, taxdump, k)

                    subprocess.call(command, shell=True)
                    Sampleorder=[]
                    with open(outTaxpasta, "r") as taxpasta:
                        header = taxpasta.readline().split("\t")
                        for i in header[4:]:
                            samplename=i.split(".unmapped.krakenuniq.report")[0]
                            Sampleorder.append([i,samplename,header.index(i)])
                                
                    for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                        samplename=i[1]
                        sampleindex=i[2]
                        outforplot="KrakenUniq/"+samplename+"_CountsForplotting.txt"

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
                                        if not speciesid in Annotation:
                                            Annotation[int(speciesid)]=speciesname


                        # Extract the reads, outputed from the bowtie directoty
                    if Fastqfiles:
                        classifiedreads=glob.glob(k+"/*.classified.txt")
                        for c in classifiedreads:
                            detectedReads={}
                            samplename=c.split(".unmapped.krakenuniq.classified.txt")[0].split("/")[-1].strip("_")
                            samplename="_".join(samplename.rsplit("_pe_", 1)) # Remove the PE that krakenUniq adds to the name
                            with open(c, "r") as inf: # get the readname from the classified reads in krakenuniq
                                for l in inf:
                                    if l.split("\t")[0] == "C": # If classified extract read id
                                        readname=l.split("\t")[1]
                                        taxid=int(l.split("\t")[2])
                                        if taxid in Annotation: # It wont be in the annotation if the read is at a level we are not targeting, say we go for species we wont get Genus  
                                            Anno=Annotation[taxid]
                                            key=Anno.replace(" ","").replace("(","").replace(")","").replace("/","")+"_"+str(taxid) # Remove space, remove parantesis, remove / from the species names
                                            if key in detectedReads:
                                                detectedReads[key].append(readname)
                                            else:
                                                detectedReads[key]=[readname]
                                                outfoldersspecies="KrakenUniq/Classified_Reads/"+key
                                                try:
                                                    os.makedirs(outfoldersspecies) # Create one output folder per species 
                                                except FileExistsError: # As we are looping throught the classified reads files there is one per sample, i only want to create one folder per species in the kraken2 folder. If we have one species in more than one sample we need to capture the error. 
                                                    continue
                                                    
                            for f in Fastqfiles: # for all classifiers the samples get the _pe_ addition. check if we have this in the sample name, in that case give a warning, if not just replace pe in the name.
                                if samplename in f:  
                                    if f.endswith(".gz"): # If the files are gziped you need to use gzip open, save record to dict and get a basename from the fastq
                                        Records=SeqIO.to_dict(SeqIO.parse(gzip.open(f, "rt"),'fastq'))
                                        fname=f.split("/")[-1].replace(".unmapped","").split(".fastq.gz")[0]
                                    else:
                                        Records=SeqIO.to_dict(SeqIO.parse(f,'fastq'))
                                        fname=f.split("/")[-1].replace(".unmapped","").split(".fastq")[0]
                                    for key, values in detectedReads.items(): # Loop through detected reads, for each species for that sample we are extracting from the fastq file
                                        outfq="KrakenUniq/Classified_Reads/"+key+"/"+key+"_"+fname+".fastq" # Out fastq filename
                                        with open(outfq, "w") as o: 
                                            for i in values:
                                                rec=Records[i].format("fastq").strip()                                                        
                                                print(rec, file=o)                                            
                               # else: 
                                #    print("We cannot have pe in the readname, need to think about this!  ")
                                 #   continue                                                
                    else:
                        print("Warning, no fastqfiles Available. No read extraction")                
    
                                    
def main(taxprofdict, taxdump, dptresh):
    ParseKrakenUniq(taxprofdict, taxdump,dptresh)

if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict,args.taxdump, args.dptresh)

