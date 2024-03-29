#!/usr/bin/py

import argparse
from Bio import SeqIO
from collections import Counter
from ete3 import NCBITaxa
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
    parser.add_argument("--taxdumpfile", dest = 'taxdump', help ="Path to taxdump (required)", required=True)
    parser.add_argument("--DepthTresh", dest = 'dptresh',default=10,type=int,help ="Minimum depth required to be reported (Default 10)")
    parser.add_argument("--IgnoreReadExtraction", dest = 'IgnoreReadExtraction',help ="If used the reads will not be extracted (Optinal)", action='store_true')
    arguments = parser.parse_args(argv)
    return arguments


def ParseDiamond(taxprofdict, taxdump,dptresh, IgnoreReadExtraction):    
    ncbi = NCBITaxa(taxdump_file=taxdump)
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="diamond"
    Fastqfiles=[]
    for i in subfolders:
        if "bowtie2" in i: # The unmapped fastq should be in the bowtie2 dir
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ]
            for k in subfolders_2:
                if "align" in k:
                    Fastqfiles=glob.glob(k+"/*.fastq.gz") # We have FastqFiles then we can extract
        if tool in i: # We can extract reads using diamond
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in diamond dir
            for k in subfolders_2:
                if "Diamond_" in k:
                    try:
                         os.mkdir("Diamond")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    diamondtsv=glob.glob(k+"/*tsv")
                    for d in diamondtsv: 
                        outcountsforplotting="Diamond/"+d.split("/")[-1].split("_Diamond_230321.diamond.tsv")[0]+"_CountsForplotting.txt"
                        with open(d, "r") as diamondtsv, open(outcountsforplotting, "w") as o: 
                            print("TaxID\tSpecies\tCounts", file=o)
                            SpeciesDic={}
                            for l in diamondtsv: 
                                l=l.strip()
                                readname=l.split("\t")[0]
                                taxid=float(l.split("\t")[1])
                                evalue=float(l.split("\t")[2])
                                if not taxid==0 and not evalue==0: 
                                    lineage = ncbi.get_lineage(taxid)
                                    names = ncbi.get_taxid_translator(lineage)
                                    lineage2ranks = ncbi.get_rank(names)
                                    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
                                    if 'species' in ranks2lineage.keys():
                                        taxidspecies=ranks2lineage['species']
                                        taxidspecies_name=ncbi.get_taxid_translator([taxidspecies])
                                        taxid2taxname_species=taxidspecies_name[taxidspecies]
                                        speciesandtaxid=taxid2taxname_species+"_"+str(taxidspecies)
                                        if not speciesandtaxid in SpeciesDic: 
                                            SpeciesDic[speciesandtaxid]=[readname]
                                        else:
                                            SpeciesDic[speciesandtaxid].append(readname)
                            for k, v in SpeciesDic.items():
                                if len(v) >= dptresh:
                                    print(str(k.split("_")[-1])+"\t"+k.split("_")[0]+"\t"+str(len(v)), file=o)
                            # Extract the reads, outputed from the bowtie directory
                            if Fastqfiles and not IgnoreReadExtraction:
                                outfolderclassifiedreads="Diamond/Classified_Reads/"
                                try:
                                    os.makedirs(outfolderclassifiedreads)
                                except FileExistsError: 
                                    pass
                                if "_se_" in d: 
                                    basename=d.split("/")[-1].split("_Diamond_230321.diamond.tsv")[0].replace("_se_","_") # remove _se_ to be able to link to original fastq
                                if "_pe_" in d: 
                                    basename=d.split("/")[-1].split("_Diamond_230321.diamond.tsv")[0].replace("_se_","_") # remove _pe_ to be able to link to original fastq
                                for f in Fastqfiles:                                    
                                    if basename in f: 
                                        if f.endswith(".gz"): # If the files are gziped you need to use gzip open, save record to dict and get a basename from the fastq
                                            Records=SeqIO.to_dict(SeqIO.parse(gzip.open(f, "rt"),'fastq'))
                                            fname=f.split("/")[-1].replace(".unmapped","").split(".fastq.gz")[0]
                                        else:
                                            Records=SeqIO.to_dict(SeqIO.parse(f,'fastq'))
                                            fname=f.split("/")[-1].replace(".unmapped","").split(".fastq")[0]
                                        for k, v in SpeciesDic.items():
                                            if len(v) >= dptresh: 
                                                taxa=k.split("_")[0].replace(" ","")+"_"+str(k.split("_")[-1])
                                                outfoldersspecies="Diamond/Classified_Reads/"+taxa
                                                try:
                                                    os.makedirs(outfoldersspecies)
                                                except FileExistsError: 
                                                    pass
                                                outfq=outfoldersspecies+"/"+taxa+"_"+fname+".fastq"
                                                with open(outfq, "w") as o: 
                                                    for reads in v:
                                                        rec=Records[reads].format("fastq").strip()
                                                        print(rec, file=o)
                                                



def ParseDiamond_withTaxpasta(taxprofdict, taxdump,dptresh):
    """
    """
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="diamond"
    Fastqfiles=[]
    for i in subfolders:
        if "bowtie2" in i: # The unmapped fastq should be in the bowtie2 dir 
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ]
            for k in subfolders_2:
                if "align" in k:
                    Fastqfiles=glob.glob(k+"/*.fastq.gz") # We have FastqFiles then we can extract

        if tool in i: # We can extract reads using diamond            
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in diamond dir
            for k in subfolders_2:
                if "Diamond_" in k:
                    try:
                         os.mkdir("Diamond")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    # Taxpasta
                    diamonddb=k.split("/")[-1]
                    outTaxpasta= "Diamond/"+ "Diamond_" + diamonddb + ".tsv"
                    command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p diamond -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*.diamond.tsv" %(outTaxpasta, taxdump, k)
                    subprocess.call(command, shell=True)
                    # Printing required files for plotting
                    Sampleorder=[]
                    with open(outTaxpasta, "r") as taxpasta:
                        header = taxpasta.readline().split("\t")
                        for i in header[4:]:
                            samplename=i.split(diamonddb+".diamond")[0]
                            Sampleorder.append([i,samplename,header.index(i)])

                    for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                        Annotation={}
                        samplename=i[1]
                        samplenamewithoutPE="_".join(samplename.rsplit("_pe_", 1)) # Remove the PE that diamond adds to the name

                        sampleindex=i[2]
                        outforplot="Diamond/"+samplenamewithoutPE+"CountsForplotting.txt"
                        
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
                        # Extract the reads, outputed from the bowtie directory
                        if Fastqfiles:
                            #If there is fastq files, loop through the diamond tsv to get the read name, append that to species in a dictionary, for each species generate a folder. 
                            #There will be one folder per species, within that folder there will be all detected samples (2 fastq (PE) ) for that species. 
                            detectedReads={}
                            c=k+"/"+samplename+diamonddb+".diamond.tsv"
                            with open(c, "r") as inf: # get the readname from the classified reads in Diamond
                                for l in inf:
                                    l=l.strip()
                                    readname=l.split("\t")[0]
                                    taxid=int(l.split("\t")[1])
                                    if taxid in Annotation: # It wont be in the annotation if the read is at a level we are not targeting, say we go for species we wont get Genus
                                        Anno=Annotation[taxid]
                                        key=Anno.replace(" ","").replace("(","").replace(")","").replace("/","")+"_"+str(taxid) # Remove space, remove parantesis, remove / from the species names
                                        if key in detectedReads:
                                            detectedReads[key].append(readname)
                                        else:
                                            detectedReads[key]=[readname]
                                            outfoldersspecies="Diamond/Classified_Reads/"+key
                                            try:
                                                os.makedirs(outfoldersspecies) # Create one output folder per species
                                            except FileExistsError: # As we are looping throught the classified reads files there is one per sample, i only want to create one folder per species in the kraken2 folder. If we have one species in more than one sample we need to capture the error.
                                                continue
                            for f in Fastqfiles: # for all classifiers the samples get the _pe_ addition. check if we have this in the sample name, in that case give a warning, if not just replace pe in the name, can we do this in a petter way?
                                if samplenamewithoutPE.strip("_") in f: 
                                    if f.endswith(".gz"): # If the files are gziped you need to use gzip open, save record to dict and get a basename from the fastq
                                        Records=SeqIO.to_dict(SeqIO.parse(gzip.open(f, "rt"),'fastq'))
                                        fname=f.split("/")[-1].replace(".unmapped","").split(".fastq.gz")[0]
                                    else:
                                        Records=SeqIO.to_dict(SeqIO.parse(f,'fastq'))
                                        fname=f.split("/")[-1].replace(".unmapped","").split(".fastq")[0]
                                    for key, values in detectedReads.items(): # Loop through detected reads, for each species for that sample we are extracting from the fastq file
                                        outfq="Diamond/Classified_Reads/"+key+"/"+key+"_"+fname+".fastq" # Out fastq filename
                                        with open(outfq, "w") as o:
                                            for v in values:
                                                try: # We need to handle if the PE info is within the read name, then it is just a try as the exact hit might be in R2 instead of R1. 
                                                    rec=Records[v].format("fastq").strip()
                                        
                                                    print(rec, file=o)
                                                except KeyError:
                                                    continue
                        else:
                            print("Warning, no fastqfiles Available. No read extraction")

def main(taxprofdict, taxdump, dptresh, IgnoreReadExtraction):
    ParseDiamond(taxprofdict, taxdump, dptresh, IgnoreReadExtraction)


if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict, args.taxdump, args.dptresh, args.IgnoreReadExtraction)

