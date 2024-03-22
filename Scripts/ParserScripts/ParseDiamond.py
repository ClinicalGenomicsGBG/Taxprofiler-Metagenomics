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
import random
import string


def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes the output from taxprofiler and parses it')
    parser.add_argument("--Taxprofiler_out", dest = 'taxprofdict', required=True, help ="Output folder from taxprofiler (required)")
    parser.add_argument("--taxdumpfile", dest = 'taxdump', help ="Path to taxdump (required)", required=True)
    parser.add_argument("--DepthTresh", dest = 'dptresh',default=10,type=int,help ="Minimum depth required to be reported (Default 10)")
    parser.add_argument("--Db_sheet", dest = 'dbsheet', help ="Path to Database sheet")
    parser.add_argument("--IgnoreReadExtraction", dest = 'IgnoreReadExtraction',help ="If used the reads will not be extracted (Optinal)", action='store_true')
    arguments = parser.parse_args(argv)
    return arguments

def ParseDiamond(taxprofdict, taxdump, dptresh,dbsheet, IgnoreReadExtraction):

    # Extract name of the database from the database sheet
    with open(dbsheet, "r") as db:
        next(db)
        for l in db:
            l=l.strip()
            tool="diamond"
            if tool in l.split(",")[0]:
                diamonddb=l.split(",")[1]
                
    print("Parsing Diamond...")
    
    letters = string.ascii_lowercase
    tmp=''.join(random.choice(letters) for i in range(5))

    os.mkdir(tmp)
    
    
    # We need to 
    
    ncbi = NCBITaxa(taxdump_file=taxdump, dbfile=f'{tmp}/taxa.sqlite', update=False)
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
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
                if diamonddb in k:
                    try:
                         os.mkdir("Diamond")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    diamondtsv=glob.glob(k+"/*tsv")
                    try:
                        os.mkdir("Diamond/Extras")
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                        
                    for d in diamondtsv:
                        
                        #samplename=d.split("/")[-1].split(".diamond.tsv")[0]
                        samplename=d.split("/")[-1].split(diamonddb)[0].strip("_")
                        
                        outcountsforplotting=f'Diamond/{samplename}_CountsForplotting.txt'
                        outcountsforplotting_discarded=f'Diamond/Extras/{samplename}_CountsForplotting_DiscardedSpecies.txt'
                        outcounts_extras=f'Diamond/Extras/{samplename}_CountsForplotting_Others.txt'
                        SpeciesDomainLinkage=f'Diamond/Extras/{samplename}_SpeciesDomainLinkage.txt'
                        outcounts_log=f'Diamond/Extras/{samplename}.log'
                        #TaxidsForDiamond={}

                        TaxidsForDiamond=[]
                        
                        # Read the diamond tsv to create a dictionary with taxid and readname
                        with open(d, "r") as diamondtsv:
                            for l in diamondtsv:
                                l=l.strip()
                                readname=l.split("\t")[0]
                                taxid=int(l.split("\t")[1])
                                evalue=float(l.split("\t")[2])
                                if not taxid==0 and not evalue==0.0:
                                    if not taxid in TaxidsForDiamond:
                                        TaxidsForDiamond.append(taxid)
                                        
                        # Next generate a dictionary with species information from this dict! 
                        # My taxid as key, the species name and the taxid as item, 
                        SpeciesDict={} # Incldues all species but also everything below
                        EverythingAboveSpecies={}
                        for taxid in TaxidsForDiamond:                            
                            lineage = ncbi.get_lineage(taxid)
                            names = ncbi.get_taxid_translator(lineage)
                            lineage2ranks = ncbi.get_rank(names)

                            # We want to swap so the key is the rank name and the taxid is the value, the issue is that root  from names get no rank together with the possible additional no rank. Therefore we should drop taxid 1 from lineage2ranks (the root)
                            # We want to swap so the key is the rank name and the taxid is the value, the issue is that root  from names get no rank together with the possible additional no rank. Therefore we should drop taxid 1 from lineage2ranks (the root)
                            #ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items()) # Here we swap to rankname is key and taxid is value
                            # no rank might pop up many times, need to fix that!

                            ranks2lineage ={}
                            for k, v in lineage2ranks.items():
                                ranks2lineage[v] = ranks2lineage.get(v, []) + [k]
                                                     
                            if 'species' in ranks2lineage.keys():
                                taxidspecies=ranks2lineage['species'][0] # There should be only one species, the nested lists are for the no ranks which contains more hits
                                superkingdom=ranks2lineage['superkingdom'][0] # Taxid for superkingdom
                                superkingdom_name=names[superkingdom]
                                taxidspecies_name=ncbi.get_taxid_translator([taxidspecies])
                                taxid2taxname_species=taxidspecies_name[taxidspecies]                                
                                speciesandtaxid=f'{taxid2taxname_species}_{taxidspecies}_{superkingdom_name}'
                                if not speciesandtaxid in SpeciesDict:
                                    SpeciesDict[taxid]=speciesandtaxid

                            else:
                                if not taxid in EverythingAboveSpecies:
                                    for key, values in ranks2lineage.items(): # To get the rank for our alread
                                        for v in values: 
                                            if v == taxid: # Get rank like phylum
                                                rank=key
                                    taxname=ncbi.get_taxid_translator([taxid])
                                    taxname=taxname[taxid]
                                    EverythingAboveSpecies[taxid]=(rank,taxname)
                                    
                                    
                        # Species with reads, loop throught he tsv again, 
                        with open(d, "r") as diamondtsv, open(outcountsforplotting, "w") as o, open(outcounts_extras, "w") as o2, open(outcountsforplotting_discarded, "w") as o3, open(outcounts_log, "w") as o4, open(SpeciesDomainLinkage, "w") as o5:
                            ntotreads=0
                            ClassifiedReads=0
                            nReadsinSpecies=0
                            nReadsSpecies_discarded=0
                            nReadsHighTaxLevels=0
                            print("TaxID\tSpecies\tCounts", file=o)
                            print("TaxID\tSpecies\tCounts", file=o3)
                            print("TaxID\tTaxname\tLineage\tCounts", file=o2)
                            print("TaxID\tTaxname\tDomain", file=o5)
                            Species_withReads={}
                            EverythingAboveSpecies_withReads={}
                            for l in diamondtsv: 
                                l=l.strip()
                                ntotreads+=1
                                readname=l.split("\t")[0]
                                taxid=int(l.split("\t")[1])
                                evalue=float(l.split("\t")[2])
                                if not taxid==0 and not evalue==0.0:
                                    ClassifiedReads+=1
                                    if taxid in SpeciesDict:
                                        specieswithtaxid=SpeciesDict[taxid]
                                        #speciesandtaxid=species+"_"+str(taxid)                                        
                                        if not specieswithtaxid in Species_withReads: 
                                            Species_withReads[specieswithtaxid]=[readname]
                                        else:
                                            Species_withReads[specieswithtaxid].append(readname)

                                    else: # We are above species
                                        if not taxid in EverythingAboveSpecies_withReads:
                                            EverythingAboveSpecies_withReads[taxid]=[readname]
                                        else:
                                            EverythingAboveSpecies_withReads[taxid].append(readname)

                            #print(Species_withReads)
                            for k, v in Species_withReads.items():
                                taxid=k.split("_")[-2]
                                domain=k.split("_")[-1]
                                taxname="_".join(k.split("_")[:-2])
                                counts=len(v)
                                if counts >= dptresh:
                                    #print(str(k.split("_")[-1])+"\t"+"_".join(k.split("_")[:-1])+"\t"+str(len(v)), file=o) # Need to take everything except as some species have more _ than one!
                                    print(f'{taxid}\t{taxname}\t{counts}', file=o)
                                    print(f'{taxid}\t{taxname}\t{domain}', file=o5)
                                    nReadsinSpecies+=len(v)
                                    
                                else:
                                    print(str(k.split("_")[-1])+"\t"+"_".join(k.split("_")[:-1])+"\t"+str(len(v)), file=o3) # Need to take everything except as some species have more _ than one!
                                    nReadsSpecies_discarded+=len(v)
                            for k, v in EverythingAboveSpecies_withReads.items():
                                print(str(k)+"\t"+ EverythingAboveSpecies[k][1] +"\t"+ EverythingAboveSpecies[k][0] +"\t"+ str(len(v)), file=o2)
                                nReadsHighTaxLevels+=len(v)

                            if not (nReadsSpecies_discarded + nReadsHighTaxLevels + nReadsinSpecies) == ClassifiedReads:
                                print("warning")
                                print(nReadsSpecies_discarded, nReadsHighTaxLevels, nReadsinSpecies)
                                
                            logging.info('%s\tAmount of ClassifiedReads: %s', time.ctime(), ClassifiedReads)
                            logging.info('%s\tAmount of reads towards Species: %s', time.ctime(), nReadsinSpecies)

                            print("nReads total:\t"+str(ntotreads), file=o4)
                            print("nReads classified:\t"+str(ClassifiedReads), file=o4)
                            print("nReads to Species:\t"+str(nReadsinSpecies), file=o4)
                            print("nReads to Species_discarded (due to low counts):\t" + str(nReadsSpecies_discarded), file=o4)
                            print("nReads to Higher ranks:\t" + str(nReadsHighTaxLevels), file=o4)
                            
                            # Extract the reads, outputed from the bowtie directory
                            if Fastqfiles and not IgnoreReadExtraction:
                                outfolderclassifiedreads="Diamond/Classified_Reads/"
                                try:
                                    os.makedirs(outfolderclassifiedreads)
                                except FileExistsError: 
                                    pass
                                if "_se_" in d: 
                                    basename=d.split("/")[-1].split(f'{diamonddb}.diamond.tsv')[0].replace("_se_","_").strip("_") # remove _se_ to be able to link to original fastq
                                if "_pe_" in d: 
                                    basename=d.split("/")[-1].split(f'{diamonddb}.diamond.tsv')[0].replace("_pe_","_").strip("_") # remove _pe_ to be able to link to original fastq
                                for f in Fastqfiles:
                                    if basename in f:
                                        
                                        if f.endswith(".gz"): # If the files are gziped you need to use gzip open, save record to dict and get a basename from the fastq
                                            Records=SeqIO.to_dict(SeqIO.parse(gzip.open(f, "rt"),'fastq'))
                                            fname=f.split("/")[-1].replace(".unmapped","").split(".fastq.gz")[0]
                                        else:
                                            Records=SeqIO.to_dict(SeqIO.parse(f,'fastq'))
                                            fname=f.split("/")[-1].replace(".unmapped","").split(".fastq")[0]
                                        for k, v in Species_withReads.items():
                                            if len(v) >= dptresh: 
                                                taxa="_".join(k.split("_")[:-1]).replace(" ","")+"_"+str(k.split("_")[-1])
                                                outfoldersspecies="Diamond/Classified_Reads/"+taxa
                                                try:
                                                    os.makedirs(outfoldersspecies)
                                                except FileExistsError: 
                                                    pass
                                                outfq=outfoldersspecies+"/"+taxa+"_"+fname+".fastq"
                                                with open(outfq, "w") as o: 
                                                    for reads in v:
                                                        try: 
                                                            print(Records[reads].format("fastq").strip(), file=o)
                                                        except KeyError: # We get key errors we are grepping for R1 in R2.
                                                            continue
    # Remove the tmp where we download the sqlite database
    os.remove(f'{tmp}/taxa.sqlite')
    os.remove(f'{tmp}/taxa.sqlite.traverse.pkl')
    #os.rmdir(tmp)
    
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

def main(taxprofdict, taxdump, dptresh, dbsheet, IgnoreReadExtraction):
    ParseDiamond(taxprofdict, taxdump, dptresh, dbsheet, IgnoreReadExtraction)


if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict, args.taxdump, args.dptresh, args.dbsheet, args.IgnoreReadExtraction)

