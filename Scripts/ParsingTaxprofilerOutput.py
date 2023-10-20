#!/usr/bin/py


# The scripts takes the output folder from taxprofiler and merges the different taxa reports

# Merge with taxpasta,
# Plot per tool
# Compare the tools with venn

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


## OBS OBS

## Error message

"""
qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.
This application failed to start because no Qt platform plugin could be initialized. Reinstalling the application may fix this problem.

Available platform plugins are: eglfs, linuxfb, minimal, minimalegl, offscreen, vnc, wayland-egl, wayland, wayland-xcomposite-egl, wayland-xcomposite-glx, webgl, xcb.

Aborted (core dumped)
"""

#solution with 
#export QT_QPA_PLATFORM=offscreen


# Change taxpasta executable to later!
# /home/xabras/.conda/envs/TaxPasta/bin/taxpasta


def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes the output from taxprofiler and parses it')
    parser.add_argument("--Taxprofiler_out", dest = 'taxprofdict', required=True, help ="Output folder from taxprofiler (required)")
    parser.add_argument("--taxdumpfile", dest = 'taxdump', required=True, help ="Path to taxdump (required)")
    parser.add_argument("--Pathfinder_out", dest = 'pathfinderoutfile', help ="Pathfinder output, (optional). Only there to compare the out from pathfinder to taxprofiler")
    parser.add_argument("--DepthTresh", dest = 'dptresh',default=10,help ="Minimum depth required to be reported")
    parser.add_argument("--Output", dest = 'output', required=True,help ="Output Directory for the parsing (required)")
    arguments = parser.parse_args(argv)
    return arguments

def ParseKraken2(taxprofdict, taxdump, dptresh,output, Outputfolders, Samples):
    """
    
    """

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
                        Outputfolders.append(tool)
                        try:
                            os.mkdir("%s/%s" %(output,tool))
                        except FileExistsError:
                            logging.info('%s\tFolder already exists', time.ctime())
                        # Taxpasta
                        krakdb=k.split("/")[-1]
                        outTaxpasta= output+"/"+tool+"/"+ "kraken2_" + krakdb+".tsv"
                        command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p %s -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*.kraken2.report.txt" %(tool, outTaxpasta, taxdump, k)

                        subprocess.call(command, shell=True)
                        Sampleorder=[]
                        with open(outTaxpasta, "r") as taxpasta:
                            header = taxpasta.readline().split("\t")
                            for i in header[4:]:
                                samplename=i.split(krakdb+".kraken2.kraken2.report")[0]
                                Samples.append(samplename)
                                Sampleorder.append([i,samplename,header.index(i)])
                                
                        for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                            samplename=i[1]
                            sampleindex=i[2]
                            outforplot=output+"/"+tool+"/"+samplename+"CountsForplotting.txt"

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
                            classifiedreads=glob.glob(k+"/*.classifiedreads.txt")
                            for c in classifiedreads:
                                detectedReads={}
                                samplename=c.split(krakdb+".kraken2.kraken2.classifiedreads.txt")[0].split("/")[-1].strip("_")
                                with open(c, "r") as inf: # get the readname from the classified reads in kraken2
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
                                                    outfoldersspecies=output+"/"+tool+"/Classified_Reads/"+key
                                                    try:
                                                        os.makedirs(outfoldersspecies) # Create one output folder per species 
                                                    except FileExistsError: # As we are looping throught the classified reads files there is one per sample, i only want to create one folder per species in the kraken2 folder. If we have one species in more than one sample we need to capture the error. 
                                                        continue
                                                    
                                for f in Fastqfiles: # for all classifiers the samples get the _pe_ addition. check if we have this in the sample name, in that case give a warning, if not just replace pe in the name.
                                    if not "_pe" in f:
                                        snamewithoutpe=samplename.replace("_pe","")
                                        if snamewithoutpe in f:
                                            if f.endswith(".gz"): # If the files are gziped you need to use gzip open, save record to dict and get a basename from the fastq
                                                Records=SeqIO.to_dict(SeqIO.parse(gzip.open(f, "rt"),'fastq'))
                                                fname=f.split("/")[-1].replace(".unmapped","").split(".fastq.gz")[0]
                                            else:
                                                Records=SeqIO.to_dict(SeqIO.parse(f,'fastq'))
                                                fname=f.split("/")[-1].replace(".unmapped","").split(".fastq")[0]
                                            for key, values in detectedReads.items(): # Loop through detected reads, for each species for that sample we are extracting from the fastq file
                                                outfq=output+"/"+tool+"/Classified_Reads/"+key+"/"+key+"_"+fname+".fastq" # Out fastq filename
                                                with open(outfq, "w") as o: 
                                                    for i in values:
                                                        rec=Records[i].format("fastq").strip()                                                        
                                                        print(rec, file=o)                                            
                                    else: 
                                        print("We cannot have pe in the readname, need to think about this!  ")
                                        continue                                                
                                        
                        else:
                            print("Warning, no fastqfiles Available. No read extraction")            


    return(Samples, Outputfolders)


def ParseBracken(taxprofdict, taxdump,dptresh, output, Samples, Outputfolders):
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="bracken"

    
    for i in subfolders:
        if tool in i: # We can extract reads using krak
            subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
            for k in subfolders_2:
                    if "bracken_" in k:
                        Outputfolders.append(tool)
                        try:
                            os.mkdir("%s/%s"%(output,tool))
                        except FileExistsError:
                            logging.info('%s\tFolder already exists', time.ctime())
                        # Taxpasta
                    brackendb=k.split("/")[-1]
                    outTaxpasta= output+"/"+tool+"/"+ "bracken_" + brackendb + ".tsv"
                    command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p %s -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*.%s.tsv" %(tool, outTaxpasta, taxdump, k, tool)
                    subprocess.call(command, shell=True)
                    # Printing required files for plotting
                    Sampleorder=[]
                    with open(outTaxpasta, "r") as taxpasta:
                        header = taxpasta.readline().split("\t")
                        for i in header[4:]:
                            samplename=i.split(brackendb+".bracken")[0]
                            Samples.append(samplename)
                            Sampleorder.append([i,samplename,header.index(i)])
                    for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                        samplename=i[1]
                        sampleindex=i[2]
                        outforplot=output+"/"+tool+"/"+samplename+"CountsForplotting.txt"
                        with open(outTaxpasta,"r") as taxpasta, open(outforplot, "w") as o:
                            print("Taxonomy_nr\tTaxonomy_name\tCounts", file=o)
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
    return(Samples, Outputfolders)



def ParseKaiju(taxprofdict, taxdump,dptresh, output, Samples, Outputfolders):
    """
    """
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="kaiju"
    for i in subfolders:
        subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
        for k in subfolders_2:
            if "Kaiju_" in k:
                Outputfolders.append(tool)
                try:
                    os.mkdir("%s/%s"%(output,tool))
                except FileExistsError:
                    logging.info('%s\tFolder already exists', time.ctime())
                # Taxpasta
                kaijudb=k.split("/")[-1]
                outTaxpasta= output+"/"+tool+"/"+ "kaiju_" + kaijudb + ".tsv"
                command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p kaiju -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*.kaijutable.txt" %(outTaxpasta, taxdump, k)

                subprocess.call(command, shell=True)
                # Printing required files for plotting
                Sampleorder=[]
                with open(outTaxpasta, "r") as taxpasta:
                    header = taxpasta.readline().split("\t")
                    for i in header[4:]:
                        samplename=i.split(kaijudb+".kaijutable")[0]
                        Samples.append(samplename)
                        Sampleorder.append([i,samplename,header.index(i)])
                        
                for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                    samplename=i[1]
                    sampleindex=i[2]
                    outforplot=output+"/"+tool+"/"+samplename+"CountsForplotting.txt"
                    with open(outTaxpasta,"r") as taxpasta, open(outforplot, "w") as o:
                        print("Taxonomy_nr\tTaxonomy_name\tCounts", file=o)
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
    return(Samples, Outputfolders)

def ParseDiamond(taxprofdict, taxdump,dptresh, output, Samples, Outputfolders):
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
                    Outputfolders.append(tool)
                    try:
                         os.mkdir("%s/%s"%(output,tool))
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    # Taxpasta
                    diamonddb=k.split("/")[-1]
                    outTaxpasta= output+"/"+tool+"/"+ "Diamond_" + diamonddb + ".tsv"
                    command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p diamond -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*.diamond.tsv" %(outTaxpasta, taxdump, k)
                    subprocess.call(command, shell=True)
                    # Printing required files for plotting
                    Sampleorder=[]
                    with open(outTaxpasta, "r") as taxpasta:
                        header = taxpasta.readline().split("\t")
                        for i in header[4:]:
                            samplename=i.split(diamonddb+".diamond")[0]
                            Samples.append(samplename)
                            Sampleorder.append([i,samplename,header.index(i)])

                    for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                        Annotation={}
                        samplename=i[1]
                        sampleindex=i[2]
                        outforplot=output+"/"+tool+"/"+samplename+"CountsForplotting.txt"
                        
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
                                            outfoldersspecies=output+"/"+tool+"/Classified_Reads/"+key
                                            try:
                                                os.makedirs(outfoldersspecies) # Create one output folder per species
                                            except FileExistsError: # As we are looping throught the classified reads files there is one per sample, i only want to create one folder per species in the kraken2 folder. If we have one species in more than one sample we need to capture the error.
                                                continue
                            for f in Fastqfiles: # for all classifiers the samples get the _pe_ addition. check if we have this in the sample name, in that case give a warning, if not just replace pe in the name, can we do this in a petter way?
                                if not "_pe" in f:
                                    snamewithoutpe=samplename.replace("_pe","").strip("_") 
                                    #print(snamewithoutpe, f)
                                    if snamewithoutpe in f:
                                        if f.endswith(".gz"): # If the files are gziped you need to use gzip open, save record to dict and get a basename from the fastq
                                            Records=SeqIO.to_dict(SeqIO.parse(gzip.open(f, "rt"),'fastq'))
                                            fname=f.split("/")[-1].replace(".unmapped","").split(".fastq.gz")[0]
                                        else:
                                            Records=SeqIO.to_dict(SeqIO.parse(f,'fastq'))
                                            fname=f.split("/")[-1].replace(".unmapped","").split(".fastq")[0]
                                        for key, values in detectedReads.items(): # Loop through detected reads, for each species for that sample we are extracting from the fastq file
                                            outfq=output+"/"+tool+"/Classified_Reads/"+key+"/"+key+"_"+fname+".fastq" # Out fastq filename
                                            with open(outfq, "w") as o:
                                                for i in values:
                                                    rec=Records[i].format("fastq").strip()
                                                    print(rec, file=o)
                                else:
                                    print("We cannot have pe in the readname, need to think about this!  ")
                                    continue
                        else:
                            print("Warning, no fastqfiles Available. No read extraction")



    return(Samples, Outputfolders)


def ParseMetaphlan3(taxprofdict, taxdump,dptresh, output, Samples, Outputfolders):
    """
    """
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="metaphlan3"
    for i in subfolders:
        subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
        for k in subfolders_2:
            if "Metaphlan3_" in k:
                Outputfolders.append(tool)
                try:
                    os.mkdir("%s/%s"%(output,tool))
                except FileExistsError:
                    logging.info('%s\tFolder already exists', time.ctime())
                # Taxpasta
                metaphlandb=k.split("/")[-1]
                outTaxpasta= output+"/"+tool+"/"+ "Metaphlan3_" + metaphlandb + ".tsv"
                command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p metaphlan -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*_profile.txt" %(outTaxpasta, taxdump, k)
                subprocess.call(command, shell=True)
                # Printing required files for plotting
                Sampleorder=[]
                with open(outTaxpasta, "r") as taxpasta:
                    header = taxpasta.readline().split("\t")
                    for i in header[4:]:                        
                        samplename=i.split(metaphlandb+".metaphlan3_profile")[0]
                        Samples.append(samplename)
                        Sampleorder.append([i,samplename,header.index(i)])
                for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                    samplename=i[1]
                    sampleindex=i[2]
                    outforplot=output+"/"+tool+"/"+samplename+"CountsForplotting.txt"
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
    return(Samples, Outputfolders)

def ParseMetaphlan4(taxprofdict, taxdump,dptresh, output, Samples, Outputfolders):
    """
    """
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    tool="metaphlan"
    for i in subfolders:
        subfolders_2=[ f.path for f in os.scandir(i) if f.is_dir() ] # Check subfolders in kraken2 dir
        for k in subfolders_2:
            if "Metaphlan4" in k:
                Outputfolders.append(tool)
                try:
                    os.mkdir("%s/%s"%(output,tool))
                except FileExistsError:
                    logging.info('%s\tFolder already exists', time.ctime())
                # Taxpasta
                metaphlandb=k.split("/")[-1]
                outTaxpasta= output+"/"+tool+"/"+ "Metaphlan4_" + metaphlandb + ".tsv"
                command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p metaphlan -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*_profile.txt" %(outTaxpasta, taxdump, k)
                subprocess.call(command, shell=True)
                # Printing required files for plotting
                Sampleorder=[]
                with open(outTaxpasta, "r") as taxpasta:
                    header = taxpasta.readline().split("\t")
                    for i in header[4:]:                        
                        samplename=i.split(metaphlandb+".metaphlan_profile")[0]
                        Samples.append(samplename)
                        Sampleorder.append([i,samplename,header.index(i)])
                for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                    samplename=i[1]
                    sampleindex=i[2]
                    outforplot=output+"/"+tool+"/"+samplename+"CountsForplotting.txt"
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
    return(Samples, Outputfolders)



def ParseKrakenUniq(taxprofdict, taxdump,dptresh, output, Samples, Outputfolders):
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
                    Outputfolders.append(tool)
                    try:
                        os.mkdir("%s/%s" %(output,tool))
                    except FileExistsError:
                        logging.info('%s\tFolder already exists', time.ctime())
                    # Taxpasta
                    krakdb=k.split("/")[-1]
                    outTaxpasta= output+"/"+tool+"/"+ "krakenuniq_" + krakdb+".tsv"
                    command="/home/xabras/.conda/envs/TaxPasta/bin/taxpasta merge -p %s -o %s --add-name --add-rank --add-lineage --taxonomy %s %s/*.krakenuniq.report.txt" %(tool, outTaxpasta, taxdump, k)

                    subprocess.call(command, shell=True)
                    Sampleorder=[]
                    with open(outTaxpasta, "r") as taxpasta:
                        header = taxpasta.readline().split("\t")
                        for i in header[4:]:
                            samplename=i.split(".unmapped.krakenuniq.report")[0]
                            Samples.append(samplename)
                            Sampleorder.append([i,samplename,header.index(i)])
                                
                    for i in Sampleorder: # loop through the Taxpasta file ones per sample, save to file!
                        samplename=i[1]
                        sampleindex=i[2]
                        outforplot=output+"/"+tool+"/"+samplename+"CountsForplotting.txt"

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
                                                outfoldersspecies=output+"/"+tool+"/Classified_Reads/"+key
                                                try:
                                                    os.makedirs(outfoldersspecies) # Create one output folder per species 
                                                except FileExistsError: # As we are looping throught the classified reads files there is one per sample, i only want to create one folder per species in the kraken2 folder. If we have one species in more than one sample we need to capture the error. 
                                                    continue
                                                    
                            for f in Fastqfiles: # for all classifiers the samples get the _pe_ addition. check if we have this in the sample name, in that case give a warning, if not just replace pe in the name.
                                if not "_pe" in f:
                                    snamewithoutpe=samplename.replace("_pe","")
                                    if snamewithoutpe in f:
                                        if f.endswith(".gz"): # If the files are gziped you need to use gzip open, save record to dict and get a basename from the fastq
                                            Records=SeqIO.to_dict(SeqIO.parse(gzip.open(f, "rt"),'fastq'))
                                            fname=f.split("/")[-1].replace(".unmapped","").split(".fastq.gz")[0]
                                        else:
                                            Records=SeqIO.to_dict(SeqIO.parse(f,'fastq'))
                                            fname=f.split("/")[-1].replace(".unmapped","").split(".fastq")[0]
                                        for key, values in detectedReads.items(): # Loop through detected reads, for each species for that sample we are extracting from the fastq file
                                            outfq=output+"/"+tool+"/Classified_Reads/"+key+"/"+key+"_"+fname+".fastq" # Out fastq filename
                                            with open(outfq, "w") as o: 
                                                for i in values:
                                                    rec=Records[i].format("fastq").strip()                                                        
                                                    print(rec, file=o)                                            
                                else: 
                                    print("We cannot have pe in the readname, need to think about this!  ")
                                    continue                                                
                                        
                    else:
                        print("Warning, no fastqfiles Available. No read extraction")            

    return(Samples, Outputfolders)
    

def plotting_withinEachTool(output,Outputfolders):
    """
    Here we look at each tool each sample by themselvses which are their top 10
    """
    
    for i in Outputfolders:
        path=output+"/"+i
        for f in glob.glob(path+"/*_CountsForplotting.txt"): # loop through all _CountsForplotting.txt
            d={}
            with open(f, "r") as inf: # Check amount of lines in Counts for plotting, when they dont have any counts they only have the header, then skip the plotting for these! 
                nrowsCounts=len(inf.readlines())
            if nrowsCounts > 1: # We need data to be able to plot!
                with open(f, "r") as inf:
                    next(inf)
                    summa=0
                    for l in inf:
                        l=l.strip()
                        taxid=l.split("\t")[0]
                        taxname=l.split("\t")[1]
                        counts=int(l.split("\t")[2])
                        summa+=counts
                        d[taxid]=[taxname,counts]
                # Check top 10, the ones with the highest fraction            
                df = pd.DataFrame(d)
                df=df.transpose()
                df=df.rename(columns={0:"Taxonomy_name", 1:"Counts"})
                df["Counts"] = df["Counts"].astype('int')
                df['frac']=df['Counts']/df['Counts'].sum()
                df=df.sort_values(by=['Counts'],ascending=False).head(n=10).reset_index(drop=True)
                out=f.split("_CountsForplotting.txt")[0]+"_"+i+"_Top10.pdf"
                plot=(ggplot(df) + aes(y="Counts", x="reorder(Taxonomy_name, -Counts)", fill="Taxonomy_name") + geom_bar(stat="identity", position='dodge', show_legend=False) + theme(axis_text_x = element_text(rotation=90, hjust=1)) + xlab("Species"))            
                plot.save(out, width=8, height=8)


            

def plotting_comapareTools(dptresh, output, Outputfolders,Samples):
    """
    Compare the different tools within one sample, 
    """

    OutExcel=output+"/"+"Out_ParsedExcel.xlsx"
    writer = pd.ExcelWriter(OutExcel, engine='xlsxwriter')

    
    if not len(list(set(list(Counter(Samples).values())))) == 1: # Make sure all samples were run in all tools!
        print("Warning sample missing from taxprofiler output")

        
    for s in list(set(Samples)):
        #d={}
        s=s.strip("_")
        ToParse=[]
        df=pd.DataFrame()
        counter=0
        for tool in Outputfolders:
            print(tool, s)
            counter+=1
            path=output+"/"+tool+"/"+s+"_CountsForplotting.txt"
            if os.path.exists(path):
                d={}
                with open(path, "r") as inf: # Check amount of lines again, if it should be added to excel but empty! 
                    nrowsCounts=len(inf.readlines())
                if nrowsCounts > 1: # We need data to be able to plot!
                    with open(path, "r") as inf:
                        next(inf)
                        for l in inf:
                            l=l.strip()
                            taxnr=l.split("\t")[0]
                            taxname=l.split("\t")[1]
                            counts=int(l.split("\t")[2])
                            if counts>=dptresh:
                                d[taxname]=counts
                if counter==1:
                    df=pd.DataFrame.from_dict(d, orient='index', columns=[tool])
                else:
                    df2=pd.DataFrame.from_dict(d, orient='index', columns=[tool])
                    df=pd.merge(df,df2, how='outer', left_index=True, right_index=True)

        df=df.astype('Int64') # To convert columns with NAN to int we use int64
                
        # Print calculate the detection rate, if 1 in 4 the detection is 25%
        df["DetectionRate"]=(len(df.columns)-df.isna().sum(axis=1))/len(df.columns)
        df=df.sort_values('DetectionRate', ascending=False)


        # We cannot have sheet name longer than 30 characters in excel
        if len([*s]) > 30:
            new_s="".join([*s][0:30])

        else:
            new_s=s

        # Create a Pandas Excel writer using XlsxWriter as the engine.

        df.to_excel(writer, sheet_name=new_s)
    writer.close()

    
                                    
def main(taxprofdict,taxdump,pathfinderoutfile,dptresh,output):
    # Create the logger
    process = psutil.Process(os.getpid())
    try: 
        os.mkdir(output)
    except FileExistsError:
        logging.info('%s\tFolder already exists', time.ctime())
    
    logging.basicConfig(level=logging.INFO)
    start = time.time()
    logging.info('%s\t--Parsing the reports from taxprofiler--', time.ctime())
    
    Outputfolders=[]
    subfolders = [ f.path for f in os.scandir(taxprofdict) if f.is_dir() ]
    Samples=[]

    for i in subfolders:
        if "kraken2" in i:
            logging.info('%s\tDetected Kraken2 run, Parsing ', time.ctime())
            #(Samples,Outputfolders)=ParseKraken2(taxprofdict, taxdump,dptresh, output, Samples, Outputfolders)
        if "bracken" in i:
            logging.info('%s\tDetected Bracken run, Parsing ', time.ctime())
            #(Samples,Outputfolders)=ParseBracken(taxprofdict, taxdump, dptresh, output, Samples, Outputfolders)
        if "kaiju" in i:
            logging.info('%s\tDetected kaiju run, Parsing ', time.ctime())
            #(Samples,Outputfolders)=ParseKaiju(taxprofdict, taxdump, dptresh, output, Samples, Outputfolders)
        if "diamond" in i:
            logging.info('%s\tDetected diamond run, Parsing ', time.ctime())
            #(Samples,Outputfolders)=ParseDiamond(taxprofdict, taxdump, dptresh, output, Samples, Outputfolders)
        if "metaphlan" in i:
            logging.info('%s\tDetected Metaphlan4 run, Parsing ', time.ctime())
            (Samples,Outputfolders)=ParseMetaphlan4(taxprofdict, taxdump, dptresh, output, Samples, Outputfolders)
        if "krakenuniq" in i:
            logging.info('%s\tDetected krakenuniq run, Parsing ', time.ctime())
            (Samples,Outputfolders)=ParseKrakenUniq(taxprofdict, taxdump, dptresh, output, Samples, Outputfolders)
            
    logging.info('%s\t--Plotting within Tools--', time.ctime())
    plotting_withinEachTool(output,Outputfolders)
    logging.info('%s\tParsing Finished', time.ctime())
    plotting_comapareTools(dptresh,output, Outputfolders,Samples)

    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.taxprofdict,args.taxdump, args.pathfinderoutfile, args.dptresh, args.output)

