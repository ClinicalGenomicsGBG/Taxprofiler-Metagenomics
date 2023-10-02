#!/usr/bin/py

import argparse
from Bio import SeqIO
from collections import Counter
import glob
import gzip
import logging
#import numpy as np
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
    parser.add_argument("--Tools", dest = 'toolfolders', required=True, help ="All tools you want to parse together from the nextflow parser (Required), could be [Diamond, Kraken2, KrakenUniq, Metaphlan4], presented in a space seperatated string", nargs="+")
    arguments = parser.parse_args(argv)
    return arguments




def SummarizeInExcel(toolfolders):
    """
    One sheet per sample
    """

    OutExcel="Out_ParsedExcel.xlsx"
    writer = pd.ExcelWriter(OutExcel, engine='xlsxwriter')

    Countfiles={}

    # Generate 1 dataframe per sample! 
    for t in toolfolders:
        Countfiles[t]=[]
        countstoplot=glob.glob(t+"/*_CountsForplotting.txt")
        for c in countstoplot:
            samplename=c.split("/")[-1]
            Countfiles[t].append(samplename)



        
    uniqueCountFiles=list(set(sum(list(Countfiles.values()),[]))) # Get the unique counts for plotting name to be able to loop through them! 




    for c in uniqueCountFiles:
        #df=pd.DataFrame() # For every sample generate a dataframe 
        counter=0
        for tool, reports in Countfiles.items():
            counter+=1
            if c in reports:
                samplename=c.split("_CountsForplotting.txt")[0]
                d={}
                fullpath=tool+"/"+c
                
                with open(fullpath, "r") as inf: # Check amount of rows to make sure you can loop through it! 
                    nrowsCounts=len(inf.readlines())
                if nrowsCounts > 1:
                    with open(fullpath, "r") as inf:
                        next(inf) # Skip headerline
                        for l in inf:
                            l=l.strip()
                            print(l)
                            taxnr=l.split("\t")[0]
                            taxname=l.split("\t")[1]
                            counts=int(l.split("\t")[2])
                            d[taxname]=counts
                    if counter == 1: # First tool for the sample generate the dataframe
                        df=pd.DataFrame.from_dict(d, orient='index', columns=[tool])
                    else: # If we already had a tool append to the existing dataframe
                        df2=pd.DataFrame.from_dict(d, orient='index', columns=[tool])
                        df=pd.merge(df,df2, how='outer', left_index=True, right_index=True)

        df=df.astype('Int64') # To convert columns with NAN to int we use int64
        # Print calculate the detection rate, if 1 in 4 the detection is 25%
        df["DetectionRate"]=(len(df.columns)-df.isna().sum(axis=1))/len(df.columns)
        df=df.sort_values('DetectionRate', ascending=False)

        # We cannot have sheet name longer than 30 characters in excel
        if len([*samplename]) > 30:
            new_samplename="".join([*samplename][0:30])

        else:
            new_samplename=samplename

        # Create a Pandas Excel writer using XlsxWriter as the engine.
        df.to_excel(writer, sheet_name=new_samplename)
    writer.close()
                                    
def main(toolfolders):
    SummarizeInExcel(toolfolders)
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.toolfolders)

    
