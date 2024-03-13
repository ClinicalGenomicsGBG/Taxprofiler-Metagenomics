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
    parser = argparse.ArgumentParser(description='Generates an Excel file for the different tools that were run in taxprofiler')
    parser.add_argument("--Tools", dest = 'Tools', required=True, help ="Path to taxprofiler_parsed output directory",nargs="+")
    arguments = parser.parse_args(argv)
    return arguments



def SummarizeInExcel(Tools):
    """
    One sheet per Tool
    """

    OutExcel="Out_ParsedExcel.xlsx"
    writer = pd.ExcelWriter(OutExcel, engine='xlsxwriter')
    Countfiles={}
    for tool in Tools:
        Countfiles[tool]=[]
        countstoplot=glob.glob(tool+"/*_CountsForplotting.txt")
        for c in countstoplot:
            samplename=c.split("/")[-1]
            Countfiles[tool].append(samplename)
    for folder, files in Countfiles.items():
        counter=0
        for f in files:
            counter+=1
            samplename=f.split("_CountsForplotting.txt")[0]
            df=pd.read_csv(folder+"/"+f, sep="\t")
            df=df.rename(columns={'Counts': samplename})
            if counter==1:
                dfMerged=df.copy()
                
            else:
                dfMerged=dfMerged.merge(df, how='outer',on=['TaxID','Species']).fillna(0) # Merge and set NA to 0 
        cols=[i for i in dfMerged.columns if i not in ["TaxID","Species"]]
        for col in cols:
            dfMerged[col]=dfMerged[col].astype(int) # Convert from float to int, will be float as we introduced NAs in the merging
        tool=folder.split("/")[-1]
        dfMerged.to_excel(writer, sheet_name=tool)
        
    writer.close()
                                    
def main(Tools):
    SummarizeInExcel(Tools)
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.Tools)

    
