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

# Turning off the performance warnings, might want to fix this later! 
from warnings import simplefilter 
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


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
    #OutExcel2="Out_ParsedExcel.xlsx" 
    #writer_2 = pd.ExcelWriter(OutExcel, engine='xlsxwriter')
    workbook=writer.book
    
    Countfiles={}
    Mapping={}
    for tool in Tools:
        detec=[]
        Countfiles[tool]=[]
        countstoplot=glob.glob(tool+"/*_CountsForplotting.txt")
        for c in countstoplot:
            samplename=c.split("/")[-1]
            Countfiles[tool].append(samplename)
        mappingfiles=glob.glob(tool+"/Extras/*_SpeciesDomainLinkage.txt")
        Mapping[tool]=[]
        for m in mappingfiles:
            samplename=m.split("/")[-1]
            with open(m, "r") as inf:
                next(inf)
                for l in inf:
                    l=l.strip()
                    taxid=int(l.split("\t")[0])
                    species=l.split("\t")[1]
                    kingdom=l.split("\t")[2]
                    if not taxid in detec:
                        Mapping[tool].append([taxid, species, kingdom])
                    
        
    for folder, files in Countfiles.items():
        tool=folder.split("/")[-1]
        Mapping_df = pd.DataFrame(Mapping[tool], columns = ['TaxID', 'Species', 'Kingdom']) 
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


        # Move TaxID and species to first and second column!
        # Spits out performance warnings, fix that!
        
        #column_to_move = dfMerged.pop("Species")
        
        #dfMerged.insert(0, "Species", column_to_move)
        #column_to_move = dfMerged.pop("TaxID")
        #dfMerged.insert(0, "TaxID", column_to_move)

        dfMerged.insert(0, 'Species', dfMerged.pop('Species'))
        dfMerged.insert(0, 'TaxID', dfMerged.pop('TaxID'))
        
        #print(columns.index('TaxID'))
        #print(columns.index('Species'))

        #df2 = dfMerged[['TaxID', 'Species']].copy()
        #dfMerged=dfMerged.drop(['TaxID', 'Species'], axis=1)
        #dfMerged=pd.concat([df2,dfMerged], axis=1)

        
        #dfMerged=pd.concat([column_to_move,dfMerged.insert], axis=1)
        
        # Add normalization to the excel sheet also 

        for sample in dfMerged.columns[2:]:
            norm=f'{sample}_Percent'
            #dfMerged[norm]=(dfMerged[sample]+1)/(dfMerged[sample].sum()+len(dfMerged))*100 # Add a pseudocount to all sample and add the amount of rows to the sum
            dfMerged[norm]=dfMerged[sample]/dfMerged[sample].sum()*100
            
        dfMerged_kingdom=dfMerged.merge(Mapping_df, how='left', on=['TaxID','Species']) # Merge with kingdom
        dfMerged_kingdom.drop_duplicates(subset=['TaxID'], keep='first', inplace=True, ignore_index=True) # Drop the duplicated rows from the merge!

        #print(dfMerged_kingdom)
        
        Archaea=dfMerged_kingdom[dfMerged_kingdom['Kingdom']=='Archaea']
        Archaea=Archaea.drop(columns=['Kingdom'])
        dimArchaea = len(Archaea.columns)
        Bacteria=dfMerged_kingdom[dfMerged_kingdom['Kingdom']=='Bacteria']
        Bacteria=Bacteria.drop(columns=['Kingdom'])
        dimBacteria = len(Bacteria.columns)
        Eukaryota=dfMerged_kingdom[dfMerged_kingdom['Kingdom']=='Eukaryota']
        Eukaryota=Eukaryota.drop(columns=['Kingdom'])
        dimEukaryota=len(Eukaryota.columns)
        Viral=dfMerged_kingdom[dfMerged_kingdom['Kingdom']=='Viruses']
        Viral=Viral.drop(columns=['Kingdom'])
        dimViral=len(Viral.columns)
        
        col=0
        row=1
        
        # --- Write to excel ---, seperate them to data tables for excel when possible! 
        Archaea.to_excel(writer, sheet_name=tool,  startrow=row, startcol=col, index=True)        
        worksheet = writer.sheets[tool]
        worksheet.write(0, col, "Archaea")
        if Archaea.shape[0] > 1:
            Archaea.reset_index(inplace=True)
            header = [{'header': di} for di in Archaea.columns.tolist()]
            worksheet.add_table(row, col, Archaea.shape[0], Archaea.shape[1]-1,{'header_row': True,'first_column': True,'columns':header, 'style': 'Table Style Medium 15'})

        col+=2+dimBacteria
        Bacteria.to_excel(writer, sheet_name=tool,  startrow=row, startcol=col, index=True)
        worksheet.write(0, col, "Bacteria")
        if Bacteria.shape[0] >1: 
            Bacteria.reset_index(inplace=True)
            header = [{'header': di} for di in Bacteria.columns.tolist()]
            worksheet.add_table(row, col, Bacteria.shape[0], col+Bacteria.shape[1]-1,{'header_row': True,'first_column': True,'columns':header, 'style': 'Table Style Medium 15'})        
        
        col+=2+dimEukaryota
        Eukaryota.to_excel(writer, sheet_name=tool,  startrow=row, startcol=col, index=True)
        worksheet.write(0, col, "Eukaryota")
        if Eukaryota.shape[0] > 1:
            Eukaryota.reset_index(inplace=True)
            header = [{'header': di} for di in Eukaryota.columns.tolist()]
            worksheet.add_table(row, col, Eukaryota.shape[0], col+Eukaryota.shape[1]-1,{'header_row': True,'columns':header, 'style': 'Table Style Medium 15'})
        
        col+=2+dimViral
        Viral.to_excel(writer, sheet_name=tool,  startrow=row, startcol=col, index=True)
        worksheet.write(0, col, "Viral")
        if Viral.shape[0] > 1: 
            Viral.reset_index(inplace=True)
            header = [{'header': di} for di in Viral.columns.tolist()]
            worksheet.add_table(row, col, Viral.shape[0], col+Viral.shape[1]-1,{'header_row': True,'first_column': True,'columns':header, 'style': 'Table Style Medium 15'})
            
        dfMerged.to_excel(writer, sheet_name=f'{tool}_raw')

    writer.close()
    
def main(Tools):
    SummarizeInExcel(Tools)
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.Tools)

    
