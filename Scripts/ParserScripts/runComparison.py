#!/usr/bin/py

import argparse
from Bio import SeqIO
import glob
import logging
import os
import pandas as pd
import psutil
import subprocess
import sys
import time
import xlsxwriter
from datetime import date



def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes the output from Parsed taxprofiler and Calculates to comparison Patient vs Control')
    parser.add_argument("--Parsed_Taxprofiler_out", dest = 'parsedFolder', required=True, help ="Path to the output folder that is parsed from taxprofiler, contains the reads and the count tables")
    parser.add_argument("--Metadata", dest = 'metadata', required=True, help ="csv file with metadata, Sample (samplename), Type (DNA or RNA), Group (P (patient) or C (control))")
    arguments = parser.parse_args(argv)
    return arguments


def ReadMetadata(metadata):

    comparisons={}

    RNAControls=[]
    DNAControls=[]
    RNAPatients=[]
    DNAPatients=[]

    # Reading the Metadata
    with open(metadata, "r") as inf:
        for l in inf:
            l=l.strip()
            if "Sample,Type,Group" in l: 
                continue
            else:
                sample=l.split(",")[0]
                nucl=l.split(",")[1]
                group=l.split(",")[2]
                
                if nucl == "RNA" and group == "C":
                    RNAControls.append(sample)
                if nucl == "RNA" and group == "P":
                    RNAPatients.append(sample)
                if nucl == "DNA" and group == "C":
                    DNAControls.append(sample)
                if nucl == "DNA" and group == "P":
                    DNAPatients.append(sample)

        if RNAControls and RNAPatients:
            comparisons["RNA"]=[RNAControls,RNAPatients]
        if DNAControls and DNAPatients:
            comparisons["DNA"]=[DNAControls,DNAPatients]

    return(comparisons)

def readCountTables(parsedFolder, comparisons):

    subfolders = [ f.path for f in os.scandir(parsedFolder) if f.is_dir() ]


    outExcel="Comparisons.xlsx"
    writer=pd.ExcelWriter(outExcel, engine='xlsxwriter')
    
    for f in subfolders:
        tool=f.split("/")[-1]
        countfiles=glob.glob(f+"/*_CountsForplotting.txt")
        counter=0
        for c in countfiles:
            counter+=1
            sample=c.split("/")[-1].split("_CountsForplotting.txt")[0]
            #norm=sample+"_norm"
            if counter==1: # first df
                df=pd.read_table(c, index_col=None, header=0)
                df.columns=['TaxID','Species',sample]                
            else:
                df2=pd.read_table(c, index_col=None, header=0)
                df2.columns=['TaxID','Species',sample]
                df=df.merge(df2, left_on=['TaxID','Species'], right_on=['TaxID','Species'], how='outer')

        df.fillna(0, inplace=True)
        df[df.columns[2:]]=df[df.columns[2:]]+1

        for sample in df.columns[2:]:
            norm=sample+"_norm"
            df[norm]=df[sample]/df[sample].sum()
            df[sample]=df[sample]-1 # Remove the pseudocount from the table!
            
        normalizedcolumns=[col for col in df.columns if col.endswith('_norm')]
        
        for key, values in comparisons.items():
            sheet=key+"_"+tool
            Controls=values[0]
            Patients=values[1]
            subsetdf=df[['TaxID','Species']].copy()
            for p in Patients:
                for c in Controls:
                    patientcolumns_both=[col for col in df.columns if p in col]
                    controlcolumns_both=[col for col in df.columns if c in col]
                    patientcolumn_norm=[col for col in normalizedcolumns if p in col]
                    controlcolumn_norm=[col for col in normalizedcolumns if c in col]
                    if len(patientcolumn_norm) == 1 and len(controlcolumn_norm) == 1:
                        patientcolumn=patientcolumn_norm[0]
                        controlcolumn=controlcolumn_norm[0]

                        col=patientcolumn+"/"+controlcolumn
                        subsetdf[patientcolumns_both]=df[patientcolumns_both]
                        subsetdf[controlcolumns_both]=df[controlcolumns_both]
                        subsetdf[col]=df[patientcolumn]/df[controlcolumn]


            raw=[]
            norm=[]
            fold=[]
            
            for s in subsetdf.columns:
                if not 'norm' in s:
                    raw.append(s)
                elif not '/' in s:
                    norm.append(s)
                else:
                    fold.append(s)

            ordered=raw+norm+fold     
            subsetdf=subsetdf[ordered]
            
            subsetdf.to_excel(writer, sheet_name=sheet, index=False)

    writer.close()
        
    
def main(parsedFolder, metadata):
    comparisons=ReadMetadata(metadata)
    readCountTables(parsedFolder, comparisons)
    
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.parsedFolder, args.metadata)

