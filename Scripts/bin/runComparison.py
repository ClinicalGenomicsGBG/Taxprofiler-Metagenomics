#!/usr/bin/env python


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

pd.options.display.max_columns = None
pd.options.display.max_rows = 100



def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes the output from Parsed taxprofiler and Calculates to comparison Patient vs Control')
    parser.add_argument("--Tools", dest = 'Tools', required=True, help ="Path to the tool folders within the taxprofiler parsed output directory", nargs="+")
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

def readCountTables(Tools, comparisons):

    #subfolders = [ f.path for f in os.scandir(parsedFolder) if f.is_dir() ]
    outExcel="Comparisons.xlsx"
    writer=pd.ExcelWriter(outExcel, engine='xlsxwriter')

    
    Mapping={}

    for f in Tools:
         tool=f.split("/")[-1]
         mappingfiles=glob.glob(f'{f}/Extras/*_SpeciesDomainLinkage.txt')
         Mapping[tool]=[]
         detec=[]
         for m in mappingfiles:
            with open(m, "r") as inf:
                next(inf)
                for l in inf:
                    l=l.strip()
                    taxid=int(l.split("\t")[0])
                    species=l.split("\t")[1]
                    kingdom=l.split("\t")[2]
                    if not taxid in detec:
                        Mapping[tool].append([taxid, species, kingdom])
                        detec.append(taxid)

                        
    for f in Tools:
        tool=f.split("/")[-1]
        countfiles=glob.glob(f'{f}/*_CountsForplotting.txt')
        
        #counter=0

        AllPossibleTypes=['RNA','DNA']

        Type=[]
        
        for t in list(comparisons.keys()): # Check for if we have only RNA or only DNA
            if t in AllPossibleTypes:
                Type.append(t)

        
        Mapping_df=pd.DataFrame(Mapping[tool], columns = ['TaxID','Species','Kingdom']) 
            
        for t in Type:
            counter=0
            for c in countfiles:
                if "_se_" in c: 
                    sample=c.split("/")[-1].split("_se_")[0] # Splitting like this and you remove run id! 
                elif "_pe_" in c:
                    sample=c.split("/")[-1].split("_pe_")[0] # Splitting like this and you remove run id!
                else: # KrakenUniq has _RNA_CountsForplotting.txt or _DNA_CountsForplotting.txt, now... It only has that 
                    if "_RNA_CountsForplotting.txt" in c.split("/")[-1]:
                        sample=c.split("/")[-1].split("_RNA_CountsForplotting.txt")[0]
                    elif "_DNA_CountsForplotting.txt" in c.split("/")[-1]:
                        sample=c.split("/")[-1].split("_DNA_CountsForplotting.txt")[0]
                    else: # We have named it something else, but we dont know the run id!
                        sample=c.split("/")[-1].split("_CountsForplotting.txt")[0]
                        sample = "\t".join(sample.rsplit("_", 1)).split("\t")[0] # Remove everything after the last occurence of _
                if sample in [x for xs in comparisons[t] for x in xs]: # Flatten the list
                    #print(f, t, c)
                    counter+=1
                    df=pd.read_table(c, index_col=None, header=0)
                    df.columns=['TaxID','Species',sample]
                    norm=f'{sample}_Percent'
                    df[norm]=(df[sample]+1)/(df[sample].sum()+len(df))*100 # Perform normalization directly instead at the matrix
                    if counter==1: # first df
                        dfMerged=df.copy()

                    else: 
                        dfMerged=dfMerged.merge(df, how='outer',on=['TaxID','Species']).fillna(0) # Merge and set NA to 0

                    #else:
                    #    df2=pd.read_table(c, index_col=None, header=0)
                    #    df2.columns=['TaxID','Species',sample]
                     #   df=df.merge(df2, left_on=['TaxID','Species'], right_on=['TaxID','Species'], how='outer')


            for col in dfMerged.columns:
                if not col.endswith("_Percent") and not col == "TaxID" and not col == "Species":
                    dfMerged[col]=dfMerged[col].astype(int) # Convert from float to int, will be float as we introduced NAs in the merging

            #print(t,f)
            #print(dfMerged)
                     
            #df.fillna(0, inplace=True)
            #df[df.columns[2:]]=df[df.columns[2:]]+1 # add the + to avoid division with the 0s 
            
            #for sample in df.columns[2:]: # Normalization is the value if percentage of total
             #   norm=f'{sample}_Percent'
             #   df[norm]=df[sample]/df[sample].sum()*100
             #   df[sample]=df[sample]-1


            #if t == "RNA":
            #    if f == "Kraken2":
            #        print(dfMerged)
             
            #print(dfMerged.columns)
             
            normalizedcolumns=[col for col in dfMerged.columns if col.endswith('_Percent')]
            #countcolumns=[col for col in dfMerged.columns if not col.endswith('_Percent')]
            # Write the excel report

            #print(normalizedcolumns)

            
            sheet=f'{t}_{tool}'
            sheet_old=f'{t}_{tool}_old' # This is only kept for the validation
            
            Controls=comparisons[t][0]
            Patients=comparisons[t][1]
            subsetdf=dfMerged[['TaxID','Species']].copy()
            
            for p in Patients:
                for c in Controls:
                    patientcolumns_both=[col for col in dfMerged.columns if p in col]
                    controlcolumns_both=[col for col in dfMerged.columns if c in col]
                    patientcolumn_norm=[col for col in normalizedcolumns if p in col]
                    controlcolumn_norm=[col for col in normalizedcolumns if c in col]
                    if len(patientcolumn_norm) == 1 and len(controlcolumn_norm) == 1:
                        patientcolumn=patientcolumn_norm[0]
                        controlcolumn=controlcolumn_norm[0]
                        
                        col=f'{patientcolumn}+1/{controlcolumn}+1'
                        #patientcolumn+"/"+controlcolumn

                        print(patientcolumns_both)
                        
                        subsetdf[patientcolumns_both]=dfMerged[patientcolumns_both]
                        subsetdf[controlcolumns_both]=dfMerged[controlcolumns_both]

                        
                        
                        subsetdf[col]=(dfMerged[patientcolumn]+1)/(dfMerged[controlcolumn]+1)


            raw=[]
            norm=[]
            fold=[]

                        
            # Reporder so first is the raw values, next if the normalized values and last is the foldchange
            
            for s in subsetdf.columns:
                print(s)
                if not 'Percent' in s:
                    raw.append(s)
                elif not '/' in s:
                    norm.append(s)
                else:
                    fold.append(s)


            
                    
            ordered=raw+norm+fold     
            subsetdf=subsetdf[ordered]

            
            # Split based on kingdom
            
            dfMerged_kingdom=subsetdf.merge(Mapping_df, how='left', on=['TaxID','Species']) 

            Archaea=dfMerged_kingdom[dfMerged_kingdom['Kingdom']=='Archaea']
            Archaea=Archaea.drop(columns=['Kingdom'])
            Archaea=Archaea.sort_values(by=Archaea.columns[-1], ascending=False)
            dimArchaea = len(Archaea.columns)

            Bacteria=dfMerged_kingdom[dfMerged_kingdom['Kingdom']=='Bacteria']
            Bacteria=Bacteria.drop(columns=['Kingdom'])
            Bacteria=Bacteria.sort_values(by=Bacteria.columns[-1], ascending=False)
            dimBacteria = len(Bacteria.columns)
            Eukaryota=dfMerged_kingdom[dfMerged_kingdom['Kingdom']=='Eukaryota']
            Eukaryota=Eukaryota.drop(columns=['Kingdom'])
            Eukaryota=Eukaryota.sort_values(by=Eukaryota.columns[-1], ascending=False)
            dimEukaryota=len(Eukaryota.columns)
            Viral=dfMerged_kingdom[dfMerged_kingdom['Kingdom']=='Viruses']
            Viral=Viral.drop(columns=['Kingdom'])
            Viral=Viral.sort_values(by=Viral.columns[-1], ascending=False)
            dimViral=len(Viral.columns)
        
            col=0
            row=1
        
            # --- Write to excel ---, seperate them to data tables for excel when possible!
            
            Archaea.to_excel(writer, sheet_name=sheet,  startrow=row, startcol=col, index=False)
            
            worksheet = writer.sheets[sheet]
            worksheet.write(0, col, "Archaea")
            if Archaea.shape[0] > 0:
                #Archaea.reset_inde0x(inplace=True)
                header = [{'header': di} for di in Archaea.columns.tolist()]
                worksheet.add_table(row, col, Archaea.shape[0]+1, Archaea.shape[1]-1,{'header_row': True,'first_column': True,'columns':header, 'style': 'Table Style Medium 15'})

            col+=2+dimBacteria
            Bacteria.to_excel(writer, sheet_name=sheet, startrow=row, startcol=col, index=False)
            worksheet.write(0, col, "Bacteria")
            if Bacteria.shape[0] >0: 
                #Bacteria.reset_index(inplace=True)
                header = [{'header': di} for di in Bacteria.columns.tolist()]
                worksheet.add_table(row, col, Bacteria.shape[0]+1, col+Bacteria.shape[1]-1,{'header_row': True,'first_column': True,'columns':header, 'style': 'Table Style Medium 15'})        
        
            col+=2+dimEukaryota
            Eukaryota.to_excel(writer, sheet_name=sheet,  startrow=row, startcol=col, index=False)
            worksheet.write(0, col, "Eukaryota")
            if Eukaryota.shape[0] > 0:
                #Eukaryota.reset_index(inplace=True)
                header = [{'header': di} for di in Eukaryota.columns.tolist()]
                worksheet.add_table(row, col, Eukaryota.shape[0]+1, col+Eukaryota.shape[1]-1,{'header_row': True,'first_column': True, 'columns':header, 'style': 'Table Style Medium 15'})
        
            col+=2+dimViral
            Viral.to_excel(writer, sheet_name=sheet,  startrow=row, startcol=col, index=False)
            worksheet.write(0, col, "Viral")
            if Viral.shape[0] > 0: 
                #Viral.reset_index(inplace=True)
                header = [{'header': di} for di in Viral.columns.tolist()]
                worksheet.add_table(row, col, Viral.shape[0]+1, col+Viral.shape[1]-1,{'header_row': True,'first_column': True,'columns':header, 'style': 'Table Style Medium 15'})
            
                        
            #subsetdf.to_excel(writer, sheet_name=sheet_old, index=False)

    writer.close()
        
    
def main(Tools, metadata):
    comparisons=ReadMetadata(metadata)
    readCountTables(Tools, comparisons)
    
    
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.Tools, args.metadata)

