#!/usr/bin/py

import argparse
import os
import subprocess
import sys
import shutil
from datetime import datetime

import smtplib
from email.message import EmailMessage



def parseArgs(argv):
    '''
    Parsing the arguments
    '''
    parser = argparse.ArgumentParser(description='Takes the output from taxprofiler and parses it')
    parser.add_argument("--FolderToMove", dest = 'folder', required=True, help ="Path to outputfolder to copy to webstore (required)")
    parser.add_argument("--mail", dest = 'mail', action='store_true', help ="Send mail after run is finished (optional)")
    
    arguments = parser.parse_args(argv)
    return arguments

def CreateFolderInWebstore(folder):
    pathtowebstore="/medstore/Development/Metagenomics/TaxProfiler/InsilicoTester_SE/Transfer"
    now=datetime.now()
    time=datetime.strftime(now, '%y%m%d-%H%M%S') # Datetime required for the backup to work
    folderbasename=folder.split("/")[-1] #if full path is given we need the basename for creation
    folderinWebstore=f'{pathtowebstore}/{folderbasename}_{time}'
    os.mkdir(folderinWebstore)
    return(folderinWebstore)


def CopyFiles(folder,folderinWebstore):
    for i in os.listdir(folder):
        if ".xlsx" in i:
            # we are transfering the excels
            command = f'cp {i} {folderinWebstore}'
            shutil.copy2(f'{folder}/{i}',folderinWebstore)
        if 'Taxprofiler_Parsed' in i:
            # We are transferring the entire Parsed directory
            command = f'cp -r {i} {folderinWebstore}'
            os.mkdir(f'{folderinWebstore}/TaxProfiler_Parsed')
            shutil.copytree(f'{folder}/{i}', f'{folderinWebstore}/TaxProfiler_Parsed',dirs_exist_ok=True)
        if 'TaxProfiler_out' in i:
            os.mkdir(f'{folderinWebstore}/TaxProfiler_out')
            for taxoutfo in os.listdir(f'{folder}/{i}'):
                if taxoutfo in ['krona','multiqc','pipeline_info']: # transfer krona, multiqc and pipeline_info from taxprofiler
                    os.mkdir(f'{folderinWebstore}/TaxProfiler_out/{taxoutfo}')
                    shutil.copytree(f'{folder}/TaxProfiler_out/{taxoutfo}', f'{folderinWebstore}/TaxProfiler_out/{taxoutfo}',dirs_exist_ok=True)


def sendingmail(folder):

    folderbasename=folder.split("/")[-1] #if full path is given remove split to get to the basename
    
    subject="Results from TaxProfiler pipeline now on BDC-portal"
    body=f'TaxProfiler results from {folderbasename} is now available on the BDC-portal.\n' \
         f'\n' \
         f'GU link: http://webstore.sa.gu.se/tree/development/Taxprofiler \n'\
         f'SU link: http://webstore-gu.vgregion.se/tree/development/Taxprofiler \n'
    
    msg = EmailMessage()
    msg.set_content(f'{body}\n'
                    f'\n'
                    f'Kind regards,\n'
                    f'CGG micro')

    msg['Subject'] = f'{subject}'
    msg['From'] = "sanna.abrahamsson@gu.se"
    msg['To'] = ["sanna.abrahamsson@gu.se"]
    msg['Cc'] = ["sanna.abrahamsson@gu.se"]

    #Send the messege
    s = smtplib.SMTP('smtp.gu.se')
    s.send_message(msg)
    s.quit()
                    
                    
def main(folder, mail):
    folderinWebstore=CreateFolderInWebstore(folder)
    CopyFiles(folder,folderinWebstore)
    
    if mail:
        sendingmail(folder)
        
if __name__ == '__main__':
    args=parseArgs(sys.argv[1:])
    main(args.folder, args.mail)

