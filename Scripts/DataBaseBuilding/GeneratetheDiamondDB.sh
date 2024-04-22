#!/bin/bash -l

#$ -N DiamondDB
#$ -j y
#$ -cwd
#$ -pe mpi 6
#$ -q development.q

module load miniconda/4.14.0

source activate diamond

set -x

zcat fastas/*.faa.gz | diamond makedb --taxonmap prot.accession2taxid.FULL.gz --taxonnames taxdump/names.dmp --taxonnodes taxdump/nodes.dmp -d refseq_223_Bacteria_Virus_Archaea_Fungi_diamond -v --log --threads 6

