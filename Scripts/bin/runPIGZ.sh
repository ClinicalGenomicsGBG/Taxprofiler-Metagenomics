#!/bin/bash -l

#$ -N pigz
#$ -j y
#$ -cwd
#$ -pe mpi 20
#$ -q development.q


set -x

pigz $1 -p 20


