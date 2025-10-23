#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name random_5
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 5 # cores requested
#SBATCH -t 7-00:00:00
#SBATCH --mem=500000 # memory in Mb
#SBATCH -o outputfile_random_5.out # send stdout to outfile
#SBATCH -e errfile_random_5.out  # send stderr to errfile
module load miniforge/24.11.2-py312
source activate /cluster/tufts/lovelab/fqian03/condaenv/gap
python job_random_5.py 1>out_random_5 2>error_random_5
