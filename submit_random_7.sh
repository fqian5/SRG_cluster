#!/bin/bash
#SBATCH -p batch
#SBATCH --job-name random_7
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 5 # cores requested
#SBATCH -t 3-00:00:00
#SBATCH --mem=50000 # memory in Mb
#SBATCH -o outputfile_random_7.out # send stdout to outfile
#SBATCH -e errfile_random_7.out  # send stderr to errfile
module load miniforge/24.11.2-py312
source activate /cluster/tufts/lovelab/fqian03/condaenv/gap
python job_random_7.py 1>out_random_7 2>error_random_7
