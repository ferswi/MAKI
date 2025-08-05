#!/bin/bash
#SBATCH --job-name=kinase_job_test1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=07:00:00
#SBATCH --account=rrg-najmanov
#SBATCH --mail-user=sophie.bernard-doucet@umontreal.ca
#SBATCH --mail-type=START,END,FAIL

module load python/3.10

source /home/sbernarr/ENV/bin/activate

which python

python MAKI.py