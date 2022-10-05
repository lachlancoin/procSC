#!/bin/bash
#SBATCH --job-name=scp_proc
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=7000 #mb
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=2
echo "User: $USER"
echo "Datetime: $(date)"
echo "Current directory: $PWD"
echo "Anaconda environment path: $CONDA_PREFIX"

#bash run_scp.sh  params.json


#module load r-bundle-bioconductor/3.12-r-4.0.4
module load   r-bundle-bioconductor/3.10-r-3.6.2
# Activate an environment:
#source activate R
#rcmd="R CMD BATCH --no-save --no-restore"
#$rcmd procSC.R
#tail procSC.Rout
dat=$(date +%Y%m%d%H%M%S)


Rscript  ~/github/procSC/R/procSC.R $1 

#now prepare the zip output files
#cd $outdir
#find . -mindepth 2 -type d > todo.txt
#while read line; do a=$(basename $line) ; find $line -type f | xargs -I {} zip -R $a.zip {} ; done < todo.txt

#find . -maxdepth 1 -mindepth 1 -type d > todo.txt
#while read line; do a=$(basename $line) ; find $line -type f | xargs -I {} zip -R $a.zip {} ; done < todo.txt

