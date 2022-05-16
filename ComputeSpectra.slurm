#!/bin/bash
#
#SBATCH --job-name=BispectrumTemplates
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ophilcox@princeton.edu
#SBATCH --cpus-per-task=28
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=1gb
#SBATCH --array=1-47
#SBATCH -N 1
#SBATCH --output=/home/ophilcox/out/BispectrumTemplates_%a.log

# SLURM Submission script for the one-loop bispectrum calculations

# Load modules
module load mathematica/11.0.1

# Compute the full one-loop bispectrum templates
echo "Computing full spectra for angular index $SLURM_ARRAY_TASK_ID"
./Spectra.wls $SLURM_ARRAY_TASK_ID -print all

# Compute the flat one-loop bispectrum templates, if required
if [ $SLURM_ARRAY_TASK_ID -lt 8 ]
then
	echo "Computing flattened spectra for angular index $SLURM_ARRAY_TASK_ID"
	./SpectraFlat.wls $SLURM_ARRAY_TASK_ID -print all
fi

echo "Calculations complete!"
