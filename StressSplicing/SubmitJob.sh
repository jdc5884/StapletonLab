#!/bin/bash


#SBATCH -J BreedvQTL	# Job name

#SBATCH -e myjob.e%j    # Name of stderr error file

#SBATCH -p skx-normal   # Queue (partition) name

#SBATCH -N 1            # Total # of nodes (must be 1 for serial)
#SBATCH -n 1            # Total # of mpi tasks (should be 1 for serial)

#SBATCH -t 48:00:00     # Run time (hh:mm:ss)

#SBATCH --mail-user=jdc5884@uncw.edu

#SBATCH --mail-type=all  # Send email at begin and end of job

# Other commands must follow all 

mkdir output
Rscript --vanilla --verbose ./BreedvQTL.R > ./output.Rout
