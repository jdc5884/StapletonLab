#!/bin/bash

#SBATCH -J vQTL          # Job name
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p skx-normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=jdc5884@uncw.edu
#SBATCH --mail-type=end    # Send email at begin and end of job


# Other commands must follow all #SBATCH directives...

module load Rstats # Load the R module along with some popular packages so it will run R files
Rscript --verbose ./vQTL.R > ./output.Rout
