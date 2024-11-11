#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --job-name=IVSA_SALINE_PLEIO
#SBATCH --chdir=/projects/csna/rta/
#SBATCH --output=/projects/csna/rta/logs/ivsa_sens_pleio.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=michael.saul@jax.org
#SBATCH --account=chesler-lab
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=700GB
#SBATCH --export=ALL
### SLURM HEADER
localcores=${SLURM_CPUS_PER_TASK}
simg_path="/projects/csna/rta/singularity/rstudio_rta.sif"
work_dir="/projects/csna/rta/"
module load singularity
cd ${work_dir}
singularity exec \
    ${simg_path} R --slave -e 'source("/projects/csna/rta/code/DO_IVSA_SALINE_pleiotropy.R")'
