#!/bin/bash

. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
#conda init scRNA4
conda activate scRNA4

export TMPDIR=/lscratch/$SLURM_JOB_ID

module purge
module load R/4.2

if [[ ! -d /data/NCI_VB/$USER/R_libs/4.2 ]]; then
                  mkdir /data/NCI_VB/$USER/R_libs/4.2
fi
export R_LIBS_USER="/data/NCI_VB/$USER/R_libs/4.2"

#module load R/4.2
#module load python

Rscript workflow/scripts/harmony.R filtered/batch123_38samps integration/batch123_38samps_merged.rds integration/batch123_38samps_harmonySample.rds integration/batch123_38samps_harmonyGroup.rds hg38 batch123_38samps 0.2,0.4,0.6,0.8,1.0,1.2,1.4
