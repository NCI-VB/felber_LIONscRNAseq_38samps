#!/bin/bash

. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
#conda init scRNA4
#conda activate scRNA4

export TMPDIR=/lscratch/$SLURM_JOB_ID

module purge
module load R/4.2

if [[ ! -d /data/NCI_VB/$USER/R_libs/4.2 ]]; then
                  mkdir /data/NCI_VB/$USER/R_libs/4.2
fi
export R_LIBS_USER="/data/NCI_VB/$USER/R_libs/4.2"

module load python

Rscript workflow/scripts/integrationReport.R batch123_38samps_merged.rds batch123_38samps_integrated_cca.rds batch123_38samps_harmonyGroup.rds batch123_38samps_harmonySample.rds batch123_38samps_rpca.rds $PWD/integration/images_batch123_38samps $PWD/integration 0.2,0.4,0.6,0.8,1.0,1.2,1.4 batch123_38samps YES NO
