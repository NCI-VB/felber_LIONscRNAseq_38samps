#!/bin/bash
#SBATCH -J vsAll_rpca_fullwithSubs_tt           # Job name
#SBATCH -n 1                # Run on a single node
#SBATCH --cpus-per-task=48               # number of CPUs
#SBATCH --mem 200G          # memory (RAM) request
#SBATCH -o ccrvb47_18samps_vsAll_rpca_fullwithSubs_tt_out.log        # output (.%j)
#SBATCH --gres=lscratch:500
export TMPDIR=/lscratch/$SLURM_JOB_ID

module purge
module load R/4.4

if [[ ! -d /data/NCI_VB/$USER/R_libs/4.4 ]]; then
		  mkdir /data/NCI_VB/$USER/R_libs/4.4
fi
export R_LIBS_USER="/data/NCI_VB/$USER/R_libs/4.4"


### move to the dir with the R script
cd /data/NCI_VB/goldfarbmurenkc/github/VaccineBranch/Felber/Rcall_scripts/

### call R script 
# study_name_to_use <- as.character(args[1]) ##includes batch correction batch123_38samps_rpca_ann, batch12_18samps_rpca_tiss_trimmed_ann, 
# res_to_use <- as.character(args[2]) #respt4_rpca_broadClusters, rpca_b12_full_res1pt2_withSubs
# tissue_sep ##yes/no


Rscript ccrvb47_LIONscRNAseq_ConservedClusterDEGs_38samps_KCG.R batch12_18samps_rpca_tiss_trimmed_ann rpca_b12_full_res1pt2_withSubs yes