#!/bin/bash
#SBATCH -J Time_rpca_withSubs_TT           # Job name
#SBATCH -n 1                # Run on a single node
#SBATCH --cpus-per-task=48               # number of CPUs
#SBATCH --mem 200G          # memory (RAM) request
#SBATCH -o ccrvb47_38samps_Time_rpca_withSubs_TT_out.log        # output (.%j)
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
# study_name_to_use <- as.character(args[1]) ##includes batch correction batch123_38samps_rpca_tiss_trimmed_ann, ccrvb47_LIONscRNAseq_batch123_38samps_SubclusteringBroad_RPCA_Bcellonly_reSCT
# res_to_use <- as.character(args[2]) #respt4_rpca_broadClusters, respt4_rpca_ClustersWithSubs, rpca_SubBcellClusters_respt4, rpca_SubLymphClusters_respt6, respt4_rpca_broadClusters_myeonly
# tissue_sep ##yes/no


Rscript ccrvb47_LIONscRNAseq_ConservedTimeDEGs_38samps_rpca_KCG.R batch123_38samps_rpca_tiss_trimmed_ann respt4_rpca_ClustersWithSubs yes