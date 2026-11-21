#!/bin/bash
#SBATCH -J LN_DAY2vsPRE_cellchatsec_tt            # Job name
#SBATCH -n 1                # Run on a single node
#SBATCH --cpus-per-task=48               # number of CPUs
#SBATCH --mem 200G          # memory (RAM) request
#SBATCH -o cellchatsecDiffComp_batch123_tissTrimed_LN_DAY2vsPRE_out.log        # output (.%j)
#SBATCH --gres=lscratch:500
export TMPDIR=/lscratch/$SLURM_JOB_ID

module purge
module load R/4.3

if [[ ! -d /data/NCI_VB/$USER/R_libs/4.3 ]]; then
		  mkdir /data/NCI_VB/$USER/R_libs/4.3
fi
export R_LIBS_USER="/data/NCI_VB/$USER/R_libs/4.3"


### move to the dir with the R script
cd /data/NCI_VB/goldfarbmurenkc/github/VaccineBranch/Felber/Rcall_scripts/

### call R script (same for LN or PBMC, each contrast)
# study_name_to_use <- as.character(args[1]) ##batch123_38samps_rpca or batch123_38samps_rpca_tiss_trimmed # 
# tissue_to_use <- as.character(args[2]) ##LN or PBMC # 
# res_to_use <- as.character(args[3]) ## respt4_rpca_ClustersWithSubs
# contrast_to_use <- as.character(args[4]) ## DAY2vsPRE, DAY8vsPRE, DAY8vsDAY2
# species <- as.character(args[5]) ## human (or could be mouse)
# both_or_second_only <- as.character(args[6]) ## "both" or "second_only" if had issues with second half



Rscript ccrvb47_CellChatPrep_DiffComp_batch123.R batch123_38samps_rpca_tiss_trimmed LN respt4_rpca_ClustersWithSubs DAY2vsPRE human second_only