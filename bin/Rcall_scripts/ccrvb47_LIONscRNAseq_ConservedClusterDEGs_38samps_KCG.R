
############## 1vsAll DEGs conserved across 6 animals, for each tissue separately (or combined) ##########



######## get args from shell script
args <- commandArgs(trailingOnly = TRUE)

study_name_to_use <- as.character(args[1]) ##includes batch correction ccrvb47_LIONscRNAseq_batch123_38samps_rpca_ann, ccrvb47_LIONscRNAseq_batch123_38samps_SubclusteringBroad_RPCA_Bcellonly_reSCT
res_to_use <- as.character(args[2]) #respt4_rpca_broadClusters, rpca_SubBcellClusters_respt4
tissue_sep <- as.character(args[3]) ## yes or no


############################################################
# packages & functions
############################################################
# package list
# list.of.packages=c("openxlsx","readxl", "tidyverse","cowplot", "pheatmap","ggpubr","eulerr", "enrichR","knitr","broom","kableExtra","gtsummary","org.Hs.eg.db","Matrix",
#                    "Seurat","SingleR","SingleCellExperiment",
#                    "gridExtra","DT","downloadthis","ggbeeswarm","ggrepel","plotly","ggplotify","flextable","PCAtools","officer","janitor","dittoSeq",
#                    "DESeq2","RColorBrewer","colorspace","filesstrings","geepack","glmGamPoi","sctransform","paletteer","R.utils","data.table", "devtools","ggridges")
list.of.packages=c(
  "Seurat","openxlsx","readxl","tidyverse", "ggpubr","strex")

#install as needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

## MAST directly from github
devtools::install_github("RGLab/MAST")


# load packages
invisible(lapply(list.of.packages, library, character.only = TRUE))
library(MAST)


##### to be run on biowulf, called by ccrvb47_LIONscRNAseq_ConservedClusterDEGs_LN_KCG.sh & ccrvb47_LIONscRNAseq_ConservedClusterDEGs_PBMC_KCG.sh
params <- list(
  "study_name" = study_name_to_use
)

input_dir = "/data/NCI_VB/felber/ccrvb47_LIONscRNAseq_KCG/RData" ### the same as output for running on biowulf, be sure the ann R objects are up to date here
outputRData_dir = "/data/NCI_VB/felber/ccrvb47_LIONscRNAseq_KCG/RData"
outputTABLES_dir = "/data/NCI_VB/felber/ccrvb47_LIONscRNAseq_KCG/TABLES"
code_path = "/data/NCI_VB/goldfarbmurenkc/github/VaccineBranch/Felber/ccrvb47_LIONscRNAseq_KCG.Rmd" ### for reference only, not used


######## load seurat object
seuratObj_to_use_full <- readRDS(file.path(input_dir,paste0(params$study_name,".rds")))
  
  if(grepl("Subclustering",params$study_name)){

    ##### per https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette, merged won't recorrect and gives error if recorrect_umi=T which is default
    recorrect_umi_to_use <- FALSE
    
  }else{
    recorrect_umi_to_use <- TRUE
    
    
    
    # 
    # #### add additional factorized metadata
    # seuratObj_to_use_full$Tissue <- factor(seuratObj_to_use_full$Tissue, c("LN","PBMC"))
    # 
    # seuratObj_to_use_full$Animal_tissue_time <- factor(paste0(as.character(seuratObj_to_use_full$Animal),"_",as.character(seuratObj_to_use_full$Tissue),"_",as.character(seuratObj_to_use_full$Time)),
    #                                               c("MB59_LN_PRE","MC90_LN_PRE", ### combines method3+1 = 2 total samples for each
    #                                                 "MC98_LN_PRE", ## only method1
    #                                                 "MF30_LN_PRE", ## only method3
    #                                                 
    #                                                 "MB59_LN_DAY2", ## only method1
    #                                                 "MC90_LN_DAY2","MC98_LN_DAY2", ### combines method3+1 = 2 total samples for each
    #                                                 "MF30_LN_DAY2", ## only method3
    #                                                 
    #                                                 "MB59_LN_DAY8","MC90_LN_DAY8","MC98_LN_DAY8",#"MF30_LN_DAY8", ## only method3
    #                                                 
    #                                                 "MB59_PBMC_PRE",##2 replicate method3 samples combined + 1 method1 = 3 total samples
    #                                                 "MC90_PBMC_PRE", ### combines method3+1 = 2 total samples
    #                                                 "MC98_PBMC_PRE", ### only method1
    #                                                 "MC73_PBMC_PRE","MF07_PBMC_PRE","MF30_PBMC_PRE",## only method2
    #                                                 
    #                                                 "MB59_PBMC_DAY2","MC90_PBMC_DAY2", ### combines method3+1 = 2 total samples for each
    #                                                 "MC98_PBMC_DAY2", ## only method1
    #                                                 "MC73_PBMC_DAY2", ## combines method3+2 = 2 total samples
    #                                                 "MF07_PBMC_DAY2","MF30_PBMC_DAY2", ## only method2
    #                                                 
    #                                                 
    #                                                 "MB59_PBMC_DAY8","MC90_PBMC_DAY8", ## only method 3
    #                                                 "MC73_PBMC_DAY8","MF07_PBMC_DAY8","MF30_PBMC_DAY8" ## only method 3
    #                                               ))
    # 
    # ###### set names for broad clusters
    # ### for testing
    # # seuratObj_to_use_full <- data.frame("SLM_int_snn_res.0.4"=as.character(seq(0,30,1))) %>%
    # #   dplyr::mutate(SLM_int_snn_res.0.4 = factor(SLM_int_snn_res.0.4,sort(SLM_int_snn_res.0.4)))
    # 
    # #### first relevel
    # seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4 <- factor(seuratObj_to_use_full$SLM_int_snn_res.0.4, as.character(seq(0,length(unique(seuratObj_to_use_full$SLM_int_snn_res.0.4))-1,1)))
    # 
    # 
    # seuratObj_to_use_full$respt4_rpca_broadClusters <- NA
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "0")] <- "CD4 T F0"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "1")] <- "CD4 T naive F1"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "2")] <- "CD8 T naive F2"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "3")] <- "CD8 T F3"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "4")] <- "myeloid CL F4"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "5")] <- "CD8 T F5"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "6")] <- "B cell naive F6"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "7")] <- "CD8 T F7"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "8")] <- "B cell nonsw F8"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "9")] <- "CD4 T F9"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "10")] <- "CD4 T Treg F10"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "11")] <- "CD8 T F11"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "12")] <- "DN T F12"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "13")] <- "B cell sw F13"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "14")] <- "T NK F14"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "15")] <- "Neutrophil F15"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "16")] <- "Erythro F16"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "17")] <- "Platelet F17"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "18")] <- "DN T F18"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "19")] <- "CD8 T F19"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "20")] <- "myeloid IN NC F20"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "21")] <- "CD8 T prolif F21"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "22")] <- "B cell naive F22"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "23")] <- "B cell naive prolif F23"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "24")] <- "DN T B F24"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "25")] <- "NK F25"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "26")] <- "myeloid mDC F26"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "27")] <- "myeloid pDC F27"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "28")] <- "prog F28"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "29")] <- "B cell plasmablast F29"
    # seuratObj_to_use_full$respt4_rpca_broadClusters[which(as.character(seuratObj_to_use_full$relevel_SLM_int_snn_res.0.4) == "30")] <- "GATA3 F30"
    # #### factorize
    # seuratObj_to_use_full$respt4_rpca_broadClusters <- factor(seuratObj_to_use_full$respt4_rpca_broadClusters, c(
    #   "B cell naive F6","B cell naive F22","B cell nonsw F8","B cell sw F13","B cell naive prolif F23","B cell plasmablast F29",
    #   "CD4 T naive F1","CD4 T Treg F10","CD4 T F0","CD4 T F9",
    #   "CD8 T naive F2","CD8 T F3","CD8 T F5","CD8 T F7","CD8 T F11","CD8 T F19","CD8 T prolif F21",
    #   "DN T F12","DN T F18","T NK F14","NK F25","DN T B F24",
    #   
    #   "myeloid CL F4","myeloid IN NC F20","myeloid mDC F26","myeloid pDC F27",
    #   "Neutrophil F15","Platelet F17","Erythro F16","prog F28","GATA3 F30"
    #   
    # ))
    # 
    
  }
  


test_to_use <- "MAST"
assay_to_use <- "SCT"
minpct_to_use <- 0.1
logfc_to_use <- 0.25


##### needs multtest and metap libs
#https://hbctraining.github.io/scRNA-seq_online/lessons/09_merged_SC_marker_identification.html
### can still use SCT now post prepSCTfindmarkers
#https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
#Idents(seuratObj_to_use_full) <- paste0(res_to_use,"nospaces")
subset_res_to_use <- res_to_use
Idents(seuratObj_to_use_full) <- subset_res_to_use

## with SCT
DefaultAssay(seuratObj_to_use_full)<-"SCT"

##### if tissue separately
if(tissue_sep == "yes"){
  tissues_to_use <- list("LN"="LN","PBMC"="PBMC")
}else{
  tissues_to_use <- list("LNandPBMCtogether"=c("LN","PBMC"))
}

print(tissue_sep)
print(tissues_to_use)

########### for each tisue (or both togehter)
list_of_consMarkers <- list()
consMarkers_df <- data.frame()
for(curr_tissue in names(tissues_to_use)){
  
  ### trim to tissues
  seuratObj_to_use <- subset(seuratObj_to_use_full, subset = Tissue %in% tissues_to_use[[paste0(curr_tissue)]])
  
  ### refactorize animal
  seuratObj_to_use$Animal <- factor(seuratObj_to_use$Animal, c("MB59","MC90","MC98", "MC73","MF07","MF30"))

seuratObj_to_use <- PrepSCTFindMarkers(object = seuratObj_to_use)
#SCTmarkers <- FindAllMarkers(seuratObj_to_use, only.pos = TRUE, min.pct = 0.1, thresh.use = 0.25)
#SCTmarkers25 <- FindAllMarkers(seuratObj_to_use, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) #defaults used by nidap, thresh.use is deprecated https://github.com/satijalab/seurat/issues/273




### for testing if missing an animal 
groups_df <- as.data.frame.matrix(table(seuratObj_to_use[[c(subset_res_to_use,"Animal")]], useNA = "ifany")) %>%
  dplyr::mutate(num_groups_zero = ifelse(MB59 < 3 & MC90 < 3 & MC98 < 3 & MC73 < 3 & MF07 < 3 & MF30 < 3,
                                         "all",
                                         ifelse(MB59 < 3 & MC90 < 3 & MC98 < 3 & MC73 < 3 & MF07 < 3 & MF30 >= 3 |
                                                  MB59 < 3 & MC90 < 3 & MC98 < 3  & MC73 < 3 & MF07 >= 3 & MF30 < 3 |
                                                  MB59 < 3 & MC90 < 3 & MC98 < 3  & MC73 >= 3 & MF07 < 3 & MF30 < 3 |
                                                  MB59 < 3 & MC90 < 3 & MC98 >= 3  & MC73 < 3 & MF07 < 3 & MF30 < 3 |
                                                  MB59 < 3 & MC90 >= 3 & MC98 < 3  & MC73 < 3 & MF07 < 3 & MF30 < 3 |
                                                  MB59 >= 3 & MC90 < 3 & MC98 < 3  & MC73 < 3 & MF07 < 3 & MF30 < 3,
                                                "fiveAnimals",
                                                ifelse(MB59 < 3 & MC90 < 3 & MC98 < 3 & MC73 < 3 & MF07 >= 3 & MF30 >= 3 |
                                                         MB59 < 3 & MC90 < 3 & MC98 < 3 & MC73 >= 3 & MF07 < 3 & MF30 >= 3 |
                                                         MB59 < 3 & MC90 < 3 & MC98 >= 3 & MC73 < 3 & MF07 < 3 & MF30 >= 3 |
                                                         MB59 < 3 & MC90 >= 3 & MC98 < 3 & MC73 < 3 & MF07 < 3 & MF30 >= 3 |
                                                         MB59 >= 3 & MC90 < 3 & MC98 < 3 & MC73 < 3 & MF07 < 3 & MF30 >= 3 |
                                                         
                                                         MB59 < 3 & MC90 < 3 & MC98 < 3  & MC73 >= 3 & MF07 >= 3 & MF30 < 3 |
                                                         MB59 < 3 & MC90 < 3 & MC98 >= 3  & MC73 < 3 & MF07 >= 3 & MF30 < 3 |
                                                         MB59 < 3 & MC90 >= 3 & MC98 < 3  & MC73 < 3 & MF07 >= 3 & MF30 < 3 |
                                                         MB59 >= 3 & MC90 < 3 & MC98 < 3  & MC73 < 3 & MF07 >= 3 & MF30 < 3 |
                                                         
                                                         MB59 < 3 & MC90 < 3 & MC98 >= 3  & MC73 >= 3 & MF07 < 3 & MF30 < 3 |
                                                         MB59 < 3 & MC90 >= 3 & MC98 < 3  & MC73 >= 3 & MF07 < 3 & MF30 < 3 |
                                                         MB59 >= 3 & MC90 < 3 & MC98 < 3  & MC73 >= 3 & MF07 < 3 & MF30 < 3 |
                                                         
                                                         MB59 < 3 & MC90 >= 3 & MC98 >= 3  & MC73 < 3 & MF07 < 3 & MF30 < 3 |
                                                         MB59 >= 3 & MC90 < 3 & MC98 >= 3  & MC73 < 3 & MF07 < 3 & MF30 < 3 |
                                                         
                                                         MB59 >= 3 & MC90 >= 3 & MC98 < 3  & MC73 < 3 & MF07 < 3 & MF30 < 3,
                                                       "fourAnimals",
                                                       ifelse(MB59 < 3 & MC90 < 3 & MC98 < 3 & MC73 >= 3 & MF07 >= 3 & MF30 >= 3 |
                                                                MB59 < 3 & MC90 < 3 & MC98 >= 3 & MC73 < 3 & MF07 >= 3 & MF30 >= 3 |
                                                                MB59 < 3 & MC90 >= 3 & MC98 < 3 & MC73 < 3 & MF07 >= 3 & MF30 >= 3 |
                                                                MB59 >= 3 & MC90 < 3 & MC98 < 3 & MC73 < 3 & MF07 >= 3 & MF30 >= 3 |
                                                                
                                                                MB59 < 3 & MC90 < 3 & MC98 >= 3 & MC73 >= 3 & MF07 < 3 & MF30 >= 3 |
                                                                MB59 < 3 & MC90 >= 3 & MC98 < 3 & MC73 >= 3 & MF07 < 3 & MF30 >= 3 |
                                                                MB59 >= 3 & MC90 < 3 & MC98 < 3 & MC73 >= 3 & MF07 < 3 & MF30 >= 3 |
                                                                
                                                                MB59 < 3 & MC90 >= 3 & MC98 >= 3 & MC73 < 3 & MF07 < 3 & MF30 >= 3 |
                                                                MB59 >= 3 & MC90 < 3 & MC98 >= 3 & MC73 < 3 & MF07 < 3 & MF30 >= 3 |
                                                                
                                                                MB59 >= 3 & MC90 >= 3 & MC98 < 3 & MC73 < 3 & MF07 < 3 & MF30 >= 3 |
                                                                
                                                                MB59 < 3 & MC90 < 3 & MC98 >= 3 & MC73 >= 3 & MF07 >= 3 & MF30 < 3 |
                                                                MB59 < 3 & MC90 >= 3 & MC98 < 3 & MC73 >= 3 & MF07 >= 3 & MF30 < 3 |
                                                                MB59 >= 3 & MC90 < 3 & MC98 < 3 & MC73 >= 3 & MF07 >= 3 & MF30 < 3 |
                                                                
                                                                MB59 < 3 & MC90 >= 3 & MC98 >= 3 & MC73 < 3 & MF07 >= 3 & MF30 < 3 |
                                                                MB59 >= 3 & MC90 < 3 & MC98 >= 3 & MC73 < 3 & MF07 >= 3 & MF30 < 3 |
                                                                
                                                                MB59 >= 3 & MC90 >= 3 & MC98 < 3 & MC73 < 3 & MF07 >= 3 & MF30 < 3 |
                                                                
                                                                MB59 < 3 & MC90 >= 3 & MC98 >= 3 & MC73 >= 3 & MF07 < 3 & MF30 < 3 |
                                                                MB59 >= 3 & MC90 < 3 & MC98 >= 3 & MC73 >= 3 & MF07 < 3 & MF30 < 3 |
                                                                MB59 >= 3 & MC90 >= 3 & MC98 >= 3 & MC73 < 3 & MF07 < 3 & MF30 < 3,
                                                              "threeAnimals",
                                                              ifelse(MB59 >= 3 & MC90 >= 3 & MC98 >= 3 & MC73 >= 3 & MF07 < 3 & MF30 < 3 |
                                                                       MB59 >= 3 & MC90 >= 3 & MC98 >= 3 & MC73 < 3 & MF07 >= 3 & MF30 < 3 |
                                                                       MB59 >= 3 & MC90 >= 3 & MC98 < 3 & MC73 >= 3 & MF07 >= 3 & MF30 < 3 |
                                                                       MB59 >= 3 & MC90 < 3 & MC98 >= 3 & MC73 >= 3 & MF07 >= 3 & MF30 < 3 |
                                                                       MB59 < 3 & MC90 >= 3 & MC98 >= 3 & MC73 >= 3 & MF07 >= 3 & MF30 < 3 |
                                                                       
                                                                       MB59 >= 3 & MC90 >= 3 & MC98 >= 3  & MC73 < 3 & MF07 < 3 & MF30 >= 3 |
                                                                       MB59 >= 3 & MC90 >= 3 & MC98 < 3  & MC73 >= 3 & MF07 < 3 & MF30 >= 3 |
                                                                       MB59 >= 3 & MC90 < 3 & MC98 >= 3  & MC73 >= 3 & MF07 < 3 & MF30 >= 3 |
                                                                       MB59 < 3 & MC90 >= 3 & MC98 >= 3  & MC73 >= 3 & MF07 < 3 & MF30 >= 3 |
                                                                       
                                                                       MB59 >= 3 & MC90 >= 3 & MC98 < 3  & MC73 < 3 & MF07 >= 3 & MF30 >= 3 |
                                                                       MB59 >= 3 & MC90 < 3 & MC98 >= 3  & MC73 < 3 & MF07 >= 3 & MF30 >= 3 |
                                                                       MB59 < 3 & MC90 >= 3 & MC98 >= 3  & MC73 < 3 & MF07 >= 3 & MF30 >= 3 |
                                                                       
                                                                       MB59 >= 3 & MC90 < 3 & MC98 < 3  & MC73 >= 3 & MF07 >= 3 & MF30 >= 3 |
                                                                       MB59 < 3 & MC90 >= 3 & MC98 < 3  & MC73 >= 3 & MF07 >= 3 & MF30 >= 3 |
                                                                       
                                                                       MB59 < 3 & MC90 < 3 & MC98 >= 3  & MC73 >= 3 & MF07 >= 3 & MF30 >= 3,
                                                                     "twoAnimals",
                                                                     ifelse(MB59 >= 3 & MC90 >= 3 & MC98 >= 3 & MC73 >= 3 & MF07 >= 3 & MF30 < 3 |
                                                                              MB59 >= 3 & MC90 >= 3 & MC98 >= 3  & MC73 >= 3 & MF07 < 3 & MF30 >= 3 |
                                                                              MB59 >= 3 & MC90 >= 3 & MC98 >= 3  & MC73 < 3 & MF07 >= 3 & MF30 >= 3 |
                                                                              MB59 >= 3 & MC90 >= 3 & MC98 < 3  & MC73 >= 3 & MF07 >= 3 & MF30 >= 3 |
                                                                              MB59 >= 3 & MC90 < 3 & MC98 >= 3  & MC73 >= 3 & MF07 >= 3 & MF30 >= 3 |
                                                                              MB59 < 3 & MC90 >= 3 & MC98 >= 3  & MC73 >= 3 & MF07 >= 3 & MF30 >= 3,
                                                                            "oneAnimal",
                                                                            "none")))))),
                MB59_zero = ifelse(MB59 < 3,"MB59","_"),
                MC90_zero = ifelse(MC90 < 3,"MC90","_"),
                MC98_zero = ifelse(MC98 < 3,"MC98","_"),
                MC73_zero = ifelse(MC73 < 3,"MC73","_"),
                MF07_zero = ifelse(MF07 < 3,"MF07","_"),
                MF30_zero = ifelse(MF30 < 3,"MF30","_"),
                
                which_groups_zero = paste0(MB59_zero,MC90_zero,MC98_zero,MC73_zero,MF07_zero,MF30_zero)) %>%
  rownames_to_column("cluster")

############### conserved across both treatments

for (l in levels(unique(seuratObj_to_use[[subset_res_to_use]])[,subset_res_to_use])){
  
  print(l) ### for debug
  
  #### check have at least one animal
  print(groups_df$num_groups_zero[which(groups_df$cluster == l)])
  if(groups_df$num_groups_zero[which(groups_df$cluster == l)] != "all"){
    


  
  
  curr_df_prep <- FindConservedMarkers(seuratObj_to_use, only.pos = TRUE,
                                  min.pct = minpct_to_use, logfc.threshold = logfc_to_use,
                                  test.use = test_to_use, assay = assay_to_use, recorrect_umi = recorrect_umi_to_use,
                                  ident.1 = l, grouping.var = "Animal") 
  ### only if have anything
  if(dim(curr_df_prep)[1]>0){
    
    curr_df <- curr_df_prep %>%
    rownames_to_column("gene") %>%
    dplyr::mutate(cluster = l,
                  clust_tissue = paste0(l,"_",curr_tissue)) %>%
    ## add groups_df num_groups_zero
    dplyr::left_join(groups_df[,c("cluster","num_groups_zero","which_groups_zero")], by = "cluster") %>%
    ### add back in all the columns in case of missing animals
    dplyr::mutate(MB59_p_val = ifelse(grepl("MB59",which_groups_zero),NA,MB59_p_val),
                  MB59_avg_log2FC = ifelse(grepl("MB59",which_groups_zero),NA,MB59_avg_log2FC),
                  MB59_pct.1 = ifelse(grepl("MB59",which_groups_zero),NA,MB59_pct.1),
                  MB59_pct.2 = ifelse(grepl("MB59",which_groups_zero),NA,MB59_pct.2),
                  MB59_p_val_adj = ifelse(grepl("MB59",which_groups_zero),10,MB59_p_val_adj), ### make a number
                  
                  MC90_p_val = ifelse(grepl("MC90",which_groups_zero),NA,MC90_p_val),
                  MC90_avg_log2FC = ifelse(grepl("MC90",which_groups_zero),NA,MC90_avg_log2FC),
                  MC90_pct.1 = ifelse(grepl("MC90",which_groups_zero),NA,MC90_pct.1),
                  MC90_pct.2 = ifelse(grepl("MC90",which_groups_zero),NA,MC90_pct.2),
                  MC90_p_val_adj = ifelse(grepl("MC90",which_groups_zero),10,MC90_p_val_adj),
                  
                  MC98_p_val = ifelse(grepl("MC98",which_groups_zero),NA,MC98_p_val),
                  MC98_avg_log2FC = ifelse(grepl("MC98",which_groups_zero),NA,MC98_avg_log2FC),
                  MC98_pct.1 = ifelse(grepl("MC98",which_groups_zero),NA,MC98_pct.1),
                  MC98_pct.2 = ifelse(grepl("MC98",which_groups_zero),NA,MC98_pct.2),
                  MC98_p_val_adj = ifelse(grepl("MC98",which_groups_zero),10,MC98_p_val_adj),
                  
                  MC73_p_val = ifelse(grepl("MC73",which_groups_zero),NA,MC73_p_val),
                  MC73_avg_log2FC = ifelse(grepl("MC73",which_groups_zero),NA,MC73_avg_log2FC),
                  MC73_pct.1 = ifelse(grepl("MC73",which_groups_zero),NA,MC73_pct.1),
                  MC73_pct.2 = ifelse(grepl("MC73",which_groups_zero),NA,MC73_pct.2),
                  MC73_p_val_adj = ifelse(grepl("MC73",which_groups_zero),10,MC73_p_val_adj),
                  
                  MF07_p_val = ifelse(grepl("MF07",which_groups_zero),NA,MF07_p_val),
                  MF07_avg_log2FC = ifelse(grepl("MF07",which_groups_zero),NA,MF07_avg_log2FC),
                  MF07_pct.1 = ifelse(grepl("MF07",which_groups_zero),NA,MF07_pct.1),
                  MF07_pct.2 = ifelse(grepl("MF07",which_groups_zero),NA,MF07_pct.2),
                  MF07_p_val_adj = ifelse(grepl("MF07",which_groups_zero),10,MF07_p_val_adj),
                  
                  MF30_p_val = ifelse(grepl("MF30",which_groups_zero),NA,MF30_p_val),
                  MF30_avg_log2FC = ifelse(grepl("MF30",which_groups_zero),NA,MF30_avg_log2FC),
                  MF30_pct.1 = ifelse(grepl("MF30",which_groups_zero),NA,MF30_pct.1),
                  MF30_pct.2 = ifelse(grepl("MF30",which_groups_zero),NA,MF30_pct.2),
                  MF30_p_val_adj = ifelse(grepl("MF30",which_groups_zero),10,MF30_p_val_adj),
                  
                  
                  
                  
                  max_pval = ifelse(num_groups_zero %in% c("none","oneAnimal","twoAnimals","threeAnimals","fourAnimals"),max_pval,NA),
                  minimump_p_val = ifelse(num_groups_zero %in% c("none","oneAnimal","twoAnimals","threeAnimals","fourAnimals"),minimump_p_val,NA),
                  
                  num_groups_zero = NULL,
                  
                  overlap = ifelse(MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05,
                                   "all",
                                   ifelse(MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                            MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                            MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                            MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                            MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                            MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05,
                                          "fiveAnimals",
                                          ifelse(MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                   MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                   
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                   MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                   
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                                   MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                                   
                                                   MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                                   MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                                   
                                                   MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05,
                                                 "fourAnimals",
                                                 ifelse(MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                          
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                          
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                          
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                          
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                          
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                          
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                          
                                                          MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                                          MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05,
                                                        "threeAnimals",
                                                        
                                                        
                                                        
                                                        ifelse(MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj < 0.05 |
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                                 MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                                 
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 
                                                                 MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                                 
                                                                 MB59_p_val_adj < 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05,
                                                               "twoAnimals",
                                                               ifelse(MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05 & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj < 0.05 |
                                                                        MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj < 0.05 & MF30_p_val_adj >= 0.05 |
                                                                        MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj < 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                                        MB59_p_val_adj >= 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj < 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                                        MB59_p_val_adj >= 0.05 & MC90_p_val_adj < 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05 |
                                                                        MB59_p_val_adj < 0.05 & MC90_p_val_adj >= 0.05 & MC98_p_val_adj >= 0.05  & MC73_p_val_adj >= 0.05 & MF07_p_val_adj >= 0.05 & MF30_p_val_adj >= 0.05,
                                                                        "oneAnimal",
                                                                      "none"))))))
    ) %>%
    dplyr::mutate(overlap = factor(overlap, rev(c("none","oneAnimal","twoAnimals","threeAnimals","fourAnimals","fiveAnimals","all")))) %>%
    dplyr::arrange(overlap, max_pval)

  
  
  list_of_consMarkers[[paste0(l,"_",curr_tissue)]] <- curr_df[,c("gene",
                                         "MB59_p_val", "MB59_avg_log2FC", "MB59_pct.1",	"MB59_pct.2",	"MB59_p_val_adj",
                                         "MC90_p_val", "MC90_avg_log2FC", "MC90_pct.1",	"MC90_pct.2",	"MC90_p_val_adj",
                                         "MC98_p_val", "MC98_avg_log2FC", "MC98_pct.1",	"MC98_pct.2",	"MC98_p_val_adj",
                                         "MC73_p_val", "MC73_avg_log2FC", "MC73_pct.1",	"MC73_pct.2",	"MC73_p_val_adj",
                                         "MF07_p_val", "MF07_avg_log2FC", "MF07_pct.1",	"MF07_pct.2",	"MF07_p_val_adj",
                                         "MF30_p_val", "MF30_avg_log2FC", "MF30_pct.1",	"MF30_pct.2",	"MF30_p_val_adj",
                                         "max_pval","minimump_p_val","cluster","clust_tissue","which_groups_zero","overlap")] 
  consMarkers_df <- rbind(consMarkers_df,curr_df[,c("gene",
                                                    "MB59_p_val", "MB59_avg_log2FC", "MB59_pct.1",	"MB59_pct.2",	"MB59_p_val_adj",
                                                    "MC90_p_val", "MC90_avg_log2FC", "MC90_pct.1",	"MC90_pct.2",	"MC90_p_val_adj",
                                                    "MC98_p_val", "MC98_avg_log2FC", "MC98_pct.1",	"MC98_pct.2",	"MC98_p_val_adj",
                                                    "MC73_p_val", "MC73_avg_log2FC", "MC73_pct.1",	"MC73_pct.2",	"MC73_p_val_adj",
                                                    "MF07_p_val", "MF07_avg_log2FC", "MF07_pct.1",	"MF07_pct.2",	"MF07_p_val_adj",
                                                    "MF30_p_val", "MF30_avg_log2FC", "MF30_pct.1",	"MF30_pct.2",	"MF30_p_val_adj",
                                                    "max_pval","minimump_p_val","cluster","clust_tissue","which_groups_zero","overlap")])
  #### if find markers has no lines
  }else{
  curr_df_prep <- NULL
  }
  
  ##### if had no animals with at least 3 cells
  }else{
    curr_df_prep <- NULL
  }
  
}

}

#### save first in case gets stuck
saveRDS(list_of_consMarkers,file.path(outputRData_dir,paste0(params$study_name,
                                                             "_FindConservedMarkers_",
                                                             test_to_use,"_",assay_to_use,
                                                             "_minpct",str_replace(minpct_to_use,"\\.","pt"),
                                                             "_logfc",str_replace(logfc_to_use,"\\.","pt"),
                                                             "_",subset_res_to_use,".rds")))
saveRDS(consMarkers_df,file.path(outputRData_dir,paste0(params$study_name,
                                                        "_FindConservedMarkers_DF_",
                                                        test_to_use,"_",assay_to_use,
                                                        "_minpct",str_replace(minpct_to_use,"\\.","pt"),
                                                        "_logfc",str_replace(logfc_to_use,"\\.","pt"),
                                                        "_",subset_res_to_use,".rds")))


#### add summary tabs


overlaps_have <- sort(unique(consMarkers_df$overlap))
degs_per_cluster_tab_prep <- data.frame("clust_tissue"=factor(unique(consMarkers_df$clust_tissue),as.numeric(unique(consMarkers_df$clust_tissue)))) %>%
  dplyr::arrange(clust_tissue) %>%
  dplyr::mutate(clust_tissue=as.character(clust_tissue))
for(ov in overlaps_have){
  curr_ov_df <- as.data.frame(table(consMarkers_df$clust_tissue[which(consMarkers_df$overlap == ov)]))
  if(dim(curr_ov_df)[1] > 0){
    colnames(curr_ov_df) <- c("clust_tissue",paste0("Num_",ov,"consDEGs"))
  }else{
    curr_ov_df <- data.frame(matrix(nrow = 0, ncol = 2))
    colnames(curr_ov_df) <- c("clust_tissue",paste0("Num_",ov,"consDEGs"))
  }
  degs_per_cluster_tab_prep <- full_join(degs_per_cluster_tab_prep, curr_ov_df, by = "clust_tissue")
}

list_of_consMarkers[["degs_per_cluster_tab"]] <- degs_per_cluster_tab_prep

saveRDS(list_of_consMarkers,file.path(outputRData_dir,paste0(params$study_name,
                                                             "_FindConservedMarkers_",
                                                             test_to_use,"_",assay_to_use,
                                                             "_minpct",str_replace(minpct_to_use,"\\.","pt"),
                                                             "_logfc",str_replace(logfc_to_use,"\\.","pt"),
                                                             "_",subset_res_to_use,".rds")))
saveRDS(consMarkers_df,file.path(outputRData_dir,paste0(params$study_name,
                                                        "_FindConservedMarkers_DF_",
                                                        test_to_use,"_",assay_to_use,
                                                        "_minpct",str_replace(minpct_to_use,"\\.","pt"),
                                                        "_logfc",str_replace(logfc_to_use,"\\.","pt"),
                                                        "_",subset_res_to_use,".rds")))

### save as excel last, so rds can be used to rename tabs if > 31 characters
# can't handle / in names
list_of_consMarkers_for_xlsx <- list_of_consMarkers
names(list_of_consMarkers_for_xlsx) <- gsub("[/]"," ",names(list_of_consMarkers_for_xlsx)) %>%
  str_remove_all(" ") 

openxlsx::write.xlsx(list_of_consMarkers_for_xlsx,file.path(outputTABLES_dir,paste0(params$study_name,
                                                                                    "_FindConservedMarkers_",
                                                                                    test_to_use,"_",assay_to_use,
                                                                                    "_minpct",str_replace(minpct_to_use,"\\.","pt"),
                                                                                    "_logfc",str_replace(logfc_to_use,"\\.","pt"),
                                                                                    "_",subset_res_to_use,".xlsx")))

# must have MAST installed https://github.com/RGLab/MAST, a bit faster than refault wilcoxon above
#https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#identify-differential-expressed-genes-across-conditions-1
#https://github.com/satijalab/seurat/issues/6489 bug
##### install from github
#install_github("RGLab/MAST")


