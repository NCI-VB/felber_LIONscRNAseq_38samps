################# generates cellchat objects ###################

######## get args from shell script
args <- commandArgs(trailingOnly = TRUE)

study_name_to_use <- as.character(args[1]) ##batch123_38samps_rpca # 
tissue_to_use <- as.character(args[2]) ##LN or PBMC # 
res_to_use <- as.character(args[3]) ## respt4_rpca_ClustersWithSubs
contrast_to_use <- as.character(args[4]) ## DAY2vsPRE, DAY8vsPRE, DAY8vsDAY2
species <- as.character(args[5]) ## human (or could be mouse)
prep_both_or_second_only <- as.character(args[6]) ## "prep" = for generating the individual cellchat only,"both" = for running it all (takes forever with large dateste), or "second_only"=for after running merge of pairs on my computer

############################################################
# packages & functions
############################################################
# package list
list.of.packages=c(
  "tidyverse","strex","Seurat"#"reticulate"
  )


#install as needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
### needs to be installed directly from github
devtools::install_github("jinworks/CellChat")
devtools::install_github('immunogenomics/presto')

# load packages
invisible(lapply(list.of.packages, library, character.only = TRUE))
library(CellChat)
library(presto)
reticulate::py_install(packages = 'umap-learn')
#library(umap-learn)

##### to be run on biowulf, called by ccrvb47_batch123_CellChatPrep_KCG.sh
#study_name_to_use <- "LN"
params <- list(
  "study_name" = study_name_to_use,
  "curr_tissue" = tissue_to_use
)

#### for biowulf
input_dir <- "/data/NCI_VB/felber/ccrvb47_LIONscRNAseq_KCG/RData" ### the same as output for running on biowulf, be sure the ann R objects are up to date here
outputRData_dir <- "/data/NCI_VB/felber/ccrvb47_LIONscRNAseq_KCG/RData"
code_path <- "/data/NCI_VB/goldfarbmurenkc/github/VaccineBranch/Felber/ccrvb47_LIONscRNAseq_batch123_38samps_RPCA_Subclustering_KCG.Rmd"
lib_loc_to_use <- "/data/goldfarbmurenkc/R/rhel8/4.4/"




#### load full dataset
seuratObj_to_use <- readRDS(file.path(input_dir,paste0(params$study_name,"_ann.rds"))) 
### print to see original size
seuratObj_to_use

### set timepoints to use tp2 vs tp1
tp1 <- str_after_first(contrast_to_use,"vs") #"PRE"
tp2 <- str_before_first(contrast_to_use,"vs") #"POST"

#### trim to tissue & time
seuratObj_to_use <- subset(seuratObj_to_use, subset = Tissue == tissue_to_use)
seuratObj_to_use <- subset(seuratObj_to_use, subset = Time %in% c(tp1,tp2))
### also relevel timepoints to only be these 2
seuratObj_to_use$Time <- factor(as.character(seuratObj_to_use$Time), c(tp1, tp2))

### print to see that size changed
seuratObj_to_use

######### remove levels not present in both pre and post
clusters_by_group_df <- as.data.frame.matrix(table(seuratObj_to_use[[c(res_to_use,"Time")]])) %>%
  ### rename to generic
  dplyr::rename(g1=1,g2=2) %>%
  ### add column with if pre/post or both
  dplyr::mutate(helper = ifelse(g1 == 0 | g2 == 0, "not_both","both")) %>%
  ## remove not_both
  dplyr::filter(helper == "both")


################## for both or prep, generate separate cellchat objects ##################
if(prep_both_or_second_only == "both" | prep_both_or_second_only == "prep" ){

#### first separate the POST and PRE
PRE_so <- subset(x = seuratObj_to_use, subset = Time == tp1)
POST_so <- subset(x = seuratObj_to_use, subset = Time == tp2)
#### next make each into cellchat object https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html

###### PRE
### set idents
Idents(PRE_so) <- res_to_use
###### for different composition of celltypes (i.e. some cell pops that are not present in both timepoints), relevel to ony the cell pops present in this dataset
PRE_so[[res_to_use]][,res_to_use] <- factor(PRE_so[[res_to_use]][,res_to_use],
                                                      intersect(levels(PRE_so[[res_to_use]][,res_to_use]),
                                                                as.character(unique(PRE_so[[res_to_use]][,res_to_use]))))

data.input_PRE <- PRE_so[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels_PRE <- Idents(PRE_so)
meta_PRE <- data.frame(labels = labels_PRE, row.names = names(labels_PRE)) # create a dataframe of the cell labels
cellchat_PRE <- createCellChat(object = data.input_PRE, meta = meta_PRE, group.by = "labels")
### warning no samples column, not in tutorial

#####POST
Idents(POST_so) <- res_to_use
###### for different composition of celltypes (i.e. some cell pops that are not present in both dataset), relevel to ony the cell pops present in this dataset
POST_so[[res_to_use]][,res_to_use] <- factor(POST_so[[res_to_use]][,res_to_use],
                                            intersect(levels(POST_so[[res_to_use]][,res_to_use]),
                                                      as.character(unique(POST_so[[res_to_use]][,res_to_use]))))


data.input_POST <- POST_so[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels_POST <- Idents(POST_so)
meta_POST <- data.frame(labels = labels_POST, row.names = names(labels_POST)) # create a dataframe of the cell labels
cellchat_POST <- createCellChat(object = data.input_POST, meta = meta_POST, group.by = "labels")



### save levels for pre and post and all
saveRDS(levels(PRE_so[[res_to_use]][,res_to_use]), file.path(outputRData_dir,paste0(params$study_name,"_",params$curr_tissue,"_",res_to_use,"_diffComp_",contrast_to_use,"_",tp1,"_LEVELS.rds")))
saveRDS(levels(POST_so[[res_to_use]][,res_to_use]), file.path(outputRData_dir,paste0(params$study_name,"_",params$curr_tissue,"_",res_to_use,"_diffComp_",contrast_to_use,"_",tp2,"_LEVELS.rds")))
saveRDS(levels(seuratObj_to_use[[res_to_use]][,res_to_use]), file.path(outputRData_dir,paste0(params$study_name,"_",params$curr_tissue,"_",res_to_use,"_diffComp_",contrast_to_use,"_LEVELS.rds")))



##### set DB
if(species == "human"){
  CellChatDB <- CellChatDB.human
}
if(species == "mouse"){
  CellChatDB <- CellChatDB.mouse
}

#showDatabaseCategory(CellChatDB)
### can subset to just particular signaling, but use all for now
# set the used database in the object
cellchat_PRE@DB <- CellChatDB
cellchat_POST@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
# This step is necessary even if using the whole database
cellchat_PRE <- subsetData(cellchat_PRE) 
cellchat_POST <- subsetData(cellchat_POST)


options(future.globals.maxSize = 400000 * 1024^2) ### set to 200G
future::plan("multisession", workers = 12) # do parallel
### using wilcoxon for overexpressed genes, load presto
cellchat_PRE <- identifyOverExpressedGenes(cellchat_PRE) 
cellchat_POST <- identifyOverExpressedGenes(cellchat_POST) 
cellchat_PRE <- identifyOverExpressedInteractions(cellchat_PRE)
#The number of highly variable ligand-receptor pairs used for signaling inference is 2209 
cellchat_POST <- identifyOverExpressedInteractions(cellchat_POST)
#The number of highly variable ligand-receptor pairs used for signaling inference is 2202


##### also save the cellchat objects
saveRDS(cellchat_PRE, file.path(outputRData_dir,paste0(params$study_name,"_",params$curr_tissue,"_",res_to_use,"_diffComp_",contrast_to_use,"_",tp1,"_cellchat.rds")))
saveRDS(cellchat_POST, file.path(outputRData_dir,paste0(params$study_name,"_",params$curr_tissue,"_",res_to_use,"_diffComp_",contrast_to_use,"_",tp2,"_cellchat.rds")))
}

########################### only for both, this step is usually prohibitively long on biowulf ###################


#### only if running both from slurm
if(prep_both_or_second_only == "both"){
# Compute the communication probability and infer cellular communication network
##### this takes a while (~15 min) for each
ptm = Sys.time()
cellchat_PRE <- computeCommunProb(cellchat_PRE, type = "triMean")
cellchat_POST <- computeCommunProb(cellchat_POST, type = "triMean")
#### not explicitly filtering communication, default is 10
#cellchat <- filterCommunication(cellchat, min.cells = 10)
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
# df.net <- subsetCommunication(cellchat) returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# 
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

# Infer the cell-cell communication at a signaling pathway level
cellchat_PRE <- computeCommunProbPathway(cellchat_PRE)
cellchat_POST <- computeCommunProbPathway(cellchat_POST)

#Calculate the aggregated cell-cell communication network
cellchat_PRE <- aggregateNet(cellchat_PRE)
cellchat_POST <- aggregateNet(cellchat_POST)

# Compute the network centrality scores
cellchat_PRE <- netAnalysis_computeCentrality(cellchat_PRE, slot.name = "netP")
cellchat_POST <- netAnalysis_computeCentrality(cellchat_POST, slot.name = "netP")

##### now merge the 2 cellchat objects https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
### but first "lift up" to have the same cell populations
group.new = union(levels(cellchat_PRE@idents), levels(cellchat_POST@idents))
cellchat_PRE <- liftCellChat(cellchat_PRE, group.new)
cellchat_POST <- liftCellChat(cellchat_POST, group.new)


object.list <- list(PRE = cellchat_PRE, POST = cellchat_POST)
names(object.list) <- c(tp1,tp2) ### rename to reflect contrast
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#### save these initially
save(object.list, file = file.path(input_dir,paste0(params$study_name,"_",params$curr_tissue,"_diffComp_cellchat_object.list_",contrast_to_use,".RData")))
save(cellchat, file = file.path(input_dir,paste0(params$study_name,"_",params$curr_tissue,"_diffComp_cellchat_merged_",contrast_to_use,".RData")))


}




################### for both or second only, complete the cellchat on the merged objects

if(prep_both_or_second_only != "prep"){

  ### if prep or second only
  ### load merged cellchat objects saved locally since runs way faster
  load(file = file.path(input_dir,paste0(params$study_name,"_",params$curr_tissue,"_diffComp_cellchat_object.list_",contrast_to_use,".RData")))
  load(file = file.path(input_dir,paste0(params$study_name,"_",params$curr_tissue,"_diffComp_cellchat_merged_",contrast_to_use,".RData")))
  



########
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
#### for the next step - had to run reticulate::py_install(packages = 'umap-learn') - didn't work
### pip install umap-learn in terminal, said successful, needed to restart rstudio
cellchat <- netEmbedding(cellchat, type = "functional") ### this
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")


##### identify dysfunctional signaling by OE analysis (uses wilcoxon)
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset <- tp2
# define a char name used for storing the results of differential expression analysis
features.name <- paste0(pos.dataset, ".merged")
# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

########## this can only be done for clusters with both groups

cellchat <- identifyOverExpressedGenes(object = cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       idents.use = rownames(clusters_by_group_df)[which(clusters_by_group_df$helper == "both")], #### only perform for clusters with both groups
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.05,
                                       thresh.p = 0.05, 
                                       group.DE.combined = FALSE) 

######### for running source code
# data.use = NULL
# group.by = NULL
# idents.use = NULL
# invert = FALSE
# 
# features = NULL
# return.object = TRUE
# do.DE = TRUE
# do.fast = TRUE 
# min.cells = 10
# 
# object = cellchat 
# group.dataset = "datasets" 
# pos.dataset = pos.dataset
# features.name = features.name
# only.pos = FALSE
# thresh.pc = 0.1 
# thresh.fc = 0.05
# thresh.p = 0.05
# group.DE.combined = FALSE
# 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)


#### save
save(object.list, file = file.path(input_dir,paste0(params$study_name,"_",params$curr_tissue,"_diffComp_cellchat_object.list_",contrast_to_use,".RData"))) 
save(cellchat, file = file.path(input_dir,paste0(params$study_name,"_",params$curr_tissue,"_diffComp_cellchat_merged_",contrast_to_use,".RData")))
saveRDS(net, file = file.path(input_dir,paste0(params$study_name,"_",params$curr_tissue,"_diffComp_cellchat_merged",contrast_to_use,"_updownLRpairs.rds")))


}
