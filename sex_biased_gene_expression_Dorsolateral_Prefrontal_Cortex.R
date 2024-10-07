library(Seurat)
library(stringr)
library(presto)                         
library(clusterProfiler)                #for GO enrichment analysis
library(org.Hs.eg.db)                   #for human gene annotation database

#files used here were downloaded from the Human Cell Atlas Data Portal
#Based on the availability of adequate male + female scRNA data, I chose data for 
#the Dorsolateral Prefrontal Cortex, or DLPFC [Jorstad et al., 2023]
#LINK: https://data.humancellatlas.org/hca-bio-networks/nervous-system/atlases/cortex-v1-0

#This is an exploratory analysis, aimed to find if there's significant sex differences in gene expression in the DLPFC
#I perform differential expression over the 3 most frequent cell types, which incidentally are the same for males and females


###################################################################################
#load files, subset data and do some pre-processing before differential expression#
###################################################################################
#load rds file for dorsolateral prefrontal cortex (DFC)
dfc <- readRDS("C:/Users/DELL/Desktop/scRNA/Data/dorsolateral_prefrontal_cortex.rds")
dfc_female <- subset(dfc, subset = development_stage == "43-year-old human stage")  #5591 female samples
dfc_male <- subset(dfc, subset = development_stage == "42-year-old human stage")    #30804 male samples

#remove original .rds file to save memory
rm(dfc)

#subset the male dataframe to have same number of cells as female df
#randomly select 5591 cells from the male dataset
sampled_male_cellnames <- sample(Cells(dfc_male), size = 5591, replace = F)
dfc_male <- subset(dfc_male, cells = sampled_male_cellnames)

#check dimensions to make sure
print(dim(dfc_female))
print(dim(dfc_male))

#identify the 3 most frequent cell types in male and female dfs
f_cell_type_counts <- table(dfc_female$cell_type)
m_cell_type_counts <- table(dfc_male$cell_type)

f_common_celltypes <- names(sort(f_cell_type_counts, decreasing = T)[1:3])
m_common_celltypes <- names(sort(m_cell_type_counts, decreasing = T)[1:3])
print(f_common_celltypes)    #the top 3 frequent cell types are the same in males and females
print(m_common_celltypes)

#subset the male/female dfs to only include these common cell types
dfc_female <- subset(dfc_female, subset = cell_type %in% f_common_celltypes)
dfc_male <- subset(dfc_male, subset = cell_type %in% m_common_celltypes)

#normalise and scale the data 
dfc_female <- NormalizeData(dfc_female)
dfc_male <- NormalizeData(dfc_male)

dfc_female <- ScaleData(dfc_female)
dfc_male <- ScaleData(dfc_male)

#merge male and female dataframes; set 'sex' as an active identity class
#free up memory - delete earlier datasets
merged_dfc <- merge(dfc_male, y = dfc_female, add.cell.ids = c("male", "female"), project = "combined_sex")
Idents(merged_dfc) <- "sex"
rm(dfc_female, dfc_male)

#NOTE: alternately, merge data, then normalise + scale the joint male/female data (shouldn't make a difference, however)

##############################################################
#differential expression for each cell type
data_celltype1 <- subset(merged_dfc, subset = cell_type == f_common_celltypes[1])
DE_celltype1 <- FindMarkers(data_celltype1, ident.1 = "female", ident.2 = "male")

data_celltype2 <- subset(merged_dfc, subset = cell_type == f_common_celltypes[2])
DE_celltype2 <- FindMarkers(data_celltype2, ident.1 = "female", ident.2 = "male")

data_celltype3 <- subset(merged_dfc, subset = cell_type == f_common_celltypes[3])
DE_celltype3 <- FindMarkers(data_celltype3, ident.1 = "female", ident.2 = "male")

#genes with significant differential expression (p-adjusted < 0.05)
sig_genes_celltype1 <- rownames(DE_celltype1[which(DE_celltype1$p_val_adj < 0.05),])
sig_genes_celltype2 <- rownames(DE_celltype2[which(DE_celltype2$p_val_adj < 0.05),])
sig_genes_celltype3 <- rownames(DE_celltype3[which(DE_celltype3$p_val_adj < 0.05),])


##############################################################
#GO enrichment analysis (make sure the human genome annotation is installed, or install through bioconductor)
enrich_sig_genes_celltype1 <- enrichGO(gene = sig_genes_celltype1, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
enrich_sig_genes_celltype2 <- enrichGO(gene = sig_genes_celltype2, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
enrich_sig_genes_celltype3 <- enrichGO(gene = sig_genes_celltype3, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")

#visualise enriched GO categories
dotplot(enrich_sig_genes_celltype1)
dotplot(enrich_sig_genes_celltype2)
dotplot(enrich_sig_genes_celltype3)
