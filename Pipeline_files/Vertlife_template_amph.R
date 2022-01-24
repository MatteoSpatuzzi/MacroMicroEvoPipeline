### Data import/creation ####

#Import packages
install.packages("ape", "stringr",  "phylotools")

WorkingDirectory <- "Data/VertLife_amph"
setwd(WorkingDirectory)

library(phylotools)
library(ape)
library(stringr) 
library(seqinr) 
library(readr)

# phylo file
Phy_dat <- read.phylip("amph_shl_new.phy")
# info
classification <- read.csv("classification_tab.csv")
# tree
Tree_1 <- read.tree("amph_shl_new.tre")
dated_Tree <- read.tree("Data/VertLife_amph/amph_shl_new_Posterior_7238.1000-10000.trees/amph_shl_new_Posterior_7238.1000.trees")

# gene mat to split the phylo file correctly
gene_mat <- matrix(c(
  "BDNF", 3282, 4003,
  "CXCR4",  4005, 4765,
  "cytb", 4766, 5905,
  "H3A", 5907, 6232,
  "NCX1", 6234, 7567,
  "ND1", 7568, 8521,
  "ND2", 8522, 9493,
  "POMC", 9495, 10207,
  "RAG1", 10209, 12637,
  "RHOD", 12639, 12952,
  "SIA", 12954, 13348,
  "SLC8A3", 13350, 14488,
  "TYR", 14490, 15091), ncol = 3, byrow = TRUE)
#concat list to concatenate the genes correctly

con_list <- list(
  mit = c("cytb","ND1","ND2"),
  nuc = c("BDNF", "CXCR4","H3A","NCX1","POMC", "RAG1", "RHOD", "RHOD", "SIA", "SLC8A3", "TYR")
)

#### Pipeline ####

# Splits data MA string in gene columns
Phylo_genes <- Edit_phylo(file = Phy_dat, split_mat = gene_mat)

MAFFT_run("phy_genes_fasta", ,FALSE)

# Simply returns information about monophyly, not necessary
#Tree_monophyletic_results <- Tree_analyse(Tree_1, classification)

# Picks sister clades based on family
sister_clades <- Pick_sister_Clades(Tree_1, classification)
# Returns information about sister pairs, necessary for later functions as well
sister_clade_data <- Sister_clades_stats(sister_clades, classification)
# Returns information about clades,can be switched to family ("F") or genus ("G")
Clade_data <- Clade_stats(Phylo_genes, classification, "F")

#Reduces larger clade size based on the availability of mitochondrial genes
Phy_genes_reduced_mit <- Reduce_tips(Phylo_genes, group = "mit", concat_list = con_list, sister_clades, sister_data = sister_clade_data, classification, Tree_1)
#Finds outgroups based on the mitochondrial genes
Subtrees_note_list_mit <- Find_Outgroups(Phy_genes_reduced_mit, Phylo_genes, c(2,2), classification, Tree_1, sister_clades )
WorkingDirectory <- "/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/Pipeline_data"
setwd(WorkingDirectory)
system2("mkdir", args = "subtrees")
setwd("subtrees")
#Exports fasta files for mitochondiral genes, if necessary trims the strings length to be divisible by 3 (x %%3 == 0)
write_subtrees(Phy_genes_reduced_mit, Subtrees_note_list_mit[[1]], Subtrees_note_list_mit[[2]], "mit", con_list[1], trim = TRUE)
setwd(WorkingDirectory)

WorkingDirectory <- "/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/Pipeline_data"
setwd(WorkingDirectory)
#Reduces larger clade size based on the availability of nuclear genes
Phy_genes_reduced_nuc <- Reduce_tips(Phylo_genes, group = "nuc", concat_list = con_list, sister_clades, sister_data = sister_clade_data, classification, Tree_1)
#Finds outgroups based on the nuclear genes
Subtrees_note_list_nuc <- Find_Outgroups(Phy_genes_reduced_nuc, Phylo_genes, c(2,2), classification, Tree_1, sister_clades )
WorkingDirectory <- "/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/Pipeline_data"
setwd(WorkingDirectory)
setwd("subtrees")
#Exports fasta files for nuclear genes, if necessary trims the strings length to be divisible by 3 (x %%3 == 0)
write_subtrees(Phy_genes_reduced_nuc, Subtrees_note_list_nuc[[1]], Subtrees_note_list_nuc[[2]], "nuc",con_list[2], trim = TRUE)
setwd(WorkingDirectory)


### PAMML Results


# Nuclear data

setwd("/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/PAML")
New_branch_tree_files_base <- branch_length_treefiles("base_nuc", dated_Tree, Subtrees_note_list_nuc, Phy_genes_reduced_nuc)
Nuc_base_contrast <- contrast_calc(New_branch_tree_files_base, TRUE)

setwd("/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/PAML")
setwd("code_nuc")
New_branch_tree_files_base<- branch_length_treefiles("dN",dated_Tree, Subtrees_note_list_nuc, Phy_genes_reduced_nuc)
Nuc_dN_contrast <- contrast_calc(New_branch_tree_files_base, TRUE)

setwd("/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/PAML")
setwd("code_nuc")
New_branch_tree_files_base<- branch_length_treefiles("dS",dated_Tree, Subtrees_note_list_nuc, Phy_genes_reduced_nuc)
Nuc_dS_contrast <- contrast_calc(New_branch_tree_files_base, TRUE)


# Mitochondrial data

setwd("/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/PAML")
New_branch_tree_files_base<- branch_length_treefiles("base_mit",dated_Tree, Subtrees_note_list_mit, Phy_genes_reduced_mit)
Mit_base_contrast <- contrast_calc(New_branch_tree_files_base, TRUE)


setwd("/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/PAML")
setwd("code_mit")
New_branch_tree_files_base<- branch_length_treefiles("dN",dated_Tree, Subtrees_note_list_mit, Phy_genes_reduced_mit)
Mit_dN_contrast <- contrast_calc(New_branch_tree_files_base, TRUE)

setwd("/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/PAML")
setwd("code_mit")
New_branch_tree_files_base<- branch_length_treefiles("dS",dated_Tree, Subtrees_note_list_mit, Phy_genes_reduced_mit)
Mit_dS_contrast <- contrast_calc(New_branch_tree_files_base, TRUE)

##### Substitution rate contrast vs. number of species constrast #####

Results <- matrix(c(
  summary(Nuc_base_contrast[[3]])$coefficients[2,],
  summary(Nuc_dN_contrast[[3]])$coefficients[2,],
  summary(Nuc_dS_contrast[[3]])$coefficients[2,],
  summary(Mit_base_contrast[[3]])$coefficients[2,],
  summary(Mit_dN_contrast[[3]])$coefficients[2,],
  summary(Mit_dS_contrast[[3]])$coefficients[2,]),
  ncol = 4, nrow = 6, byrow = TRUE
)

rownames(Results) <- c("Nuc_base", "Nuc_dN", "Nuc_dS", "Mit_base","Mit_dN","Mit_dS")
colnames(Results) <- colnames(summary(Nuc_base_contrast[[3]])$coefficients)

plot(Nuc_base_contrast[[1]], Nuc_base_contrast[[2]])
plot(Nuc_dN_contrast[[1]], Nuc_dN_contrast[[2]])
plot(Nuc_dS_contrast[[1]], Nuc_dS_contrast[[2]])
plot(Mit_base_contrast[[1]], Mit_base_contrast[[2]])
plot(Mit_dN_contrast[[1]], Mit_dN_contrast[[2]])
plot(Mit_dS_contrast[[1]], Mit_dS_contrast[[2]])

##### Substitution rate contrast vs. pair age#####

Results2 <- matrix(c(
  summary(Nuc_base_contrast[[4]])$coefficients[2,],
  summary(Nuc_dN_contrast[[4]])$coefficients[2,],
  summary(Nuc_dS_contrast[[4]])$coefficients[2,],
  summary(Mit_base_contrast[[4]])$coefficients[2,],
  summary(Mit_dN_contrast[[4]])$coefficients[2,],
  summary(Mit_dS_contrast[[4]])$coefficients[2,]),
  ncol = 4, nrow = 6, byrow = TRUE
)

rownames(Results2) <- c("Nuc_base", "Nuc_dN", "Nuc_dS", "Mit_base","Mit_dN","Mit_dS")
colnames(Results2) <- colnames(summary(Nuc_base_contrast[[3]])$coefficients)

##### Number of species contrast vs. pair age#####

Results3 <- matrix(c(
  summary(Nuc_base_contrast[[5]])$coefficients[2,],
  summary(Nuc_dN_contrast[[5]])$coefficients[2,],
  summary(Nuc_dS_contrast[[5]])$coefficients[2,],
  summary(Mit_base_contrast[[5]])$coefficients[2,],
  summary(Mit_dN_contrast[[5]])$coefficients[2,],
  summary(Mit_dS_contrast[[5]])$coefficients[2,]),
  ncol = 4, nrow = 6, byrow = TRUE
)

rownames(Results3) <- c("Nuc_base", "Nuc_dN", "Nuc_dS", "Mit_base","Mit_dN","Mit_dS")
colnames(Results3) <- colnames(summary(Nuc_base_contrast[[3]])$coefficients)




