 ### Data import/creation ####

#Import packages
# install.packages("ape", "stringr",  "phylotools")

 proj_wd <- getwd()
 setwd(proj_wd)

library(phylotools)
library(ape)
library(stringr) 
library(seqinr) 
library(readr)

WorkingDirectory = "Scripts/functions"
setwd(WorkingDirectory)

files <- system2("ls", stdout = TRUE)
lapply(files, FUN = source)

setwd(proj_wd)

WorkingDirectory = "Data/Pipe_amph"
setwd(WorkingDirectory)
 
# phylo file
Phy_dat <- read.phylip("amph_shl_new.phy")
# info
classification_amph <- read.csv("classification_tab.csv", row.names = 1)
classification_amph <- classification_amph[-7239,]
# tree
Tree_1_amph <- read.tree("amph_shl_new_Consensus_7238.tre")
dated_Tree <- read.tree("amph_shl_new_Posterior_7238.1000.trees")

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

Phylo_genes <- edit_phylo(phy = Phy_dat, split_mat = gene_mat, out = "Results", name = "Phy_Genes_1.csv")

dataframe_to_fasta(phy = Phylo_genes, out = "Results", newdir = "Posterior_fasta")

# MAFFT_run(sequ = "Results/Posterior_fasta", new_sequ = "Results/Trial_py" , filewise = TRUE)

# Simply returns information about monophyly, not necessary
Tree_monophyletic_results <- tree_analyse(tree_file = Tree_1_amph, info_table = classification_amph, tab_col = "Family")

# Returns information about clades,can be switched to various columns
Clade_data <- clade_stats(Phylo_genes, classification_amph, "Family")

# Picks sister clades based on family
sister_clades <- pick_sister_Clades(Tree_1_amph, classification_amph, "Family")
# Returns information about sister pairs, necessary for later functions as well
sister_clade_data <- sister_clades_stats(sister_clades, classification_amph, "Family")

#Reduces larger clade size based on the availability of mitochondrial genes
Phy_genes_reduced_mit <- reduce_tips(phy = Phylo_genes, criteria = con_list$mit, sister_mat = sister_clades, sister_data = sister_clade_data, classification_amph, "Family", Tree_1_sqm)
#Finds outgroups based on the mitochondrial genes
Subtrees_note_list_mit <- find_Outgroups(phy = Phylo_genes, phy_reduced = Phy_genes_reduced_mit, threshold = c(2,2), info_table = classification_amph, tab_col = "Family", tree_file = Tree_1_amph, sister_mat = sister_clades )
#Exports fasta files for mitochondiral genes, if necessary trims the strings length to be divisible by 3 (x %%3 == 0)
setwd("Results")
write_subtrees(Phy_genes_reduced_mit, Subtrees_note_list_mit[[1]], Subtrees_note_list_mit[[2]], con_list[1], trim = TRUE)
setwd("..")


#Reduces larger clade size based on the availability of nuclear genes
Phy_genes_reduced_nuc <- reduce_tips(phy = Phylo_genes, criteria = con_list$nuc, sister_mat = sister_clades, sister_data = sister_clade_data, classification_amph, "Family", Tree_1_sqm)
#Finds outgroups based on the nuclear genes
Subtrees_note_list_nuc <- find_Outgroups(phy = Phylo_genes, phy_reduced = Phy_genes_reduced_nuc, threshold = c(2,2), info_table = classification_amph, tab_col = "Family", tree_file = Tree_1_amph, sister_mat = sister_clades )
#Exports fasta files for nuclear genes, if necessary trims the strings length to be divisible by 3 (x %%3 == 0)
setwd("Results")
write_subtrees(Phy_genes_reduced_nuc, Subtrees_note_list_nuc[[1]], Subtrees_note_list_nuc[[2]], con_list[2], trim = TRUE)
setwd("..")

### PAMML Results

# Nuclear data

setwd("PAML_amph")
New_branch_tree_files_base <- branch_length_treefiles("base_nuc", dated_Tree, info_tab = classification_amph, Subtrees_note_list_nuc, Phy_genes_reduced_nuc)
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




