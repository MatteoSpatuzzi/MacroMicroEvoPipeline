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

WorkingDirectory = "Data/Pipe_sqm"
setwd(WorkingDirectory)

# phylo file
Phy_dat_sqm <- read.phylip("Vertlife/squam_shl_new.phy")

# info
classification_sqm <- read.csv("classification_tab.csv", row.names = 1)

# trees
Tree_1_sqm <- read.tree("Vertlife/squam_shl_new_Consensus_9755.tre")
dated_Tree_sqm <- read.tree("Vertlife/squam_shl_new_Posterior_9755.1000.trees")

# gene mat to split the phylo file correctly
GenBank_matrix_sqm <- read.csv("GenBank_matrix_sqm.csv", row.names = 1)

gene_mat <- matrix(c(
  "12S",
  "16S",
  "AMEL",
  "BDNF",
  "BMP2", 
  "CMOS",
  "COI",
  "CYTB",
  "ND1",
  "ND2",
  "ND4",
  "NT3",
  "PDC",
  "PRLR",
  "R35",
  "RAG1",
  "RAG2"), ncol = 1)

gene_mat <- cbind(gene_mat, find_gene_loci(phy = Phy_dat_sqm, genbank_dat = GenBank_matrix_sqm))
#concat list to concatenate the genes correctly

con_list <- list(
  mit = c("CYTB",
          "ND1",
          "ND2",
          "ND4"),
  nuc = c("AMEL",
          "BDNF",
          "BMP2", 
          "CMOS","NT3",
          "PDC",
          "PRLR",
          "R35",
          "RAG1",
          "RAG2")
)

#### Pipeline ####

# Splits data MA string in gene columns
Phylo_genes <- edit_phylo(phy = Phy_dat_sqm, split_mat = gene_mat, out = "Results", name = "Phy_Genes_1.csv")

dataframe_to_fasta(phy = Phylo_genes, out = "Results", newdir = "Posterior_fasta")

# MAFFT_run(sequ = "Results/Posterior_fasta", new_sequ = "Results/Trial_py" , filewise = TRUE)

# Simply returns information about monophyly, not necessary
Tree_monophyletic_results <- tree_analyse(tree_file = Tree_1_sqm, info_table = classification_sqm, tab_col = "Family")

# Returns information about clades,can be switched to various columns
Clade_data <- clade_stats(Phylo_genes, classification_sqm, "Family")

# Picks sister clades based on family
sister_clades <- pick_sister_Clades(Tree_1_sqm, classification_sqm, "Family")
# Returns information about sister pairs, necessary for later functions as well
sister_clade_data <- sister_clades_stats(Phylo_genes, sister_clades, classification_sqm, "Family")

#Reduces larger clade size based on the availability of mitochondrial genes
Phy_genes_reduced_mit <- reduce_tips(phy = Phylo_genes, criteria = con_list$mit, tip_threshold = 10, sister_mat = sister_clades, sister_data = sister_clade_data, classification_sqm, "Family", Tree_1_sqm)
#Finds outgroups based on the mitochondrial genes
Subtrees_note_list_mit <- find_Outgroups(phy = Phylo_genes, phy_reduced = Phy_genes_reduced_mit, criteria = con_list$mit, threshold = c(2,2), tip_threshold = 10, info_table = classification_sqm, tab_col = "Family", tree_file = Tree_1_sqm, sister_mat = sister_clades )

#Exports fasta files for mitochondiral genes, if necessary trims the strings length to be divisible by 3 (x %%3 == 0)
setwd("Results")
write_subtrees(Phylo_genes, Phy_genes_reduced_mit, Subtrees_note_list_mit[[1]], Subtrees_note_list_mit[[2]], con_list[1], trim = TRUE)
setwd("..")

#Reduces larger clade size based on the availability of nuclear genes
Phy_genes_reduced_nuc <- reduce_tips(phy = Phylo_genes, criteria = con_list$nuc, tip_threshold = 10, sister_mat = sister_clades, sister_data = sister_clade_data, classification_sqm, "Family", Tree_1_sqm)
#Finds outgroups based on the nuclear genes
Subtrees_note_list_nuc <- find_Outgroups(phy = Phylo_genes, phy_reduced = Phy_genes_reduced_nuc, criteria = con_list$nuc, threshold = c(2,2), tip_threshold = 10, info_table = classification_sqm, tab_col = "Family", tree_file = Tree_1_sqm, sister_mat = sister_clades )

#Exports fasta files for nuclear genes, if necessary trims the strings length to be divisible by 3 (x %%3 == 0)
setwd("Results")
write_subtrees(Phylo_genes, Phy_genes_reduced_nuc, Subtrees_note_list_nuc[[1]], Subtrees_note_list_nuc[[2]], con_list[2], trim = TRUE)
setwd("..")

### PAMML Results

# Nuclear data

setwd("Results")
PAML_res_dir <- "PAML/baseml"
setwd(PAML_res_dir)

New_branch_tree_files_base <- branch_length_treefiles("nuc", dated_Tree_sqm, info_tab = classification_sqm, Subtrees_note_list_nuc, Phy_genes_reduced_nuc)
Nuc_base_contrast <- contrast_calc(New_branch_tree_files_base, TRUE)
setwd("..")
setwd("..")
evaluate_linear_models("squamata_baseml_nuclear", Nuc_base_contrast, New_branch_tree_files_base)
setwd("..")

setwd("baseml")
New_branch_tree_files_base <- branch_length_treefiles("mit", dated_Tree_sqm, info_tab = classification_sqm , Subtrees_note_list_mit, Phy_genes_reduced_mit)
Mit_base_contrast <- contrast_calc(New_branch_tree_files_base, TRUE)
setwd("..")
setwd("..")
evaluate_linear_models("squamata_baseml_mitochondrial", Mit_base_contrast, New_branch_tree_files_base)




