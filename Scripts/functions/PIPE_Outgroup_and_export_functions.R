#### Outgroups ####

# This function finds an outgroup for each sister clade based on the phylogenetic tree and the 
# sequences provided (Outcome may change between the "mit" and "nuc" options from the function beforehand)
# 
# file is the output of the Reduce_tips() function
# info_tab is the supplementary information table with the phylogenetic classification of the tree
# tree_file is the phylognetic tree object
# sister_list is the list of sister pairs in form of a nx2 matrix
#
# It returns a list with two elements:
# a list of subtrees with all the tips from the sister clades and the outgroup
# an expanded verison of the sister_list object with a third column for the outgroup

#'\code{find_Outgroups}
#'
#'@details This function finds an outgroup for each sister clade based on the phylogenetic tree and the 
#' sequences provided 
#'
#'@param phy An object of the phylo class or a dataframe where columns 
#'are the genes and the rows are species. In the latter case it's 
#'important that the column names correspond to the gene names if 
#'one wishes to export the sequences as fasta files for each gene.
#'@param phy_reduced is the same dataframe after being filtered so that 
#'there are an equal amount of species per sister clade in each pair
#'@param threshold a vector of two indicating how many genes must be 
#'contained in the outgroup for it to be chosen, the first number is 
#'the threshold for genes per tip and the second for how many tips 
#'must fulfill the first criteria
#'@param @param info_table a table with as many rows as the tips in the tree file
#'and corresponding names and at least one column to group them
#'@param tab_col The column of the info_table used to group the tree tips
#'@param tree a tree object
#'@param sister_mat a nx2 matrix that lists the sister pairs
#'
#'@return a nx3 matrix analogous to the sister_mat object with 
#'an additional column indicating the outgroup
#'
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import stringr
#'@import seqinr
#'
#'

find_Outgroups <-  function(phy, phy_reduced, threshold = c(3,3), info_table, tab_col, tree_file, sister_mat){
  
  Tree <- tree_file
  classification <- info_table
  # Phy_genes_reduced_mit <- file
  # Phy_genes_reduced <- Phy_genes_reduced_mit
  Phy_genes <- phy
  Phy_genes_reduced <- phy_reduced
  
  Family_vec <- unique(classification[,tab_col])
  sister_list <- sister_mat
  # Phy_names_inTree <- rownames(Phy_genes_reduced)[which(rownames(Phy_genes_reduced) %in% Tree$tip.label)]
  # Tree_phy <- keep.tip(Tree,Phy_names_inTree)
  tips_list_fam <- c()
  
  for (fam in Family_vec) {
    
    Family_names_vec <- rownames(classification)[which(classification[,tab_col] == fam)]
    Family_names_vec <- Family_names_vec[which(Family_names_vec %in% Tree$tip.label)]
    tips_list_fam <- append(tips_list_fam, Family_names_vec[1])
  }
  Tree_1_per_Family <- keep.tip(Tree, tips_list_fam)
  Tree_1_per_Family$tip.label <- Family_vec
  
  
  ##### Find outgroup #####
  
  Subtrees <- c()
  Subtrees_note <- matrix(ncol= 3, nrow = dim(sister_list)[1])
  
  for(pair in 1:dim(sister_list)[1]){
    
    families = sister_list[pair,]
    search_complete <- FALSE
    coverage_threshold <- threshold[1]
    amount_threshold <- threshold[2]
    iter = 1
    
    while(search_complete == FALSE){
      
      outgroup_index <- which(Tree_1_per_Family$tip.label %in% families) + c(-1*iter, iter)
      outgroups <- Tree_1_per_Family$tip.label[outgroup_index[which(outgroup_index > 0)]]
      if(any(is.na(outgroups))){
        outgroups <- outgroups[-which(is.na(outgroups) ==  TRUE)]
      }
      cophenetic_outgroups <- cophenetic.phylo(keep.tip(Tree_1_per_Family, tip = append(families, outgroups)))
      dist_outgroups <- which(cophenetic_outgroups[families[1],] == min(cophenetic_outgroups[outgroups,families[1]]))
      candidate <- names(dist_outgroups)
      
      Phy_Genes_reduce <- data.frame(Phy_genes[rownames(classification)[which(classification[,tab_col] == candidate & rownames(classification) %in% rownames(Phy_genes))],colnames(Phy_genes_reduced)[colnames(Phy_genes_reduced)%in%colnames(Phy_genes)]])
      Phy_Genes_reduce[] <- lapply(Phy_Genes_reduce, function(x) ceiling(str_length(str_remove_all(x,"-"))/(str_length(str_remove_all(x,"-"))+1)))
      Phy_Genes_sum_col <- apply(Phy_Genes_reduce, 2, function(x) sum(x, na.rm = TRUE))
      
      if(sum(Phy_Genes_sum_col > 2) > 2){
        outgroup <- candidate
        search_complete <- TRUE
      }else{
        iter = iter +1
      }
    }

    subtree_tips <- rownames(classification)[which(classification[,tab_col] %in% families)]
    subtree_tips <- subtree_tips[subtree_tips %in% rownames(Phy_genes_reduced)]
    subtree_tips <- append(subtree_tips, rownames(classification)[which(classification[,tab_col] %in% outgroup)])
    subtree_tips <- subtree_tips[subtree_tips %in% rownames(Phy_genes)]
    
    Subtrees_note[pair,] <- append(families, outgroup)
    Subtrees[[pair]] <- keep.tip(Tree, subtree_tips)

  }
  output_dat <- list(Subtrees, Subtrees_note)
  return(output_dat)
}


#'\code{write_subtrees}
#'
#'@details This function writes tree files that represent the subtrees 
#'that each include one sister pqir and the corresponding outgroup
#'
#'@param phy_reduced the dataframe with the multiple alignemnt 
#'sequences split by genes after the removal of species to achieve 
#'equality in the number of tips per sister clade
#'@param subtrees_list the lsit of subtrees generated by the 
#'find_Outgroups() function
#'@param subtrees_note a matrix of information about the sister 
#'pairs and outgroups generated by the find_Outgroups() function
#'@param concat_list a list of vectors, each vector contains the names of 
#'sequences to concatenate in the phy object under the name of the 
#'corresponding list name, if a gene is not present in any vector 
#'it will NOT be present in the result
#'
#'@return exports fasta files in the wd 
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import stringr
#'@import seqinr# This function exports the sequences of a subtree containing two sister clades and an outgroup
# 
# file is the output 
# subtrees_list is the first element of the Find_Outgroups() function output
# subtrees_list is the second element of the Find_Outgroups() function output
# group indicates wether "mit" or "nuc" sequences are being exported
# Concat_list is a list with only one element, a vector with the names of the genes in question
# (is convoluted, will be improved)
#
# It Exports fasta files in the wd 

write_subtrees <- function(phy, subtrees_list, subtree_note, concat_list, trim = TRUE){
  
  Phy_genes_reduced <- phy
  wd <- getwd()
  
  system2("mkdir", args = "subtrees")
  setwd("subtrees")
  
  for(subtree in 1: length(subtrees_list)){
    
    if(length(which(Phy_genes_reduced$Pair ==subtree)) == 0 ) next
    
    subtree_tips <- subtrees_list[[subtree]]$tip.label
    Phy_subtree <- Phy_genes_reduced[which(rownames(Phy_genes_reduced) %in% subtree_tips) ,]
    Phy_subtree_seq <- phy_concat(Phy_subtree, concat_list)

    rownames(Phy_subtree_seq) <- rownames(Phy_subtree)
    
    if(trim == TRUE & (str_length(Phy_subtree_seq[1,1]) %% 3)!= 0){
      
      residue <- str_length(Phy_subtree_seq[1,1]) %% 3
      Phy_subtree_seq[,1] <- substr(Phy_subtree_seq[,1],1,str_length(Phy_subtree_seq[1,1])- residue)
      
    }
    
    dir <- paste("subtree_", subtree_note[subtree,1], "_", subtree_note[subtree,2],"_out_", subtree_note[subtree,3], sep = "")
    system2("mkdir", args = dir)
    setwd(dir)
    
    write.csv(Phy_subtree, file = dir, row.names = TRUE)
    write.tree(subtrees_list[[subtree]], file = paste(dir, "_tree.tre", sep = ""))
    
    # Nuclear fasta
    string <- ""
    for(org in 1:dim(Phy_subtree_seq)[1]){
      string <- paste(string, ">",rownames(Phy_subtree_seq)[org], "\n",Phy_subtree_seq[org,names(concat_list)], "\n", sep = "" )
    }
    
    filename <- paste(dir,"_", names(concat_list),".fasta",sep = "" )
    cat(string, file =  filename)
    setwd("..")
  }
  
}
