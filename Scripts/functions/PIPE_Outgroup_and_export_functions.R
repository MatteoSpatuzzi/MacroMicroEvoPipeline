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
# phy = Phylo_genes
# phy_reduced = Phy_genes_reduced_mit
# criteria = con_list$mit
# threshold = c(2,2)
# info_table = classification_sqm
# tab_col = "Family"
# tree_file = Tree_1_sqm
# sister_mat = sister_clades
# 
# phy = Phylo_genes
# phy_reduced = Phy_genes_reduced_nuc
# threshold = c(2,2)
# tip_threshold = 10
# criteria <- con_list$nuc
# info_table = classification_sqm
# tab_col = "Family"
# tree_file = Tree_1_sqm
# sister_mat = sister_clades 

# phy = Phylo_genes
# phy_reduced = Phy_genes_reduced_mit
# criteria = con_list$mit
# threshold = c(2,2)
# tip_threshold = 10
# info_table = classification_sqm
# tab_col = "Family"
# tree_file = Tree_1_amph
# sister_mat = sister_clades

find_Outgroups <-  function(phy, phy_reduced, criteria = FALSE,  threshold = c(3,3), tip_threshold = FALSE, info_table, tab_col, tree_file, sister_mat){
  
  Tree <- tree_file
  classification <- info_table
  Phy_genes <- phy
  Phy_genes_reduced <- phy_reduced
  
  Family_vec <- unique(classification[,tab_col])
  sister_list <- sister_mat
  group = criteria 
  
  if( group == FALSE){
    Phy_genes <- Phy_genes
  }else{
    if(all(group %in% colnames(Phy_genes))){
      
      Phy_genes <- Phy_genes[,group]
      Phy_genes <- Phy_genes[which(apply(Phy_genes, 1, function(x) !all(str_remove_all(x, "-") == ""))),]
      print("removed")
    }else{
      stop("Wrong value for attribute \"criteria\"")
    }
  }
  
  Subtrees <- c()
  Subtrees_note <- matrix(ncol= 3, nrow = dim(sister_list)[1])
  
  Coph_tree_mat <- cophenetic.phylo(Tree)
  Coph_tree_fam <- classification[,tab_col][match(rownames(Coph_tree_mat), rownames(classification))]
  Coph_tree_fam <- Coph_tree_fam[1:dim(Coph_tree_mat)[1]]
  Coph_tree_mat_add <- cbind(Coph_tree_fam, Coph_tree_mat)
  
  for(pair in 1:dim(sister_list)[1]){
    
    families = sister_list[pair,]
    search_complete <- FALSE
    coverage_threshold <- threshold[1]
    amount_threshold <- threshold[2]
    Coph_tree_mat_names <- Coph_tree_mat_add
    
    iter = 1
    
    while(search_complete == FALSE){
     
      up_index <- min(max(which(Coph_tree_mat_names[,1] %in% families)) + 1, dim(Coph_tree_mat_names)[1])
      down_index <- max(min(which(Coph_tree_mat_names[,1] %in% families)) -1 , 1)
     
      Families_candidates <- Coph_tree_mat_names[c(up_index,down_index),1]
      
      if(any(Families_candidates %in% families)){
        Families_candidates <- Families_candidates[-which(Families_candidates %in% families)]
      }
      sum_mat <- matrix(ncol = dim(Phy_genes_reduced[,-c(1,2)])[2], nrow = 0)
      
      for(i in 1:length(Families_candidates)){
        
        names_candidate_1 <- rownames(Coph_tree_mat_names)[which(Coph_tree_mat_names[,1] == Families_candidates[i])]
        phy_candidate_1 <- phy[which(rownames(phy) %in% names_candidate_1), criteria]
        if(class(phy_candidate_1) == "character"){
          
          phy_candidate_1 <- as.data.frame(matrix(phy_candidate_1, nrow = 1, byrow = TRUE))
        }
        
        phy_candidate_1_bool <- apply(phy_candidate_1, c(1,2), function(x) as.numeric(str_remove_all(x, "-") != ""))
        phy_candidate_1_bool_sum <- apply(phy_candidate_1_bool,2,FUN = sum)
        sum_mat <- rbind(sum_mat, phy_candidate_1_bool_sum)
      }
      
       if(any(apply(sum_mat, 1, function(x) sum(sum_mat >= 2) >= 2))){
         
         sum_vec <- apply(sum_mat, 1, function(x) sum(x >= 2))
         candidate_nr <- which(sum_vec == max(sum_vec))
         
         if(length(candidate_nr) == 2){
           
           sum_vec_2 <- c(sum(sum_mat[1,]),sum(sum_mat[1,]))
           candidate_nr <- which(sum_vec_2 == max(sum_vec_2))[1]
         }
         
       candidate = Families_candidates[candidate_nr]
       search_complete <- TRUE
         
       }else{
         
         Coph_tree_mat_names <- Coph_tree_mat_names[-which(Coph_tree_mat_names[,1] %in% Families_candidates),]
         iter = iter +1
         print(pair)
         print(iter)
       }
    } 
    
    subtree_tips <- rownames(Phy_genes_reduced)[which(Phy_genes_reduced$Family %in% families)]
    subtree_tips <- subtree_tips[which(subtree_tips %in% rownames(Phy_genes_reduced))]
    
    outgroup_tips <- rownames(Coph_tree_mat_names)[which(Coph_tree_mat_names[,1] %in% candidate)]
    
    if(class(tip_threshold) == "numeric" & length(outgroup_tips) > tip_threshold){
      
      Phy_Genes_guide <- data.frame(Phy_genes_reduced[which(Phy_genes_reduced[,1] == pair),-c(1,2)])
      Phy_Genes_reduce <- data.frame(Phy_genes[which(rownames(Phy_genes) %in% outgroup_tips),criteria])
      output_ordered_tips <- order_by_gene(Phy_Genes_reduce, Phy_Genes_guide)
      outgroup_tips <- names(output_ordered_tips[1:min(length(output_ordered_tips),tip_threshold)])
    }
    
    subtree_tips <- append(subtree_tips,outgroup_tips)
    subtree_tips <- subtree_tips[which(subtree_tips %in% rownames(Phy_genes))]
    
    Subtrees_note[pair,] <- append(families, candidate)
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

write_subtrees <- function(phy, phy_reduced, subtrees_list, subtree_note, criteria = FALSE, trim = TRUE){
  
  Phy_genes <- phy
  Phy_genes_reduced <- phy_reduced
  wd <- getwd()
  system2("mkdir", args = "subtrees")
  setwd("subtrees")
  group = criteria
  print(length(subtrees_list))
  
  for(subtree in 1: length(subtrees_list)){
    
    if(length(which(Phy_genes_reduced$Pair == subtree)) == 0 ) next
    
    subtree_tips <- subtrees_list[[subtree]]$tip.label
    pair_tips <- rownames(Phy_genes_reduced)[which(Phy_genes_reduced$Family %in% subtree_note[subtree, c(1,2)])]
    outgroup_tips <- subtree_tips[-which(subtree_tips %in% pair_tips)]
    
    Phy_pairs <- Phy_genes_reduced[which(rownames(Phy_genes_reduced) %in% pair_tips),-c(1,2)]
    Phy_outgroup <- Phy_genes[which(rownames(Phy_genes) %in% outgroup_tips),colnames(Phy_genes_reduced)[-c(1,2)]]
    Phy_subtree_seq <- phy_concat(rbind(Phy_pairs, Phy_outgroup), group)
    str_remove_all
    
    if(trim == TRUE & (str_length(Phy_subtree_seq[1,1]) %% 3)!= 0){
      residue <- str_length(Phy_subtree_seq[1,1]) %% 3
      Phy_subtree_seq[,1] <- substr(Phy_subtree_seq[,1],1,str_length(Phy_subtree_seq[1,1])- residue)
    }
    
    dir <- paste("subtree_", subtree_note[subtree,1], "_", subtree_note[subtree,2],"_out_", subtree_note[subtree,3], sep = "")
    system2("mkdir", args = dir)
    setwd(dir)
    
    if(class(criteria) == "list"){
      dir <- paste(dir, "_", names(criteria))
    }
    
    write.csv(Phy_subtree_seq, file = paste(dir, ".csv", sep = ""), row.names = TRUE)
    write.tree(subtrees_list[[subtree]], file = paste(dir, "_tree.tre", sep = ""))
    write.fasta(as.list(Phy_subtree_seq[,1]), names = rownames(Phy_subtree_seq), file.out = paste(dir, ".fasta", sep = ""), as.string =  TRUE )
    
    setwd("..")
  }
  setwd(wd)
}

