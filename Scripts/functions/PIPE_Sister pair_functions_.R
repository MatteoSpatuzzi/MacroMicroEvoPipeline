#### Tree File #### 

#'\code{tree_analyse}
#'
#'@details Analyses monophyly of clades on a tree object according to a given classification
#'
#'@param tree_file a tree object 
#'@param info_table a table with as many rows as the tips in the tree file
#'and corresponding names and at least one column to group them
#'@param tab_col The column of the info_table used to group the tree tips, 
#'species that belong to the same group should have the same name
#'
#'@return The function will return a list of the groups defined in the 
#'tab_col column and whether they are monophyletic or not
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import ape
#'

tree_analyse <- function(tree_file, info_table, tab_col){
  
  Tree_1 <- tree_file
  classification <- info_table
  
  tip_labels <- Tree_1$tip.label
  
  
  Family_vec <- unique(classification[,tab_col])
  Families_monophyletic <- c()
  Families_non_monophyletic <- c()
  
  # Iterate over every family
  
  for (fam in Family_vec) {
    
    Family_names_vec <- rownames(classification)[which(classification[,tab_col] == fam)]
    if(is.monophyletic(Tree_1, Family_names_vec)){
      Families_monophyletic <- append(Families_monophyletic, fam)
    }else{
      Families_non_monophyletic <- append(Families_non_monophyletic, fam)
    }
    
  }

  Monophyly <- list(Families_monophyletic = Families_monophyletic, 
                    Families_non_monophyletic = Families_non_monophyletic)
  
  return(Monophyly)
}


#'\code{clade_stats}
#'
#'@details This function compiles information about groups of tips
#'
#'@param tree_file a tree object 
#'@param info_table a table with as many rows as the tips in the tree file
#'and corresponding names and at least one column to group them
#'@param tab_col The column of the info_table used to group the tree tips, 
#'species that belong to the same group should have the same name
#'
#'@return TIt returns a table of information on the given groups
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import ape
#'

clade_stats <- function(file, info_table, tab_col){
  
  Phy_Genes <- file
  classification <- info_table
  
  Clade_vec <- unique(classification[,tab_col])
  
  genes_per_tip <- c()
  for (i in 1:dim(Phy_Genes)[1]){
    genes_per_tip[i] <- dim(Phy_Genes)[2] - length(which(str_remove_all(Phy_Genes[i,], "-") == ""))
  }
  names(genes_per_tip) <- rownames(Phy_Genes)
  
  Familywise_genes_per_tip <- data.frame(matrix(NA, nrow = length(Clade_vec), ncol = 6), row.names = Clade_vec)
  for (fam in Clade_vec) {
    
  Clade_names_vec <- rownames(classification)[which(classification[,tab_col] == fam)]
   
  Familywise_genes_per_tip[fam,1] <- length(Clade_names_vec)
  Familywise_genes_per_tip_vec <- genes_per_tip[Clade_names_vec]
  Familywise_genes_per_tip_vec <- Familywise_genes_per_tip_vec[which(!is.na(Familywise_genes_per_tip_vec))]
  Familywise_genes_per_tip[fam,2] <- min(Familywise_genes_per_tip_vec, na.rm = TRUE)
  Familywise_genes_per_tip[fam,3] <- max(Familywise_genes_per_tip_vec, na.rm = TRUE)
  Familywise_genes_per_tip[fam,4] <- mean(Familywise_genes_per_tip_vec, na.rm = TRUE)
  Familywise_genes_per_tip[fam,5] <- sd(Familywise_genes_per_tip_vec, na.rm = TRUE)
  Familywise_genes_per_tip[fam,6] <- length(Familywise_genes_per_tip_vec)/length(Clade_names_vec)
  }
  
  colnames(Familywise_genes_per_tip) <- c("tips", "min_seq_per_tip", "max_seq_per_tip", "mean_seq_per_tip", "sd_seq_per_tip", "ratio seq/no_seq tips")
  return(Familywise_genes_per_tip)
}


#'\code{pick_sister_Clades}
#'
#'@details find sister clades in a tree according to the classification of a column 
#'
#'@param tree_file a tree object 
#'@param info_table a table with as many rows as the tips in the tree file
#'and corresponding names and at least one column to group them
#'@param tab_col The column of the info_table used to group the tree tips, 
#'species that belong to the same group should have the same value
#'
#'@return The function will return a 2 column matrix with each row representing 
#'a sister pair
#'@author Matteo Spatuzzi, 2022
#'
#'@import ape
#'
#'

pick_sister_Clades <- function(tree_file, info_table, tab_col){
  
  library(ape)
  library(stringr)
  
  Tree_1 <- tree_file
  classification <- info_table
  
  Phy_names <- rownames(info_table)
  Phy_names_inTree <- Phy_names[which(Phy_names %in% Tree_1$tip.label)]
  
  # if(length(Phy_names_inTree) = 0){
  #   stop("The info_table names do not match the tree_file tips")
  # }

  Tree_1_phy <- keep.tip(Tree_1,Phy_names_inTree)
  
  Family_vec <- unique(classification[,tab_col])
  
  # Dropping tips to 1 per family
  
  # tips_list_fam <- c()
  # 
  # for (fam in Family_vec) {
  #   
  #   Family_names_vec <- rownames(classification)[which(classification[,tab_col] == fam)]
  #   Family_names_vec <- Family_names_vec[which(Family_names_vec %in% Tree_1_phy$tip.label)]
  #   tips_list_fam <- append(tips_list_fam, Family_names_vec[1])
  # }
  # Tree_1_per_Family <- keep.tip(Tree_1_phy, tips_list_fam)
  # Tree_1_per_Family$tip.label <- Family_vec
  Tree_1_per_Family <- Tree_1_phy
  sister_clades <- matrix( ncol = 2)
  
  #Randomise
  
  viable_families <- sample(Family_vec, length(Family_vec))
  
  while ((length(viable_families) > 1)) {
    
    iter = 2
    search_complete <- FALSE
    family_1 <- viable_families[1]
    #print(paste("Family_1:", family_1))
    
    while ((search_complete == FALSE)) {
      #print(viable_families[iter])
      if(is.monophyletic(Tree_1_per_Family, c(rownames(classification)[classification[,tab_col] == family_1], rownames(classification)[classification[,tab_col] == viable_families[iter]]))){
        
        sister_clades <- rbind(sister_clades, c(family_1, viable_families[iter])[order(c(family_1, viable_families[iter]))])
        search_complete <- TRUE
        #print("success")
        viable_families <- viable_families[-which(viable_families == family_1 |viable_families == viable_families[iter])]
        
      }else{
        
        if(iter == length(viable_families)){
          
          search_complete <- TRUE
          viable_families <- viable_families[-which(viable_families == family_1)]
          #print("abort")
          
          
        }else{
          
          iter = iter +1
          
          #print("increase")
        }
      }
    }
    
  }
  
  sister_clades <- sister_clades[-1,]
  return(sister_clades)
}

#'\code{sister_clades_stats}
#'
#'@details returns a table of information on the given sister pairs
#'
#'@param sister_mat a nx2 matrix that lists the sister pairs
#'@param info_table a table with as many rows as the tips in the tree file
#'and corresponding names and at least one column to group them
#'@param tab_col The column of the info_table used to group the tree tips, 
#'species that belong to the same group should have the same name
#'
#'@return The function will return a 8 column matrix indicating for 
#'each sister pair the amount of tips, tips with at least one 
#'sampled gene and amount of total genes for each sister clade
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import ape

sister_clades_stats <- function(phy, sister_mat, info_table, tab_col){
  
  phy_dat <- phy
  sister_list <- sister_mat
  data <- data.frame(matrix(NA, ncol = 8, nrow = dim(sister_list)[1])) 
  classification <- info_table
  
  
  classification <- classification[which(rownames(classification) %in% rownames(phy_dat)),]
  
  Sister_clades_data <- c() 
  for (pair in 1:dim(sister_list)[1]){
    
    SistercladeN <- sister_list[pair,]
    
    SistercladeN <- append(SistercladeN, length(which(classification[,tab_col] == SistercladeN[1])))
    SistercladeN <- append(SistercladeN, length(which(classification[,tab_col] == SistercladeN[2])))
    SistercladeN <- append(SistercladeN, sum(sign(classification$Genes[which(classification[,tab_col] == SistercladeN[1])]), na.rm = TRUE))
    SistercladeN <- append(SistercladeN, sum(sign(classification$Genes[which(classification[,tab_col] == SistercladeN[2])]), na.rm = TRUE))
    SistercladeN <- append(SistercladeN, sum(classification$Genes[which(classification[,tab_col] == SistercladeN[1])], na.rm = TRUE))
    SistercladeN <- append(SistercladeN, sum(classification$Genes[which(classification[,tab_col] == SistercladeN[2])], na.rm = TRUE))
    
    
    Sister_clades_data <- rbind(Sister_clades_data, SistercladeN)
  }
  
  colnames(Sister_clades_data) = c("Sister 1", "Sister 2", "N tips 1", "N tips 2", "N sampled tips 1", "N sampled tips 2", "N genes 1", "N genes 2")
  row.names(Sister_clades_data) <- seq(1,dim(Sister_clades_data)[1])
  return(Sister_clades_data)
} 



#'\code{reduce_tips}
#'
#'@details for each sister pair reduces the amount of tips in the tree in the sister 
#'pair with more tips based on amount lost information to circumvent the 
#'node density effect
#'@param phy a dataframe with the multiple alignment sequences, 
#'rownames must correspond to the names of the tree tips
#'@param criteria vector of genes to use as criteria to remove tips, 
#'if FALSE all genes are used
#'@param sister_mat a nx2 matrix that lists the sister pairs
#'@param sister_data a nx8 matrix generated by the 
#'sister_clades_stats function from the sister_mat argument
#'@param info_table a table with as many rows as the tips in the tree file
#'and corresponding names and at least one column to group them
#'@param tab_col The column of the info_table used to group the tree tips
#'@param tree_file a tree object
#'
#'@return The function will return a 8 column matrix indicating for 
#'each sister pair the amount of tips, tips with at least one 
#'sampled gene and amount of total genes for each sister clade
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import ape



# phy = Phylo_genes
# criteria = con_list$mit
# sister_mat = sister_clades
# sister_data = sister_clade_data
# info_table = classification_amph
# tab_col = "Family"
# tree_file = Tree_1_amph

reduce_tips <- function(phy, 
                        criteria = FALSE, 
                        tip_threshold = FALSE, 
                        sister_mat, 
                        sister_data, 
                        info_table, 
                        tab_col, 
                        tree_file){
  
  Phy_Genes <- phy
  classification <- info_table
  Tree_1 <- tree_file
  group <- criteria
  sister_list <- sister_mat
  
  if( group == FALSE){
    
    Phy_Genes <- Phy_Genes
    
  }else{
    if(group %in% colnames(Phy_Genes)){
      
      Phy_Genes <- Phy_Genes[,group]
      Phy_Genes <- Phy_Genes[which(apply(Phy_Genes, 1, function(x) any(!(str_remove_all(x, "-") == "")))),]
      print("removed")
      
    }else{
      stop("Wrong value for attribute \"criteria\"")
    }
  }
  
  Phy_Genes_equal <- data.frame(matrix(NA, nrow = 0, ncol = dim(Phy_Genes)[2] + 2))
  
  for (npair in 1:dim(sister_list)[1]){
    
    
    sister_pair <- sister_list[npair,]
    tip_list <- list(
      rownames(classification)[which(classification[,tab_col] == sister_pair[1] & rownames(classification) %in% rownames(Phy_Genes))],
      rownames(classification)[which(classification[,tab_col] == sister_pair[2] & rownames(classification) %in% rownames(Phy_Genes))]
    )
    names(tip_list) <- sister_pair
    
    Ntips <- c(length(tip_list[[1]]), length(tip_list[[2]]))
    names(Ntips) = names(tip_list)
    
    max_tips <- max(Ntips)
    min_tips <- min(Ntips)
    if(min_tips == 0) next
    
    Phy_Genes_guide <- data.frame(Phy_Genes[which(rownames(Phy_Genes) %in% append(tip_list[[1]],tip_list[[2]])),])
    if(class(tip_threshold) == "numeric" & min_tips > tip_threshold){
     
      min_fam <- names(Ntips)[which(Ntips == min_tips)]
      Phy_Genes_reduce <- data.frame(Phy_Genes[which(rownames(Phy_Genes) %in% tip_list[[min_fam]]),])
      min_fam_ordered_tips <- order_by_gene(Phy_Genes_reduce, Phy_Genes_guide)
      tip_list[[min_fam]] <- tip_list[[min_fam]][tip_list[[min_fam]] %in% names(min_fam_ordered_tips)[1:tip_threshold]]
      min_tips <- tip_threshold
      Ntips[min_fam] <- tip_threshold
    }
    
    if(max_tips == min_tips){
      
      Family_reduce <- sister_pair[1]
      Family_keep <- sister_pair[2]
    
      tag_mat <- matrix(c(rep(npair, 2*min_tips), rep(Family_reduce, min_tips), rep(Family_keep, min_tips)), ncol = 2)
      Phy_Genes_equal <- rbind(Phy_Genes_equal, cbind(tag_mat,Phy_Genes[c(tip_list[[1]][1:max_tips],tip_list[[2]][1:min_tips]),]))
      
      
    }else{
      
      names(max_tips) <- names(Ntips)[which(Ntips == max_tips)]
      names(min_tips) <- names(Ntips)[which(Ntips == min_tips)]
      Family_reduce <- names(max_tips)
      tips_reduce <- tip_list[[Family_reduce]]
      Family_keep <- names(min_tips)
      tip_keep <- tip_list[[Family_keep]]
    
    
      Phy_Genes_reduce <- data.frame(Phy_Genes[which(rownames(Phy_Genes) %in% tips_reduce),])
      
      tip_scores_ordered <- order_by_gene(Phy_Genes_reduce, Phy_Genes_guide)
      tip_keep <- append(names(tip_scores_ordered[1:min(Ntips)]), tip_keep)
      
      tag_mat <- matrix(c(rep(npair, 2*min_tips), rep(Family_reduce, min_tips), rep(Family_keep, min_tips)), ncol = 2)
      
      diff_names <- setdiff(tip_keep, rownames(Phy_Genes))
      Phy_Genes_equal <- rbind(Phy_Genes_equal, cbind(tag_mat,Phy_Genes[tip_keep,]))
      row.names(Phy_Genes_equal)[which(str_detect(row.names(Phy_Genes_equal), "NA"))] <- diff_names
    }
  
  }
  
  colnames(Phy_Genes_equal)[1:2] <- c("Pair", "Family")
  return(Phy_Genes_equal)
}



order_by_gene <- function(dat, guide_dat){
  
  Guide_mat <- guide_dat
  Order_mat <- dat
  
  Guide_mat[] <- lapply(Guide_mat, function(x) as.numeric(str_remove_all(x, "-") != ""))
  Guide_vec <- apply(Guide_mat, 2, sum)
  
  Order_mat[] <- lapply(dat, function(x) as.numeric(str_remove_all(x, "-") != ""))
  tip_scores <- as.matrix(Order_mat) %*% Guide_vec
  
  return(tip_scores[order(tip_scores, decreasing = TRUE),1])
}


