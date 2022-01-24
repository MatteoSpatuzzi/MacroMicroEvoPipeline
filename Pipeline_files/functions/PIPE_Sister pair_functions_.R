#### Tree File #### 

tree_file <- Tree_1
info_table <- classification

Tree_analyse <- function(tree_file, info_table){
  
  Tree_1 <- tree_file
  classification <- info_table
  
  tip_labels <- Tree_1$tip.label
  
  ##### Monophyly #####
  # checking whether families are monophyletic
  # Establish a list of families and exclude Homo sapiens and BLANK from the list
  
  Family_vec <- unique(classification["Family"])
  Family_vec <- Family_vec[1:(dim(Family_vec)[1]-2),]
  Families_monophyletic <- c()
  Families_non_monophyletic <- c()
  
  # Iterate over every family
  
  for (fam in Family_vec) {
    
    Family_names_vec <- classification[which(classification[,"Family"] == fam),"Scientific.Name"]
    if(is.monophyletic(Tree_1, Family_names_vec)){
      Families_monophyletic <- append(Families_monophyletic, fam)
    }else{
      Families_non_monophyletic <- append(Families_non_monophyletic, fam)
    }
    
  }
  

  # checking whether genuses are monophyletic
  # Establish a list of genuses and exclude Homo sapiens and BLANK from the list
  Genus_vec <- unique(classification["Constraint"])
  Genus_vec<- Genus_vec[1:(dim(Genus_vec)[1]-2),]
  Genus_monophyletic <- c()
  Genus_non_monophyletic <- c()
  
  # Iterate over every genus
  for (genus in Genus_vec) {
    
    Genus_names_vec <- classification[which(classification[,"Constraint"] == genus),"Scientific.Name"]
    if(is.monophyletic(Tree_1, Genus_names_vec)){
      Genus_monophyletic <- append(Genus_monophyletic, genus)
    }else{
      Genus_non_monophyletic <- append(Genus_non_monophyletic, genus)
    }
    
  }

  Monophyly <- list(Genus_monophyletic = Genus_monophyletic, 
                    Genus_non_monophyletic = Genus_non_monophyletic, 
                    Families_monophyletic = Families_monophyletic, 
                    Families_non_monophyletic = Families_non_monophyletic)
  
  return(Monophyly)
  
  rm(genus, Genus_names_vec)
}

# This function determines sister clades according to family classificatio
# 
# tree_file is the phylognetic tree object
# info_table is the supplementary information table with the phylogenetic classsification of the tree
# 
# It returns a nx2 matrix of sister clades
Pick_sister_Clades <- function(tree_file, info_table){
  
  library(ape)
  library(stringr)
  
  Tree_1 <- tree_file
  classification <- info_table
  
  Phy_names <- classification$Scientific.Name
  Phy_names_inTree <- Phy_names[which(Phy_names %in% Tree_1$tip.label)]
  
  # if(length(Phy_names_inTree) = 0){
  #   stop("The info_table names do not match the tree_file tips")
  # }
  # 
  Tree_1_phy <- keep.tip(Tree_1,Phy_names_inTree)
  
  Family_vec <- unique(classification["Family"])
  Family_vec<- Family_vec[1:(dim(Family_vec)[1]-2),]
  
  # Dropping tips to 1 per family
  
  tips_list_fam <- c()
  
  for (fam in Family_vec) {
    
    Family_names_vec <- classification[which(classification[,"Family"] == fam),"Scientific.Name"]
    Family_names_vec <- Family_names_vec[which(Family_names_vec %in% Tree_1_phy$tip.label)]
    tips_list_fam <- append(tips_list_fam, Family_names_vec[1])
  }
  Tree_1_per_Family <- keep.tip(Tree_1_phy, tips_list_fam)
  Tree_1_per_Family$tip.label <- Family_vec

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
      if(is.monophyletic(Tree_1_per_Family, c(family_1, viable_families[iter]))){
        
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


# This function determines sister clades according to family classification
# 
# tree_file is the phylognetic tree object
# info_table is the supplementary information table with the phylogenetic classification of the tree
# 
# It returns a table of information on the sister pairs

sister_list <- sister_clades
info_table <- classification

Sister_clades_stats<- function(sister_list, info_table){
  
  data <- data.frame(matrix(NA, ncol = 8, nrow = dim(sister_list)[1])) 
  classification <- info_table
  
  Sister_clades_data <- c() 
  for (pair in 1:dim(sister_list)[1]){
    
    SistercladeN <- sister_list[pair,]
    
    SistercladeN <- append(SistercladeN, length(classification$Scientific.Name[which(classification$Family == SistercladeN[1])]))
    SistercladeN <- append(SistercladeN, length(classification$Scientific.Name[which(classification$Family == SistercladeN[2])]))
    SistercladeN <- append(SistercladeN, sum(sign(classification$Genes[which(classification$Family == SistercladeN[1])]), na.rm = TRUE))
    SistercladeN <- append(SistercladeN, sum(sign(classification$Genes[which(classification$Family == SistercladeN[2])]), na.rm = TRUE))
    SistercladeN <- append(SistercladeN, sum(classification$Genes[which(classification$Family == SistercladeN[1])], na.rm = TRUE))
    SistercladeN <- append(SistercladeN, sum(classification$Genes[which(classification$Family == SistercladeN[2])], na.rm = TRUE))
    
    
    Sister_clades_data <- rbind(Sister_clades_data, SistercladeN)
  }
  
  colnames(Sister_clades_data) = c("Sister 1", "Sister 2", "N tips 1", "N tips 2", "N sampled tips 1", "N sampled tips 2", "N genes 1", "N genes 2")
  row.names(Sister_clades_data) <- seq(1,dim(Sister_clades_data)[1])
  return(Sister_clades_data)
} 

# This function compiles information about each family or Genus
# 
# tree_file is the phylognetic tree object
# info_table is the supplementary information table with the phylogenetic classsification of the tree
# Mode determines whether the information is gathered about Families of Genuses
#   F: family
#   G: Genus
#
# It returns a table of information on the families or genuses
Clade_stats <- function(file, info_table, mode = "F"){

  Phy_Genes <- file
  classification <- info_table
  
  if(mode == "F"){
    Clade_vec <- unique(classification["Family"])[,1]
    
  }else{
    if(mode == "G"){
      Clade_vec <- unique(classification["Constraint"])[,1]
      
    }else{
      stop("Wrong mode")
    }
  }
  
  genes_per_tip <- c()
  for (i in 1:dim(Phy_Genes)[1]){
    genes_per_tip[i] <- dim(Phy_Genes)[2] - length(which(str_remove_all(Phy_Genes[i,], "-") == ""))
  }
  names(genes_per_tip) <- rownames(Phy_Genes)
  
  Familywise_genes_per_tip <- data.frame(matrix(NA, nrow = length(Clade_vec), ncol = 6), row.names = Clade_vec)
  for (fam in Clade_vec) {
    
    if(mode == "F"){
      Clade_names_vec <- classification[which(classification[,"Family"] == fam),"Scientific.Name"]
    }else{
      if(mode == "G"){
        Clade_names_vec <- classification[which(classification[,"Constraint"] == fam),"Scientific.Name"]
      }
    }
    
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


# This function even the amount of tips within clades by removing the excessive tips based on seqeunce 
# information to minimize the loss of sequences
#
# file is the output of the Edit_phylo() function
# group indicates whethere this should be done for every gene or only certain genes according to the concat_list argument
# concat_list must be a list vectors, the name of each element will become the colname of 
# the concatenated sequence, the values correspont to the colnames of the input genes to group under 
# said name
# E.g. concat_list <- list(
# mit = c("cytb","ND1","ND2"),
# nuc = c("BDNF", "CXCR4","H3A","NCX1","POMC", "RAG1", "RHOD", "RHOD", "SIA", "SLC8A3", "TYR")
# )
# sister list is the output of the Pick_sister_Clades() function
# sister_data is the output of the Sister_clades_stats() function
# info_table is the supplementary information table with the phylogenetic classsification of the tree
# tree_file is the phylognetic tree object
# 
# It returns a dataframe similar to the "file" argument but with fewer rows and potentially fewer colums

file <- Phylo_genes
group <- "mit" 
concat_list <- con_list 
sister_list <- sister_clades
sister_data <- sister_clade_data
info_table <- classification
tree_file <- Tree_1

Reduce_tips <- function(file, group = FALSE, concat_list = FALSE, sister_list, sister_data, info_table, tree_file){
  
  Phy_Genes <- file
  classification <- info_table
  Tree_1 <- tree_file
  
  if( group == FALSE){
    
    Phy_Genes <- Phy_Genes
    
  }else{
    if(group %in% names(concat_list)){
      
      Phy_Genes <- Phy_Genes[,concat_list[[group]]]
      Phy_Genes <- Phy_Genes[which(apply(Phy_Genes, 1, function(x) any(!str_remove_all(x, "-") == ""))),]
      
    }else{
      stop("Wrong value for attribute \"group\"")
    }
  }
  
  Phy_Genes_equal <- data.frame(matrix(NA, nrow = 0, ncol = dim(Phy_Genes)[2] + 2))
  for (npair in 1:dim(sister_list)[1]){
    
    sister_pair <- sister_list[npair,]
    
    tip_list <- list(
      
      classification$Scientific.Name[which(classification$Family == sister_pair[1] & classification$Scientific.Name %in% rownames(Phy_Genes))],
      classification$Scientific.Name[which(classification$Family == sister_pair[2] & classification$Scientific.Name %in% rownames(Phy_Genes))]
    )
    names(tip_list) <- sister_pair
    
    tips <- c(length(tip_list[[1]]), length(tip_list[[2]]))
    names(tips) = names(tip_list)
    max_tips <- max(tips)
    min_tips <- min(tips)
    
    
    if(max(tips) == min(tips)){
      
      Family_reduce <- sister_pair[1]
      tips_reduce <- classification$Scientific.Name[which(classification$Family %in% Family_reduce)]
      Family_keep <- sister_pair[2]
      tip_keep <- classification$Scientific.Name[which(classification$Family %in% Family_keep)]
      
    }else{
      
      names(max_tips) <- names(tips)[which(tips == max_tips)]
      names(min_tips) <- names(tips)[which(tips == min_tips)]
      Family_reduce <- names(max_tips)
      tips_reduce <- tip_list[[Family_reduce]]
      Family_keep <- names(min_tips)
      tip_keep <- tip_list[[Family_keep]]
    }
    
    Phy_Genes_reduce <- data.frame(Phy_Genes[which(rownames(Phy_Genes) %in% tips_reduce),])
    Phy_Genes_reduce[] <- lapply(Phy_Genes_reduce, function(x) ceiling(str_length(str_remove_all(x,"-"))/(str_length(str_remove_all(x,"-"))+1)))
    Phy_Genes_sum_col <- apply(Phy_Genes_reduce, 2, sum)
    
    tag_mat <- matrix(c(rep(npair, 2*min_tips), rep(Family_reduce, min_tips), rep(Family_keep, min_tips)), ncol = 2)
    
    tip_scores <- c(rep(0, length(tips_reduce)))
    names(tip_scores) <- tips_reduce
    
    for (i in rownames(Phy_Genes_reduce)) {
      tip_scores[i] <- sum(Phy_Genes_sum_col[Phy_Genes_reduce[which(rownames(Phy_Genes_reduce) == i),] == 1])
    } 
    
    tip_scores_ordered <- tip_scores[order(tip_scores, decreasing = TRUE)]
    tip_keep <- append(tip_keep, names(tip_scores_ordered[1:min(tips)]))
    
    diff_names <- setdiff(tip_keep, rownames(Phy_Genes))
    Phy_Genes_equal <- rbind(Phy_Genes_equal, cbind(tag_mat,Phy_Genes[tip_keep,]))
    row.names(Phy_Genes_equal)[which(str_detect(row.names(Phy_Genes_equal), "NA"))] <- diff_names
  
  }
  
  colnames(Phy_Genes_equal)[1:2] <- c("Pair", "Family")
  return(Phy_Genes_equal)
}


