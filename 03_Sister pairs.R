#### Tree File #### 

#### Cleanup ####

library(ape)
library(stringr)
WorkingDirectory <- "/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data"
setwd(WorkingDirectory)

Tree_1 <- read.tree("Tree_1")
tip_labels <- Tree_1$tip.label

# Comparing tip names and taxonomy names
classification <- read.csv("amph_shl_new_Classification.csv")
classification_names <- classification[,"Scientific.Name"]
classification_names <- str_replace(classification_names, " ", "_")
classification[,"Scientific.Name"] <- classification_names
mismatch <- which(!tip_labels %in% classification_names)
print(mismatch)

# Comparing tip labels across the trees
diff <- 0
for (tree in 2:1000) {
  diff <- diff + length(setdiff(Tree_1$tip.label, Tree_1$tip.label))
}
print(diff)

rm(classification_names, mismatch, diff, tree)

##### Monophyly #####
# checking whether families are monophyletic
# Establish a list of families and exclude Homo sapiens and BLANK from the list
Family_vec <- unique(classification["Family"])
Family_vec<- Family_vec[1:(dim(Family_vec)[1]-2),]
Families_monophyletic <- c()
Families_non_monophyletic <- c()

# Iterate over every family
for (fam in Family_vec) {

  Family_names_vec <- classification[which(classification[,"Family"] == fam),"Scientific.Name"]
  if(is.monophyletic(Tree_1, Family_names_vec)){
    print(paste("Family ", fam, " is monophyletic.", sep = ""))
    Families_monophyletic <- append(Families_monophyletic, fam)
  }else{
    print(paste("Error: Family ", fam, " is not monophyletic.", sep = ""))
    Families_non_monophyletic <- append(Families_non_monophyletic, fam)
  }

}
Families_monophyletic
Families_non_monophyletic

rm(fam, Family_names_vec)

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
    print(paste("Genus ", genus, " is monophyletic.", sep = ""))
    Genus_monophyletic <- append(Genus_monophyletic, genus)
  }else{

    print(paste("Error: Genus ", genus, " is not monophyletic.", sep = ""))
    Genus_non_monophyletic <- append(Genus_non_monophyletic, genus)
  }

}
Genus_monophyletic
Genus_non_monophyletic

rm(genus, Genus_names_vec)


# #### Sister pair analysis #####
# #We need to create as many sets of independent pairs of monophyletic sister clades

# Reducing tree to tips that we have a sequence for
# remove scientific names that do not appear in the phylogenetic file
# This way we can know which families we can not use for sister pairs

library(phylotools)

Phy <- read.csv("Phy_Genes_update.csv")
Phy_names <- Phy[,1]
rm(Phy)

Phy_names_inTree <- Phy_names[- which(!Phy_names %in% Tree_1$tip.label)]
Tree_1_phy <- keep.tip(Tree_1,Phy_names_inTree)
length(Tree_1_phy$tip.label)

classification_phy <- classification[which(classification$Scientific.Name %in% Phy_names_inTree),]

# rm(Phy_names_inTree)

# Dropping tips to 1 per genus
tips_list_genus <- c()
for (genus in Genus_vec) {
  
  Genus_names_vec <- classification_phy[which(classification_phy[,"Constraint"] == genus),"Scientific.Name"]
  tips_list_genus <- append(tips_list_genus, Genus_names_vec[1])
}
tips_list_genus <- tips_list_genus[which(tips_list_genus %in% Tree_1_phy$tip.label)]
Tree_1_phy_per_Genus <- keep.tip(Tree_1_phy,tips_list_genus)
Tree_1_phy_per_Genus$tip.label <- Genus_vec

# Dropping tips to 1 per family

tips_list_fam <- c()

for (fam in Family_vec) {
  
  Family_names_vec <- classification_phy[which(classification_phy[,"Family"] == fam),"Constraint"]
  tips_list_fam <- append(tips_list_fam, Family_names_vec[1])
}

Tree_1_phy_per_Family <- keep.tip(Tree_1_phy_per_Genus,tips_list_fam)
Tree_1_phy_per_Family$tip.label <- Family_vec

rm(tips_list_fam, fam, Family_names_vec, Genus_names_vec, genus, tips_list_genus)

# All families are still present
# Due to aberrant tree architecture, the monophyletic sister pairs are found based on the original Tree

#Produce an nx2 matrix of pairs
#NOTE: if two families are monophyletic, they will not be monophyletic with any other family -> recursive algorithm
#Therefore there should only be 1 possible way to arrange the clades in sister pairs (assuming the clades are fixed)
#To ensure this the algorithm is run 100 times, randomising the order of the clades

viable_families <- sample(Family_vec, length(Family_vec))
sister_clades <- matrix(ncol = 2)

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
rm(search_complete, iter, family_1, viable_families)


##### Sister Clades Data #####

Sister_clades_data <- data.frame(matrix(NA, ncol = 8, nrow = 0)) 

for (pair in 1:dim(sister_clades)[1]){

  SistercladeN <- sister_clades[pair,]
  
  SistercladeN <- append(SistercladeN, length(classification$Scientific.Name[which(classification$Family == SistercladeN[1])]))
  SistercladeN <- append(SistercladeN, length(classification$Scientific.Name[which(classification$Family == SistercladeN[2])]))
  SistercladeN <- append(SistercladeN, sum(classification$Data[which(classification$Family == SistercladeN[1])]))
  SistercladeN <- append(SistercladeN, sum(classification$Data[which(classification$Family == SistercladeN[2])]))
  SistercladeN <- append(SistercladeN, sum(classification$Genes[which(classification$Family == SistercladeN[1])]))
  SistercladeN <- append(SistercladeN, sum(classification$Genes[which(classification$Family == SistercladeN[2])]))
  
  
  Sister_clades_data <- rbind(Sister_clades_data, SistercladeN)
}

colnames(Sister_clades_data) = c("Sister 1", "Sister 2", "N tips 1", "N tips 2", "N sampled tips 1", "N sampled tips 2", "N genes 1", "N genes 2")

write.csv(Sister_clades_data, "Sister_clades")

##### Family stats #####
Phy_Genes_update <- read.csv("Phy_Genes_update.csv")
rownames(Phy_Genes_update) <- Phy_Genes_update[,1]

genes_per_tip <- c()

for (i in 1:dim(Phy_Genes_update)[1]){
  genes_per_tip[i] <- 15 - length(which(str_remove_all(Phy_Genes_update[i,2:16], "-") == ""))
}
names(genes_per_tip) <- Phy_Genes_update[,1]

Familywise_genes_per_tip <- data.frame(matrix(NA, nrow = length(Family_vec), ncol = 6), row.names = Family_vec)
for (fam in Family_vec) {
  
  Family_names_vec <- classification[which(classification[,"Family"] == fam),"Scientific.Name"]
  Familywise_genes_per_tip[fam,1] <- length(Family_names_vec)
  Familywise_genes_per_tip_vec <- genes_per_tip[Family_names_vec]
  Familywise_genes_per_tip_vec <- Familywise_genes_per_tip_vec[which(!is.na(Familywise_genes_per_tip_vec))]
  Familywise_genes_per_tip[fam,2] <- min(Familywise_genes_per_tip_vec, na.rm = TRUE)
  Familywise_genes_per_tip[fam,3] <- max(Familywise_genes_per_tip_vec, na.rm = TRUE)
  Familywise_genes_per_tip[fam,4] <- mean(Familywise_genes_per_tip_vec, na.rm = TRUE)
  Familywise_genes_per_tip[fam,5] <- sd(Familywise_genes_per_tip_vec, na.rm = TRUE)
  Familywise_genes_per_tip[fam,6] <- length(Familywise_genes_per_tip_vec)/length(Family_names_vec)
}

colnames(Familywise_genes_per_tip) <- c("tips", "min_seq_per_tip", "max_seq_per_tip", "mean_seq_per_tip", "sd_seq_per_tip", "ratio seq/no_seq tips")

write.csv(Familywise_genes_per_tip, "Family_stats")

##### Reduce tips #####

Phy_Genes_update <- read.csv("Phy_Genes_update.csv")
rownames(Phy_Genes_update) <- Phy_Genes_update[,1]
Phy_Genes_update <-  Phy_Genes_update[,-c(1)]
Phy_genes_update_reduced <- data.frame(matrix(NA, nrow = 0, ncol = dim(Phy_Genes_update)[2] + 2))
sister_clades_1 <- sister_clades

for (npair in 1:dim(sister_clades)[1]){
  
  sister_pair <- sister_clades[npair,]
  tips <- as.integer(Sister_clades_data[npair,c(3,4)])
  max_tips <- max(tips)
  min_tips <- min(tips)
  
  if(max(tips) == min(tips)){
    
    Family_reduce <- sister_pair[1]
    tips_reduce <- classification$Scientific.Name[which(classification$Family %in% Family_reduce)]
    Family_keep <- sister_pair[2]
    tip_keep <- classification$Scientific.Name[which(classification$Family %in% Family_keep)]
   
    
  }else{
   
    Family_reduce <- sister_pair[which(Sister_clades_data[npair,c(3,4)] == max_tips)]
    tips_reduce <- classification$Scientific.Name[which(classification$Family %in% Family_reduce)]
    Family_keep <- sister_pair[which(Sister_clades_data[npair,c(3,4)] == min_tips)]
    tip_keep <- classification$Scientific.Name[which(classification$Family %in% Family_keep)]
  }
  
  Phy_Genes_reduce <- Phy_Genes_update[which(rownames(Phy_Genes_update) %in% tips_reduce),]
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
  
  diff_names <- setdiff(tip_keep, rownames(Phy_Genes_update))
  Phy_genes_update_reduced <- rbind(Phy_genes_update_reduced, cbind(tag_mat,Phy_Genes_update[tip_keep,]))
  row.names(Phy_genes_update_reduced)[which(str_detect(row.names(Phy_genes_update_reduced), "NA"))] <- diff_names
}

colnames(Phy_genes_update_reduced)[1:2] <- c("Pair", "Family")
write.csv(Phy_genes_update_reduced, "Phy_genes_reduced")  
  

