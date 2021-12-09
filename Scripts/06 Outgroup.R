library(ape)
library(stringr)
WorkingDirectory <- "/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data"
setwd(WorkingDirectory)

Tree <- read.tree("amph_shl_new.tre")
classification <- read.csv("amph_shl_new_Classification.csv")
classification_names <- classification[,"Scientific.Name"]
classification_names <- str_replace(classification_names, " ", "_")
classification[,"Scientific.Name"] <- classification_names

Phy_genes_reduced <- read.csv("Phy_genes_reduced", row.names = 1)

Family_vec <- unique(classification["Family"])
Family_vec <- Family_vec[-which(Family_vec[,1] %in% c("Hominidae", "")),]

Phy_names_inTree <- rownames(Phy_genes_reduced)[which(rownames(Phy_genes_reduced) %in% Tree$tip.label)]
Tree_phy <- keep.tip(Tree,Phy_names_inTree)

##### Find outgroup #####

Subtrees <- c()
Subtrees_note <- matrix(ncol= 3, nrow = dim(sister_clades)[1])

for(pair in unique(Phy_genes_reduced$Pair)){
  
  families = unique(Phy_genes_reduced$Family[which(Phy_genes_reduced$Pair == pair)])
  search_complete <- FALSE
  Seq_threshold = 3
  iter = 1
  print(families)
  
  while(search_complete == FALSE){
    
    outgroup_index <- which(Tree_1_per_Family$tip.label %in% families) + c(-1*iter, iter)
    outgroups <- Tree_1_per_Family$tip.label[outgroup_index[which(outgroup_index > 0)]]
    if(any(is.na(outgroups))){
      outgroups <- outgroups[-which(is.na(outgroups) ==  TRUE)]
    }
    cophenetic_outgroups <- cophenetic.phylo(keep.tip(Tree_1_per_Family, tip = append(families, outgroups)))
    dist_outgroups <- which(cophenetic_outgroups[families[1],] == min(cophenetic_outgroups[outgroups,families[1]]))
    candidate <- names(dist_outgroups)
    print(candidate)
    
    Phy_Genes_reduce <- (Phy_Genes_update[classification$Scientific.Name[which(classification$Family == candidate)],])
    Phy_Genes_reduce[] <- lapply(Phy_Genes_reduce, function(x) ceiling(str_length(str_remove_all(x,"-"))/(str_length(str_remove_all(x,"-"))+1)))
    Phy_Genes_sum_col <- apply(Phy_Genes_reduce, 2, function(x) sum(x, na.rm = TRUE))
    
    if(sum(Phy_Genes_sum_col > 2) > 2){
      outgroup <- candidate
      search_complete <- TRUE
    }else{
      iter = iter +1
    }
  }
  print(search_complete)
  subtree_tips <- classification$Scientific.Name[which(classification$Family %in% families)]
  subtree_tips <- append(subtree_tips, classification$Scientific.Name[which(classification$Family %in% outgroup)])
  
  subtree_tips <- subtree_tips[subtree_tips %in% Tree_1$tip.label]
  subtree_tips <- subtree_tips[subtree_tips %in% rownames(Phy_genes_reduced)]
  
  Subtrees_note[pair,] <- append(families, outgroup)
  Subtrees[[pair]] <- keep.tip(Tree_1, subtree_tips)
}

Subtrees_note
 

   
#### Write sequences datasets for each sister pair subtree ####

system2("mkdir", args = "subtrees")
setwd(WorkingDirectory)
setwd("subtrees")

for(subtree in 1: length(Subtrees)){
  
  setwd(WorkingDirectory)
  setwd("subtrees")
  
  subtree_tips <- Subtrees[[subtree]]$tip.label
  Phy_subtree <- Phy_genes_reduced[which(  rownames(Phy_genes_reduced) %in% subtree_tips) ,]
  Phy_subtree_seq <- Phy_concat(Phy_subtree)
  rownames(Phy_subtree_seq) <- rownames(Phy_subtree)
  
  dir <- paste("subtree_", Subtrees_note[subtree,1], "_", Subtrees_note[subtree,2], sep = "")
  system2("mkdir", args = dir)
  setwd(dir)
  
  write.csv(Phy_subtree, file = dir, row.names = TRUE)
  write.tree(Subtrees[[subtree]], file = paste(dir, "_tree.tre", sep = ""))
  
  # Nuclear fasta
  string <- ""
  for(org in 1:dim(Phy_subtree_seq)[1]){
    string <- paste(string, ">",rownames(Phy_subtree_seq)[org], "\n",Phy_subtree_seq[org,2], "\n", sep = "" )
  }
  
  filename <- paste(dir,"_nuc.fasta",sep = "" )
  cat(string, file =  filename)
  
  # Mitochondrial fasta
  
  string <- ""
  for(org in 1:dim(Phy_subtree_seq)[1]){
    string <- paste(string, ">",rownames(Phy_subtree_seq)[org], "\n",Phy_subtree_seq[org,1], "\n", sep = "" )
  }
  
  filename <- paste(dir,"_mit.fasta",sep = "" )
  cat(string, file =  filename)
  
}
  



  


