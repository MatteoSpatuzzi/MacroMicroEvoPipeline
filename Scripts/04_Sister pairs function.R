# Actual Tree, must be of class "phylo"
myTree <- Tree_1

#Extract relevant information from Tree object
full.depth <- dist.nodes(Tree_1)[1,(length(Tree_1$tip.label)+1)]
max_depth.edglength <- max(node.depth.edgelength(myTree))
NNodes <-  length(node.depth.edgelength(myTree))
tips_available <- myTree$tip.label
Ntips <- length(myTree$tip.label)

#Set Distance to "cut" tree at
Distance <- 50

# Determine clade roots based on distance
depth.edglength_myTree <- node.depth.edgelength(myTree)
depth.edglength_myTree <- rbind(c(1:NNodes), depth.edglength_myTree)
depth.edglength_myTree_internal <- depth.edglength_myTree[,(length(tips_available)+1):NNodes]
Clade_roots <- matrix(depth.edglength_myTree_internal[,which((full.depth - depth.edglength_myTree_internal[2,]) <= Distance)], nrow = 2)

#Prepare to assign tips to roots
Clades <- vector(mode = "list", length = length(Clade_roots[1,]))
names(Clades) <- Clade_roots[1,]
Free_tips <- c()

#Check closest root for every tip within distance
for (tip in 1:Ntips){
  
  path <- nodepath(myTree, from = tip, to = Ntips+1 )
  root <- path[max(which(path %in% Clade_roots[1,]))]
  if(is.na(root)){
    Free_tips <- append(Free_tips, tip)
  }else{
    Clades[[paste(root)]] <- append(Clades[[paste(root)]], tip)
  }
  
  print(tip)
  print(root)
}

# remove empty clades
Clades <- Clades[-which(sapply(Clades, is.null))]
#Reduce to one representative tip per clade

Clades_one_per_Clade <- lapply(Clades, min)
Clades_one_per_Clade_vec <- unlist(Clades_one_per_Clade)
Clades_one_per_Clade_mat <- matrix(as.integer(c(names(Clades_one_per_Clade_vec), Clades_one_per_Clade_vec)), byrow = FALSE, ncol = 2)
Tree_1_per_Clade <- keep.tip(Tree_1, Tree_1$tip.label[Clades_one_per_Clade_mat[,2]])
Tree_1_per_Clade$tip.label <- as.character(Clades_one_per_Clade_mat[,1])

#Prepare for sister pairs
viable_clades <- as.character(Clades_one_per_Clade_mat[,1])
sister_clades <- matrix(ncol = 2)

#Assign sister peir values 
while ((length(viable_clades) > 1)) {
  
  iter = 2
  search_complete <- FALSE
  clade_1 <- viable_clades[1]
  #print(paste("clade_1:", clade_1)) Could start 
  
  while ((search_complete == FALSE)) {
    #print(viable_clades[iter])
    if(is.monophyletic(Tree_1_per_Clade, c(clade_1, viable_clades[iter]))){
      
      sister_clades <- rbind(sister_clades, c(clade_1, viable_clades[iter])[order(c(clade_1, viable_clades[iter]))])
      search_complete <- TRUE
      
      
      #print("success")
      viable_clades <- viable_clades[-which(viable_clades == clade_1 |viable_clades == viable_clades[iter])]
      
    }else{
      
      if(iter == length(viable_clades)){
        
        search_complete <- TRUE
        viable_clades <- viable_clades[-which(viable_clades == clade_1)]
        #print("abort")
        
        
      }else{
        
        iter = iter +1
        
        #print("increase")
      }
    }
  }
  
}

print(sister_clades)
