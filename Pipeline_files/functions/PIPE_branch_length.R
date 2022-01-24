## File: PIPE_branch_length.R

## Usage: subtree.mean.branch.length <- phylo.average.brlen (subtree)

## Parameters: This function takes a tree in phylo format (from the ape package). This is usually a subtree with a root edge. 
## Value: A double representing the phylogenetic average of branch lengths for the clade.

## Example: 

## test.tree <- read.tree(read.tree(text= "(((A:1, B:1):1, C:1):1, D:1);"))
## plot(test.tree)
## test.subtree <- drop.tip(test.tree, "D", root.edge = 1)
## plot(test.subtree, root.edge=T)
## av.blen <- phylo.average.brlen(test.subtree)
## av.blen # should be 2.5


require(ape)

phylo.average.brlen <- function (phy)
{
  
  po <- reorder.phylo(phy, order = "postorder")
  
  br <- array(0, dim = Nnode(phy) + Ntip(phy))
  
  for (i in 1:Nedge(phy))
  {
    br[po$edge[i, 1]] <- br[po$edge[i, 1]] + (po$edge.length[i] + br[po$edge[i, 2]])/2
  }
  
  root.edge.length = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
  singleton = ifelse(Ntip(phy) == 1, 2, 1) 
  
  return (br[Ntip(phy) + 1] * singleton + root.edge.length)
  
}


## Usage: subtree.mean.branch.length <- branch_length_treefiles (dir, Subtrees_note_list)

## Parameters: This function takes a directory with the subtrees from PAML and the output of the Find_Outgroups() function
## Value: a list of average branch lengths of the tree files

## Example: 

## test.tree <- read.tree(read.tree(text= "(((A:1, B:1):1, C:1):1, D:1);"))
## plot(test.tree)
## test.subtree <- drop.tip(test.tree, "D", root.edge = 1)
## plot(test.subtree, root.edge=T)
## av.blen <- phylo.average.brlen(test.subtree)
## av.blen # should be 2.5
  dir <- "dS"
  Subtrees_note_list <- Subtrees_note_list_mit
  dat_tree <- dated_Tree
  
  Phy_genes_reduced <- Phy_genes_reduced_mit
  branch_length_treefiles <- function(dir = getwd(), dat_tree = getwd(), Subtrees_note_list, Phy_genes_reduced){
  
  #Set directory 
  library(ape)
  WorkingDirectory <- dir
  setwd(WorkingDirectory)
  
  #get filenames of all subtrees from PAML
  files <- system2("ls", stdout = TRUE)
  # Filter filenames and read the tree files
  tree_files <- files[which(substr(files,str_length(files)-3 ,str_length(files) ) == ".tre")]
  Tree_list <- lapply(tree_files, FUN = read.tree)
  
  # Get dated tree file
  dated_Tree <- dat_tree
  
  #Iterate over every sister pair  
  # Initialize empty table to fill with results
  br_res <- data.frame(matrix(nrow = dim(Subtrees_note_list[[2]])[1], ncol = 7))
  for (tree in 1:dim(Subtrees_note_list[[2]])[1]){
    
    #Initialize empty result vector for iteration
    br_vec <- c()
   
    #Extract sister pair names 
    sister_pair <- Subtrees_note_list[[2]][tree,c(1,2)]
    
    #Extract equivalent tree from the tree files as the indexes might be different
    tree_indexes <- lapply(sister_pair, function(x) which(!is.na(str_locate(tree_files, x))[,1]))
    tree_index <- intersect(tree_indexes[[1]],tree_indexes[[2]])
    
    if(length(tree_index) == 1){
      
      Tree <- Tree_list[[tree_index]]
      
      ### Extract general information about the whole pair
      # Tree_length
      Tree$edge.length[which( Tree$edge[,1] == length(Tree$tip.label)+1)]
      length(Tree$tip.label)
      
      # Remove outgroup from subtree
      # Establish the number of the root node in reduced subtree  
      sister_pair_tips <- rownames(Phy_genes_reduced[which(Phy_genes_reduced[,2] %in% sister_pair),])
      Tree <- keep.tip(Tree, sister_pair_tips)
      root_id <- length(Tree$tip.label)+1
      
      #Establish time at which point the two sister diverged by averaging across dated trees
      lengths <- c(rep(NA, length(dated_Tree)/100))
      for (i in 1:((length(dated_Tree))/100)) {
        
        dated_Tree_red <- keep.tip(dated_Tree[[i]], sister_pair_tips)
        lengths[i] <- dist.nodes(dated_Tree_red)[1,length(dated_Tree_red$tip.label)+1]
        print(i)
      }
      length_avg <- mean(lengths)
      
      ### Extract general information about the individual sisters
      ## Sister 1
      
      # Trim subtree to only sister 1
      sister_pair_tips1 <- rownames(Phy_genes_reduced[which(Phy_genes_reduced[,2] == sister_pair[1]),])
      Tree_red1 <- keep.tip(Tree, sister_pair_tips1)
      
      # Calculate average branch length unless there is only one tip, just take the root to tip distance
      if(length(Tree_red1$tip.label) == 1){
        
        br_len1 <- Tree$edge.length[which(Tree$edge[,1] == root_id & Tree$tip.label %in% sister_pair_tips1)]
      }else{
        root_id1 <- getMRCA(Tree,tip = sister_pair_tips1)
        root_len1 <- Tree$edge.length[which(Tree$edge[,1] == root_id & Tree$edge[,2] == root_id1)]
        br_len1 <- phylo.average.brlen(Tree_red1) + root_len1
      }
      
      pair1_data <- c(sister_pair[1], length(which(classification$Family == sister_pair[1])), br_len1)
      
      ## Sister 2
      
      # Trim subtree to only sister 2
      sister_pair_tips2 <- rownames(Phy_genes_reduced[which(Phy_genes_reduced[,2] == sister_pair[2]),])
      Tree_red2 <- keep.tip(Tree, sister_pair_tips2)
      
      # Calculate average branch length unless there is only one tip, just take the root to tip distance
      if(length(Tree_red2$tip.label) == 1){
        
        br_len2 <- Tree$edge.length[which(Tree$edge[,1] == root_id & Tree$tip.label %in% sister_pair_tips2)]
      }else{
        root_id2 <- getMRCA(Tree,tip = sister_pair_tips2)
        root_len2 <- Tree$edge.length[which(Tree$edge[,1] == root_id & Tree$edge[,2] == root_id2)]
        br_len2 <- phylo.average.brlen(Tree_red2) + root_len2
      }
      
      pair2_data <- c(sister_pair[2], length(which(classification$Family == sister_pair[2])), br_len2)
      
      # Create final result vector
      # order the two sisters pair by average branch length in decraesing order
      
      pair_data <- cbind(pair1_data, pair2_data)
      order <- order(c(br_len1, br_len2), decreasing = TRUE)
      pair_data_vec <- append(pair_data[,order[1]], pair_data[,order[2]])
      pair_data_vec <- append(pair_data_vec, length_avg)
      
      #Insert in result dataframe
      br_res[tree,] <- pair_data_vec
      rownames(br_res)[tree] <- paste(sister_pair[1], sister_pair[2], sep = "_")
    }else{
      if(length(tree_index) == 0){
        
      }else{
          
          stop("Multiple files per sister pair")
        }
      }
      
    }
    
    
  
  
  colnames(br_res) <- c("Sister_1_Name", "Tips_in_Sister_1", "Sister_1_avg_br_len", "Sister_2_Name", "Tips_in_Sister_2", "Sister_2_avg_br_len", "Pair_MRCA_AgeA" )
  br_res[,2] <- as.numeric(br_res[,2])
  br_res[,3] <- as.numeric(br_res[,3])
  br_res[,5] <- as.numeric(br_res[,5])
  br_res[,6] <- as.numeric(br_res[,6])
  br_res[,7] <- as.numeric(br_res[,7])
  return(br_res)
  
}

contrast_calc <- function(dat, subs_vs_age = FALSE ){

  contrast_Nr_species <- apply(dat, 1, function(x) (log(as.numeric(x[2])) - log(as.numeric(x[5])))/sqrt(as.numeric(x[7])) )
  contrast_branch_length <- apply(dat, 1, function(x) (log(as.numeric(x[3])) - log(as.numeric(x[6])))/sqrt(as.numeric(x[7])) )
  LM <-lm(contrast_Nr_species ~ contrast_branch_length)
 
  
  L <- list(Species_contrast = contrast_Nr_species, Substitution_contrast = contrast_branch_length, Species_vs_Substitution = LM )

  if(subs_vs_age == TRUE){
    
    LM <-lm(abs(contrast_branch_length) ~ sqrt(dat[,7]))
    L[[4]] <- LM
  }
  
  if(subs_vs_age == TRUE){
    print("Oh yeah")
    LM <-lm(abs(contrast_Nr_species) ~ dat[,7])
    L[[5]] <- LM
  }
  return(L)
}
                        