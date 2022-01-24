#### Concatenate ####

# This function concatenates gene sequences 
#
# File should be a dataframe with columns equivalent to genes and rows equivalent to organisms
#
# concat_list must be a list vectors, the name of each element will become the colname of 
# the concatenated sequence, the values correspont to the colnames of the input genes to group under 
# said name
# E.g. concat_list <- list(
# mit = c("cytb","ND1","ND2"),
# nuc = c("BDNF", "CXCR4","H3A","NCX1","POMC", "RAG1", "RHOD", "RHOD", "SIA", "SLC8A3", "TYR")
# )
#
# Colnames of File must coincide with the concat_list values


Phy_concat <- function(file, concat_list){
  
  library(stringr)
  
  Phy_Genes <- file
  Phy_concat <- data.frame(matrix(NA, nrow = dim(Phy_Genes)[1], ncol = length(concat_list)))
  
  colnames(Phy_concat) <- names(concat_list)
  rownames(Phy_concat) <- rownames(Phy_Genes)
  
  for(i in 1:length(concat_list)){
    
    Phy_concat[,i] <- apply(Phy_Genes[,concat_list[[i]]], 1, FUN = paste, collapse = "") 
   
  }
  
  return(Phy_concat)
}
