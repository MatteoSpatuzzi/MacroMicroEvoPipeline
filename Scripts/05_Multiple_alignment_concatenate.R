library(stringr)

WorkingDirectory <- "/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data"
setwd(WorkingDirectory)
Phy_Genes <- read.csv("Phy_genes_reduced")[,-c(2,3)]

Phy_concat <- function(x){
  
  Phy_concat <- data.frame(matrix(NA, nrow = dim(x)[1], ncol = 3))
  colnames(Phy_concat) <- c("mit", "nuc", "rib")
  rownames(Phy_concat) <- x$X
  Phy_concat$mit <- apply(x[,c(4,7,8)], 1, FUN = paste, collapse = "") 
  Phy_concat$nuc <- apply(x[,-c(1,2,3,4,7,8)], 1, FUN = paste, collapse = "")
  Phy_concat$rib <- apply(x[,c(2,3)], 1, FUN = paste, collapse = "")
  
  return(Phy_concat)
}
Phy_concat(Phy_Genes)

write.csv(Phy_concat, "Multiple_Alignment_concat")
