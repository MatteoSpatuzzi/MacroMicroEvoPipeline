#### Functions ####

# This function imports a phylo multiple alignemtn file and if required breaks it into genes according
# to a matrix
#
# file must be a phylo or multiphylo object
# split_mat must be a matrix with n rows and three columns, each row correspont to a gene and contains (in 
# this order) the name, start and end index of the gene in the multiple alignemnt string
# E.g. gene_mat <- matrix(c(
# "BDNF", 3282, 4003,
# "CXCR4",  4005, 4765,
# "cytb", 4766, 5905,
# "H3A", 5907, 6232,
# "NCX1", 6234, 7567,
# "ND1", 7568, 8521,
# "ND2", 8522, 9493,
# "POMC", 9495, 10207,
# "RAG1", 10209, 12637,
# "RHOD", 12639, 12952,
# "SIA", 12954, 13348,
# "SLC8A3", 13350, 14488,
# "TYR", 14490, 15091), ncol = 3, byrow = TRUE)

Edit_phylo <- function(file, split_mat = FALSE, out = getwd()){
  
  library(phylotools)
  library(stringr)  
  
  Phy <- file
  
  print(paste("Object of class", class(Phy), "with", dim(Phy)[1], "rows and", dim(Phy)[2], "columns.", sep = " "))
  
  if(all(class(split_mat) == c("matrix", "array"))){
    
    Phy_Genes <- c()
    for(i in 1:dim(split_mat)[1]){
      
      split_list<- lapply(Phy[,2], function(x) substr(x, as.integer(split_mat[i,2]),  as.integer(split_mat[i,3])))
      split_vec <- unlist(split_list)
      Phy_Genes <- cbind(Phy_Genes, split_vec)
     }
    
    colnames(Phy_Genes) <- split_mat[,1]
    rownames(Phy_Genes) <- Phy$seq.name
  
    }else{
    
    Phy_Genes <- Phy
    
  }
  
  Phy_Genes <- Phy_Genes[which(apply(Phy_Genes, 1, function(x) any(!str_remove_all(x, "-") == ""))),]
  
  write.csv(Phy_Genes, paste(out, "/", "Phy_Genes.csv", sep = ""), row.names = TRUE)
  
  system2("mkdir", args = "phy_genes_fasta")
  setwd("phy_genes_fasta")
  
  for (gene in 1:dim(Phy_Genes)[2]) {
    write.fasta(sequences = as.list(Phy_Genes[,gene]), names = (rownames(Phy_Genes)), file.out = paste("Phy_genes_", colnames(Phy_Genes)[gene], sep = ""))
  }
  
  return(Phy_Genes)
  print("Export successful")
}


Find_gene_loci <- function(phy_dat, genbank_dat){
  
  genbank_dat <- genbank_dat
  phy_dat <- phy_dat
  ngenes <- dim(genbank_dat)[2]
  geneloci <- c(1)
  
  for (gene in 1:(ngenes-1)) {
    
    NameList_NoGene <- phy_dat[which(phy_dat$seq.name %in% rownames(genbank_dat)[which(genbank_dat[,gene] == "")]),2]
    geneStart <- 200000
    
    for (species in 1:length(NameList_NoGene)) {
      geneStart <- 
        min(str_locate(NameList_NoGene[species], "A")[1],
            str_locate(NameList_NoGene[species], "T")[1],
            str_locate(NameList_NoGene[species], "G")[1],
            str_locate(NameList_NoGene[species], "C")[1],
            geneStart, na.rm = TRUE
        )
    }
    
    print(geneStart)
    phy_dat[,2] <- substr(phy_dat[,2], geneStart+1, str_length(phy_dat[1,2]))
    geneloci <- append(geneloci, geneStart)
  }
  
  geneloci <- append(geneloci, str_length(phy_dat[1,2]))
  geneloci_mat <- matrix(ncol = 2, nrow = ngenes, byrow = TRUE)
  geneloci_cumulative <- c(0)
  for (gene in 1:ngenes) {
    geneloci_cumulative <- geneloci_cumulative + geneloci[gene]
    geneloci_mat[gene,1] <- geneloci_cumulative
    geneloci_mat[gene,2] <- geneloci_cumulative + geneloci[gene+1] -1
    
  }
  
  return(geneloci_mat)
  
}



