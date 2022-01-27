#### Functions ####

#'\code{edit_phylo}
#'
#'@details Imports a phylo multiple alignment file and if required breaks it into genes according to a matrix. Such a matrix can be generated with the find_geneloci
#' function
#'
#'@param phy An object of the phylo class or a dataframe where columns are the genes and the rows are species. 
#'In the latter case it's important that the column names correspond to the gene names if one wishes to export 
#' the sequences as fasta files for each gene.
#'@param split_mat a nx3 matrix. If the file object is a phylo or multiphylo object and one wishes to split the 
#'sequence into a separate columns, is FALSE as default.
#'@param out a string indicating the directory to export the resulting csv, is the 
#'current working directory as default.
#'
#'@return The function will export a dataframe with the sequence information split across columns based on the 
#'split_mat argument if one was provided and a new directory with a fasta file for each gene
#'
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import phylotools
#'@import stringr
#'@import seqinr
#'

edit_phylo <- function(phy, split_mat = FALSE, out = getwd(), name){
  
  library(phylotools)
  library(stringr)  
  
  Phy <- phy
  
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
  
  system2("mkdir", args = out)
  write.csv(Phy_Genes, paste(out, "/", name, sep = ""), row.names = TRUE)
  
  return(Phy_Genes)
  print("Export successful")
  
}


#'\code{dataframe_to_fasta}
#'
#'@details turns a dataframe of sequences of genes for different species into fasta files in a new directory
#'
#'@param phy A dataframe where columns are the genes and the rows are species. The column names will be used 
#'to name the fasta files.
#'@param out a string indicating the directory to export the fasta files, is the 
#'current working directory as default.
#'@param newdir the name of the new directory with the fasta files, "Phy_genes_fasta" as default
#'
#'@return The function will create a new directory and write as many fasta files as columns in the phy dataframe
#'
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import phylotools
#'@import stringr
#'@import seqinr
#'
#'
#'

dataframe_to_fasta <- function(phy, out = getwd(), newdir = "Phy_genes_fasta"){
  
  wd = getwd()
  setwd(out)
  Phy_Genes <- phy
  
  system2("mkdir", args = newdir)
  setwd(newdir)
  
  for (gene in 1:dim(Phy_Genes)[2]) {
    write.fasta(sequences = as.list(Phy_Genes[,gene]), names = (rownames(Phy_Genes)), file.out = paste("Phy_genes_", colnames(Phy_Genes)[gene], sep = ""))
  }
  setwd(wd)
}


#'\code{find_gene_loci}
#'
#'@details approximates the start and end of genes in single sequences that contain multiple gene multiple alignments.
#'
#'@param phy A phylo or multiphylo object with two columns, the first containing species names and the second multiple 
#'alignment sequences as strings
#'@param genbank_dat A dataframe or matrix with the same amount of rows as phy and as many columns as genes. Each cell 
#'must contain some string (e.g. an NCBI reference code or "1") if the gene is present in the corresponding species 
#'string in the phy object and be an empty string (="") it is not.
#'
#'@return Returns a nx3 matrix where n is the amount of columns of the genbank_mat argument; the first columns contains 
#'the gene names and the other two contain the strat and end of siad gene respectively
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import phylotools
#'@import stringr
#'

find_gene_loci <- function(phy, genbank_dat){
  
  genbank_dat <- genbank_dat
  phy_dat <- phy
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

#'\code{phy_concat}
#'
#'@details this function concatenates gene sequences from a dataframe of sequences
#'
#'@param phy dataframe with multiple alignment sequences ordered by 
#'genes (columns) and species (rows).
#'@param concat_list a list of vectors, each vector contains the names of 
#'sequences to concatenate in the phy object under the name of the 
#'corresponding list name, if a gene is not present in any vector 
#'it will NOT be present in the result
#'
#'@return The function will return a dataframe analogous to the phy 
#'argument but with columns according to the concat_list
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import ape


phy_concat <- function(phy, concat_list){
  
  library(stringr)
  
  Phy_Genes <- phy
  Phy_concat <- data.frame(matrix(NA, nrow = dim(Phy_Genes)[1], ncol = length(concat_list)))
  
  colnames(Phy_concat) <- names(concat_list)
  rownames(Phy_concat) <- rownames(Phy_Genes)
  
  for(i in 1:length(concat_list)){
    
    Phy_concat[,i] <- apply(Phy_Genes[,concat_list[[i]]], 1, FUN = paste, collapse = "") 
    
  }
  
  return(Phy_concat)
}


