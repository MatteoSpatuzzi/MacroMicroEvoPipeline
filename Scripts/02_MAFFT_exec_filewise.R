library(readr)
library(stringr)
#Working directory and other important global variables are set, data is imported
WorkingDirectory <- "/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/"
setwd(WorkingDirectory)

Gene_names <- c('12S', '16S', 'BDNF', 'CXCR4', 'cytb', 'H3A', 'NCX1',
                'ND1', 'ND2', 'POMC', 'RAG1', 'RHOD', 'SIA', 'SLC8A3', 'TYR')

Phy_Genes <- read.csv("Phy_Genes.csv")
rownames(Phy_Genes) <- Phy_Genes[,1]
Phy_Genes <- Phy_Genes[,-1]
colnames(Phy_Genes)[c(1,2)] <- c("12S", "16S")

#Set modus operandi
# filewise TRUE: runs multiple alignment with entire file, much faster but does not remove sequences that align poorly, outputs a set of multiple alignemnt files
#filewise FALSE: runs multplie alignemnt individualy with each sequence, can remove poor alignments and already adds them to the Phy_genes dataframe, returns Phy_genes as csv

filewise = TRUE


if(filewise){
  
  #Iterate over every gene
  for( gene in Gene_names){
    
    #Import files from respective directories
    Alignment_File <- read_file(paste("VertLife_MAlign/Sequences_", gene, sep = ""))
    New_Alignment_File <- read_file(paste("New_MAlign/New_alignment_sequences_", gene, sep = ""))
    New_Alignment_File_mod <- New_Alignment_File
    
    #Iterate over every new sequence for the gene
    system2("mkdir", args = "MAFFT_output")
    system2("mafft", args=c("--thread -1", "--add", paste("New_MAlign/New_alignment_sequences_",gene, sep = ""),paste("VertLife_MAlign/Sequences_",gene," > MAFFT_output/",gene,"_mafft", sep = "")))
    
  }
  
}else{
  #Iterate over every gene
  for (gene in Gene_names){
    
    #Import files from respective directories
    Alignment_File <- read_file(paste("VertLife_MAlign/Sequences_", gene, sep = ""))
    New_Alignment_File <- read_file(paste("New_MAlign/New_alignment_sequences_", gene, sep = ""))
    New_Alignment_File_mod <- New_Alignment_File
    
    #Iterate over every new sequence for the gene
    while(str_length(New_Alignment_File_mod) > 0){
      
      #Isolate the first seqeunce with string manipulation methods
      Start <- str_locate(New_Alignment_File_mod, pattern = ">")[1]
      Stop <- min(str_locate_all(New_Alignment_File_mod, pattern = "\n")[[1]][min(2,str_count(New_Alignment_File_mod, "\n")),1] -2, str_length(New_Alignment_File_mod))
      Seq <- substr(New_Alignment_File_mod, Start, Stop)
      
      
      #Write temporary files in the wd to run the command line
      #Subsequently, run MAFFT from the command line 
      
      writeLines(str_to_lower(Seq), "new_Seq.txt")
      system2('mafft', args = c('--auto', '--inputorder', "--add new_Seq.txt", "VertLife_MAlign/Sequences_12S > mafft_output"), stderr = "")
      MAFFT_Alignment_File <- read_file('mafft_output')
      
      #To establish rapidly whether the alignment has been successful, the difference is length of the multiple alignemnt before and after the new seqeunce is calculated
      #A small or even negative difference is used as an indicator of a successful alignment 
      #This criteria should be revised!
      Seq_length <- str_length(substr(Seq,
                                      start = str_locate(Seq, "\n")[1,1] + 2,
                                      stop = str_length(Seq)
      )
      )
      
      MAFFT_Alignment_length <- str_length(str_remove_all(substr(MAFFT_Alignment_File,
                                                                 start = str_locate(MAFFT_Alignment_File, "\n") +2, 
                                                                 stop = str_locate_all(MAFFT_Alignment_File, ">")[[1]][2,1]-3
      ), pattern = "\n"
      )
      )
      
      Alignment_length <- str_length(str_remove_all(substr(Alignment_File,
                                                           start = str_locate_all(Alignment_File, "\n")[[1]][2,1] +2, 
                                                           stop = str_locate_all(Alignment_File, ">")[[1]][2,1]-3
      ), pattern = "\n"
      )
      )
      
      if((Seq_length - (MAFFT_Alignment_length - Alignment_length)) < 0.5 * Seq_length){
        
        print(paste("Gene ", gene, " sequence ", substr(New_Alignment_File_mod, 1, Start), " does not align with MA."))
        
        
      }else {
        
        #If the critera are met, the new seqeunce is inserted into the Phy_Genes df
        name <- substr(Seq, str_locate(Seq, ">")+2, str_locate(Seq, "\n")-3)
        seq <- substr(Seq, str_locate(Seq, "\n")+1, str_length(Seq))
        Phy_Genes[name,gene] <- seq
        print("substituted:")
        print(seq)
      }
      
      New_Alignment_File_mod <- substr(New_Alignment_File_mod, Stop + 3, str_length(New_Alignment_File_mod))
      
    }
    
  }
  
  write.csv(Phy_Genes, file = "Phy_Genes_update.csv")
}






write.csv(Phy_Genes, file = "Phy_Genes_update.csv")

