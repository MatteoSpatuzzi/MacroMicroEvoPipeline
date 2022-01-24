MAFFT_run <- function(sequ, new_sequ, filewise =  FALSE){
  
  wd <- getwd()
  
  #Import old seq
  setwd(sequ)
  files <- system2("ls",stdout = TRUE)
  setwd(wd)
  
  setwd(new_sequ)
  new_files <- system2("ls",stdout = TRUE)
  new_fasta_list  <- lapply(new_files, read.fasta, as.string = TRUE)
  
  setwd(wd)
  
  system2("mkdir", args = "new_MAlign")

  for(gene in 1:length(files)) {
    
    if(filewise){
      
      system2('mafft', args = c('--auto', '--inputorder', paste("--add ", new_sequ, "/", new_files[[gene]], sep = ""), paste("\"", sequ,"/",files[[gene]], "\" > \"new_MAlign/",files[[gene]],"\"", sep = "")), stderr = "")
      
      
    }else{
      
      for(seq in 1:length(new_fasta_list[[gene]])) {
        
        setwd(wd)
        write.fasta(new_fasta_list[[1]][[1]], names = str_remove(str_remove(attr(new_fasta_list[[1]][[1]], "Annot"), " "),">"), file.out = "temp_seq.txt")
        system2('mafft', args = c('--auto', '--inputorder', "--add temp_seq.txt", paste("\"", sequ,"/",files[[gene]], "\" > \"mafft_output\"", sep = "")), stderr = "")
        if(str_length(read_file('mafft_output'))){
          
          MAFFT_Alignment_File <- read.fasta('mafft_output', as.string = TRUE, seqonly = TRUE)
          Old_Alignment_File <- read.fasta(paste(sequ,"/",files[[gene]], sep = ""), as.string = TRUE, seqonly = TRUE)
          
          
          if(abs(str_length(MAFFT_Alignment_File[[1]]) - str_length(Old_Alignment_File[[1]])) > 0.5 * str_length(Old_Alignment_File[[1]])){
            
            print(paste("Gene ", gene, " sequence ", substr(New_Alignment_File_mod, 1, Start), " does not align with MA."))
            
            
          }else {
            
            #If the critera are met, the new seqeunce is inserted into the Phy_Genes df
            
            setwd(wd)
            setwd(sequ)
            
            MAFFT_Alignment_File <- read.fasta('mafft_output', as.string = TRUE)
            setwd(sequ)
            write.fasta(MAFFT_Alignment_File, names = names(MAFFT_Alignment_File),  file.out =  files[gene])
            
          }
        }
      }
    }
    
  }
  
}
