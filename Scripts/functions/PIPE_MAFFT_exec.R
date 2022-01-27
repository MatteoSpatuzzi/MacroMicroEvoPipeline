#'\code{MAFFT_run}
#'
#'@details Adds new sequences from fasta files to other preexisting multiple alignment fasta files with the external MAFFT software
#'
#'@param sequ Directory with the preexisting multiple alignment files
#'@param new_sequ Directory with the new sequence files to align
#'@param filewise boolean argument, if TRUE the multiple alignment is 
#'conducted on the whole files instead of sequence wise, faster but higher 
#'chance of introducing a "bad" sequence. If false the sequences are added 
#'individually and at every step the resulting alignment is compared to 
#'the previous one, if it exceeds the previous length by a given 
#'percentage the new sequence is considered a "bad" match. This criteria 
#'is very rough and can not replace a more thorough investigation of the 
#'quality of the alignemnt, thus the percentage should be always very 
#'high to only remove abvious msmatches (e.g. 30%)
#'@param reject if fileswise is FALSE, the threshold for discarding sequences,
#'default = 0.5
#'
#'@return The function will overwrite the preexisting MA files
#'
#'@author Matteo Spatuzzi, 2022
#'
#'@import phylotools
#'@import stringr
#'@import seqinr
#'



MAFFT_run <- function(sequ, new_sequ, filewise =  FALSE, reject = 0.5){
  
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
        
        mafft_cmd <- paste("\"", sequ,"/",files[[gene]], "\" > \"mafft_output\"", sep = "")
        system2('mafft', args = c('--auto', '--inputorder', "--add temp_seq.txt", mafft_cmd), stderr = "")
        
        if(str_length(read_file('mafft_output'))){
          
          MAFFT_Alignment_File <- read.fasta('mafft_output', as.string = TRUE, seqonly = TRUE)
          Old_Alignment_File <- read.fasta(paste(sequ,"/",files[[gene]], sep = ""), as.string = TRUE, seqonly = TRUE)
          
          
          if(abs(str_length(MAFFT_Alignment_File[[1]]) - str_length(Old_Alignment_File[[1]])) > reject * str_length(Old_Alignment_File[[1]])){
            
            print(paste("Gene ", gene, " sequence ", substr(New_Alignment_File_mod, 1, Start), " does not align with MA."))
            
            
          }else {
            
            #If the criteria are met, the new sequence is inserted into the Phy_Genes df
            
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
