proj_wd <- getwd()
setwd("Data/Pipe_amph")

library(stringr)

#### Classification file ####

# Change the syntax of scientific names to match the phylo file
# Rownames as Scientific names

classf <- read.csv("amph_shl_new_Classification.csv")
dim(classf)
rownames(classf) <- str_replace(classf[,"Scientific.Name"], " ", "_")

write.csv(classf, "classification_tab.csv")


#### GenBank file ####

# Change the syntax of scientific names to match the phylo file
# reduce it to a matrix of gene IDs

genbank <- read.csv("amph_shl_new_GenBank.csv")
rownames(genbank) <- str_replace(genbank[,"Scientific.Name"], " ", "_")
genbank <- genbank[,-c(1,2,3,4,5)]
write.csv(genbank, "Genbank_amph.csv")


setwd(proj_wd)
