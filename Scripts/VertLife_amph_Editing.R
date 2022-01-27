setwd("Data/VertLife_amph")

library(stringr)


#### Classification ####

classf <- read.csv("amph_shl_new_Classification.csv")
dim(classf)
rownames(classf) <- str_replace(classf[,"Scientific.Name"], " ", "_")

write.csv(classf, "classification_tab.csv")




genbank <- read.csv("amph_shl_new_GenBank.csv")
rownames(genbank) <- str_replace(genbank[,"Scientific.Name"], " ", "_")
genbank <- genbank[,-c(1,2,3,4,5)]
write.csv(genbank, "Genbank_amph.csv")
