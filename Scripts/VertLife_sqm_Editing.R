setwd("Data/Pipe_sqm")

#### Classification file ####

# Change the syntax of scientific names to match the phylo file
# Rownames as Scientific names

classf <- read.csv("squam_shl_new_Classification.csv")
colnames(classf)[7] <- "Scientific.Name"
classf <- classf[which(classf$Scientific.Name != "Homo sapiens"&classf$Scientific.Name != ""),]

dim(classf)
rownames(classf) <- str_replace(classf[,"Scientific.Name"], " ", "_")

write.csv(classf, "classification_tab.csv")

#### GenBank ####

GenBank_matrix_sqm <- read.csv("squam_shl_new_GenBank.csv")
GenBank_matrix_sqm$Species <- str_replace(GenBank_matrix_sqm$Species, " ", "_")
GenBank_matrix_sqm_mod <- GenBank_matrix_sqm[,-c(1,2,3)]
rownames(GenBank_matrix_sqm_mod) <- GenBank_matrix_sqm[,2]
GenBank_matrix_sqm_mod <- GenBank_matrix_sqm_mod[-(c((dim(GenBank_matrix_sqm_mod)[1]-3):dim(GenBank_matrix_sqm_mod)[1])),]
write.csv(GenBank_matrix_sqm_mod, "GenBank_matrix_sqm.csv")


