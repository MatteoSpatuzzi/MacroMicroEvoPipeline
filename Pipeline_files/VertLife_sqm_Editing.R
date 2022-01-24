WorkingDirectory <- "/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data/Vertlife_sqm"
setwd(WorkingDirectory)

#### Classification ####

classf <- read.csv("squam_shl_new_Classification.csv")
colnames(classf)[6] <- "Constraint"
colnames(classf)[7] <- "Scientific.Name"
classf <- classf[which(classf$Scientific.Name != "Homo sapiens"&classf$Scientific.Name != ""),]

dim(classf)
classf[,"Scientific.Name"] <- str_replace(classf[,"Scientific.Name"], " ", "_")

write.csv(classf, "classification_tab.csv")

#### GenBank ####

GenBank_matrix_sqm <- read.csv("squam_shl_new_GenBank.csv")
GenBank_matrix_sqm$Species <- str_replace(GenBank_matrix_sqm$Species, " ", "_")
GenBank_matrix_sqm_mod <- GenBank_matrix_sqm[,-c(1,2,3)]
rownames(GenBank_matrix_sqm_mod) <- GenBank_matrix_sqm[,2]
write.csv(GenBank_matrix_sqm_mod, "GenBank_matrix_sqm.csv")


