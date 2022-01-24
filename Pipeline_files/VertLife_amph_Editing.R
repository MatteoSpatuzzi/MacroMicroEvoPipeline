setwd("Data/VertLife_sqm")

library(stringr)


#### Classification ####

classf <- read.csv("squam_shl_new_Classification.csv")
colnames(classf)[6] <- "Constraint"
colnames(classf)[7] <- "Scientific.Name"
dim(classf)
classf[,"Scientific.Name"] <- str_replace(classf[,"Scientific.Name"], " ", "_")

write.csv(classf, "classification_tab.csv")


