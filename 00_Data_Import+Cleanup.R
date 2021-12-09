setwd("/Users/darthvader/Desktop/Heidelberg/1_MoBi_Master/Praktikum_Bromham/Data")

install.packages("ape", "stringr",  "phylotools", "phylotools")
yes
library(phylotools)
library(stringr)

#### Alignemnt file ####

Phy <- read.phylip("amph_shl_new.phy")
Phy_names <- Phy[,1]
dim(Phy)

species_length <- dim(Phy)[1]
Seq_length <- str_length(Phy[1,2])

# Establish positions for each gene

GenBank_mat <- read.csv("amph_shl_new_GenBank.csv")
GenBank_mat$Scientific.Name <- str_replace(GenBank_mat$Scientific.Name, " ", "_")
Phy_mod <- Phy
geneloci <- c(1)
for (gene in 1:14) {
  
  NameList_NoGene <- Phy_mod[which(Phy_mod$seq.name %in% GenBank_mat[which(GenBank_mat[, 5 + gene] == ""),4]),2]
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
  Phy_mod[,2] <- substr(Phy_mod[,2], geneStart+1, str_length(Phy_mod[1,2]))
  geneloci <- append(geneloci, geneStart)
}
geneloci <- append(geneloci, Seq_length)
geneloci_mat <- matrix(ncol = 2, nrow = 15, byrow = TRUE)
geneloci_cumulative <- c(0)
for (gene in 1:15) {
  geneloci_cumulative <- geneloci_cumulative + geneloci[gene]
  geneloci_mat[gene,1] <- geneloci_cumulative
  geneloci_mat[gene,2] <- geneloci_cumulative + geneloci[gene+1] -1
  
}

#using the amph_shl_new.models file data (mit1 and nuc2 data)

Phy_Genes <- data.frame()

for (seQ in Phy[,2]) {
  
  Phy_Genes <- rbind(Phy_Genes, 
  c(
    # X12S = substr(seQ, 1, 1380),
    # X16S = substr(seQ, 1381, 3280),
    BDNF = substr(seQ, 3282, 4003),
    CXCR4 = substr(seQ, 4005, 4765),
    cytb = substr(seQ, 4766, 5905),
    H3A = substr(seQ, 5907, 6232),
    NCX1 = substr(seQ, 6234, 7567),
    ND1 = substr(seQ, 7568, 8521),
    ND2 = substr(seQ, 8522, 9493),
    POMC = substr(seQ, 9495, 10207),
    RAG1 = substr(seQ, 10209, 12637),
    RHOD = substr(seQ, 12639, 12952),
    SIA = substr(seQ, 12954, 13348),
    SLC8A3 = substr(seQ, 13350, 14488),
    TYR = substr(seQ, 14490, 15091)
  ))
  
}
colnames(Phy_Genes) <- colnames(GenBank_mat)[8:20]
row.names(Phy_Genes) <- Phy[,1]

write.csv(Phy_Genes, "Phy_Genes.csv", col.names = FALSE, row.names = TRUE)


#### Code for multiple alignement #####

system2("mkdir", args = "VertLife_MAlign")
setwd("VertLife_MAlign")

for (Gene in 1:dim(Phy_Genes)[2]){
  
  string = ""
  for (name in 1:dim(Phy_Genes)[1]) {
   
    if(str_remove_all(Phy_Genes[name, Gene], "-") != ""){
      
      sub_string <- paste("> ", rownames(Phy_Genes)[name], "\n", Phy_Genes[name, Gene], sep = "")
      string = paste(string, "\n", sub_string, sep = "")
      }
    }
  filenames <- paste("Sequences_",colnames(Phy_Genes)[Gene],sep = "" )
  cat(string, file =  filenames)
}

#### Tree File #### 

library(ape)
Tree_1000 <- read.tree("amph_shl_new_Posterior_7238.1000.trees")
class(Tree_1000)
Tree_1 <- Tree_1000[[1]]
write.tree(Tree_1)

which(!Tree_1$edge %in% Tree_2$edge)
which(!Tree_1$edge.length %in% Tree_2$edge.length)
which(!Tree_1$Nnode %in% Tree_2$Nnode)
which(!Tree_1$tip.label %in% Tree_2$tip.label)

Tree_tips <- c(rep(0,1000))
for (tree in 1:1000) {
  Tree_tips[tree] <- length(Tree_1000[[tree]][[4]])
}
write.tree(Tree_1)
library(stringr)
tip_labels <- Tree_1$tip.label

# Comparing tip names and taxonomy names
classification <- read.csv("amph_shl_new_Classification.csv")
classification_names <- classification[,"Scientific.Name"]
classification_names <- str_replace(classification_names, " ", "_")
classification[,"Scientific.Name"] <- classification_names
mismatch <- which(!tip_labels %in% classification_names)
print(mismatch)

for (tree in 2:1000) {
  
  print(setdiff(Tree_1000[[1]]$tip.label, Tree_1000[[tree]]$tip.label))
  
}

# checking whether families are monophyletic
# Establish a list of families and exclude Homo sapiens and BLANK from the list 
Family_vec <- unique(classification["Family"])
Family_vec<- Family_vec[1:(dim(Family_vec)[1]-2),]
Families_monophyletic <- c()
Families_non_monophyletic <- c()

# Iterate over every family
for (fam in Family_vec) {
  
  Family_names_vec <- classification[which(classification[,"Family"] == fam),"Scientific.Name"]
  if(is.monophyletic(Tree_1, Family_names_vec)){
    print(paste("Family ", fam, " is monophyletic.", sep = ""))
    Families_monophyletic <- append(Families_monophyletic, fam)
  }else{
    print(paste("Error: Family ", fam, " is not monophyletic.", sep = ""))
    Families_non_monophyletic <- append(Families_non_monophyletic, fam)
  }
  
}
Families_monophyletic
Families_non_monophyletic

# checking whether genuses are monophyletic
# Establish a list of genuses and exclude Homo sapiens and BLANK from the list 
Genus_vec <- unique(classification["Constraint"])
Genus_vec<- Genus_vec[1:(dim(Genus_vec)[1]-2),]
Genus_monophyletic <- c()
Genus_non_monophyletic <- c()

# Iterate over every family
for (genus in Genus_vec) {
  
  Genus_names_vec <- classification[which(classification[,"Constraint"] == genus),"Scientific.Name"]
  if(is.monophyletic(Tree_1, Genus_names_vec)){
    print(paste("Genus ", genus, " is monophyletic.", sep = ""))
    Genus_monophyletic <- append(Genus_monophyletic, genus)
  }else{
    
    print(paste("Error: Genus ", genus, " is not monophyletic.", sep = ""))
    Genus_non_monophyletic <- append(Genus_non_monophyletic, genus)
  }
  
}
Genus_monophyletic
Genus_non_monophyletic

# Dropping tips to 1 per genus
tips_list_genus <- c()
for (genus in Genus_vec) {
  
  Genus_names_vec <- classification[which(classification[,"Constraint"] == genus),"Scientific.Name"]
  tips_list_genus <- append(tips_list_genus, Genus_names_vec[1])
}

Tree_1_per_Genus <- keep.tip(Tree_1,tips_list_genus)
Tree_1_per_Genus$tip.label <- Genus_vec

# Dropping tips to 1 per family
tips_list_fam <- c()

for (fam in Family_vec) {
  
  Family_names_vec <- classification[which(classification[,"Family"] == fam),"Constraint"]
  tips_list_fam <- append(tips_list_fam, Family_names_vec[1])
}

Tree_1_per_Family <- keep.tip(Tree_1_per_Genus,tips_list_fam)
Tree_1_per_Family$tip.label <- Family_vec

# Reducing tree to tips that we have a sequence for
#remove scientific names that do not appear in the tree
Phy_names_inTree <- Phy_names[- which(!Phy_names %in% Tree_1$tip.label)]
Tree_1_phy <- keep.tip(Tree_1,Phy_names_inTree)
length(Tree_1_phy$tip.label)

#Are there trees that contain all names from the "phy" file?
for (tree in 1:1000) {
  if(length(which(!Phy_names %in% Tree_1000[[tree]]$tip.label)) == 0){
    print("All organisms in tree")
  }
}
#Apparentlòy not

#Are the species in each tree always the same?

setdiff(Tree_1000[[1]]$tip.label,
        Tree_1000[[2]]$tip.label)

# How many genes per tree tip
LociPerTip <- matrix(ncol = 15)
for (tip in Tree_1$tip.label[-1]) {
  
  LociPerTip <- rbind(LociPerTip, c(GenBank_mat[which(GenBank_mat$Scientific.Name == tip),6:20] != ""))
  
  }
LociPerTip <- LociPerTip[-1,]
colnames(LociPerTip) <- colnames(GenBank_mat)[6:20]
row.names(LociPerTip) <- Tree_1$tip.label[-1]
LociPerTip

#Determining how many species exist in each genus or family

family_table = table(classification$Family)
genus_table = table(classification$Constraint)

hist(family_table, breaks = 50, xlab = "Species per family", main = "Distribution of Species per family")
hist(genus_table, breaks = 50, xlab = "Species per genus", main = "Distribution of Species per genus")

mean(family_table)
median(family_table)
sd(family_table)
mean(genus_table)
median(genus_table)
sd(genus_table)

 # Repeat with log transformed data

family_table_ln <- log(family_table)
genus_table_ln<- log(genus_table)

hist(family_table_ln, breaks = 50, xlab = "Species per family", main = "Distribution of Species per family")
hist(genus_table_ln, breaks = 50, xlab = "Species per genus", main = "Distribution of Species per genus")


family_table_log10 <- log10(family_table)
genus_table_log10 <- log10(genus_table)

hist(family_table_log10, breaks = 50, xlab = "Species per family", main = "Distribution of Species per family")
hist(genus_table_log10, breaks = 50, xlab = "Species per genus", main = "Distribution of Species per genus")


#SHL values

Tree_SHL <- read.tree("amph_shl_new .tre")
Tree_SHL$tip.label[(!Tree_SHL$tip.label %in% Phy[,1])]

