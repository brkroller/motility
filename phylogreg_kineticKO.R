#phylogenetic logistic regression for orthologs in limokinetic motility strategy classifier
#
#load packages
require(phylolm)
require(phytools)
#
#load phylogenetic tree
tree_full <- read.tree("dataset_fulltree.treefile",tree.names=T)
tree_full
#
#
#load limokinetic data
traitdata_kinetic <- read.csv('dataset_limokinetic_KO.csv', header = T)
#
#convert traitdata_kinetic kegg data to binary variables
dim(traitdata_kinetic)
traitdata_kinetic[,c(3:24)] <- ifelse(traitdata_kinetic[,c(3:24)] > 0,1,0)
head(traitdata_kinetic)
#
#
#data prep for phylogenetic regression
#keep only tips with trait data
tree_classifier_26tips <- keep.tip(tree_full, tip = c(traitdata_kinetic$Tipname))
tree_classifier_26tips
#
#plot trimmed tree for visual inspection
plot(tree_classifier_26tips) 
#
#tree tip labels as rownames
row.names(traitdata_kinetic) <- traitdata_kinetic$Tipname
#
#visually inspect data
head(traitdata_kinetic)
#
#
#phylogenetic logistic regression limokinetic classifier outcome vs. all features
#
#K07666
phyloglm1 <- phyloglm(class ~ ko.K07666, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 25)
summary(phyloglm1)
#phylogenetic linear regression with Pagel's lambda
phylolm1alt <- phylolm(class ~ ko.K07666, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm1alt)
#linear regression without considering phylogeny
lm1 <- lm(class ~ ko.K07666, data= traitdata_kinetic)
summary(lm1)
#Pagel's method for correlated binary trait evolution
dummy1 <- traitdata_kinetic$ko.K07666
dummy2 <- traitdata_kinetic$class
names(dummy1) <- traitdata_kinetic$Tipname
names(dummy2) <- traitdata_kinetic$Tipname
phyloglm1alt2 <- fitPagel(tree= tree_classifier_26tips, dummy1, dummy2)
phyloglm1alt2
#
#K03782
phyloglm2 <- phyloglm(class ~ ko.K03782, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm2)
#
#K01507
phyloglm3 <- phyloglm(class ~ ko.K01507, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm3)
#
#K03606
phyloglm4 <- phyloglm(class ~ ko.K03606, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 11)
summary(phyloglm4)
#
#K16554
phyloglm5 <- phyloglm(class ~ ko.K16554, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm5)
#
#K01690
phyloglm6 <- phyloglm(class ~ ko.K01690, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm6)
#
#K02302
phyloglm7 <- phyloglm(class ~ ko.K02302, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm7)
#
#K04047
phyloglm8 <- phyloglm(class ~ ko.K04047, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm8)
#
#K07645
phyloglm9 <- phyloglm(class ~ ko.K07645, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm9)
#
#K20977
phyloglm10 <- phyloglm(class ~ ko.K20977, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm10)
#
#K20978
phyloglm11 <- phyloglm(class ~ ko.K20978, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm11)
#
#K03411
phyloglm12 <- phyloglm(class ~ ko.K03411, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm12)
#
#K20920
phyloglm13 <- phyloglm(class ~ ko.K20920, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm13)
#
#K20988
phyloglm14 <- phyloglm(class ~ ko.K20988, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm14)
#
#K09797
phyloglm15 <- phyloglm(class ~ ko.K09797, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm15)
#
#K00140
phyloglm16 <- phyloglm(class ~ ko.K00140, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm16)
#
#K00523
phyloglm17 <- phyloglm(class ~ ko.K00523, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm17)
#
#K00641
phyloglm18 <- phyloglm(class ~ ko.K00641, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm18)
#
#K03410
phyloglm19 <- phyloglm(class ~ ko.K03410, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol =19)
summary(phyloglm19)
#phylogenetic linear regression with Pagel's lambda
phylolm19alt <- phylolm(class ~ ko.K03410, data= traitdata_kinetic, phy= tree_classifier_26tips, method='lambda')
summary(phylolm19alt)
#linear regression without considering phylogeny
lm19 <- lm(class ~ ko.K03410, data= traitdata_kinetic)
summary(lm19)
#Pagel's method for correlated binary trait evolution
dummy3 <- traitdata_kinetic$ko.K03410
names(dummy3) <- traitdata_kinetic$Tipname
phyloglm19alt2 <- fitPagel(tree= tree_classifier_26tips, dummy3, dummy2)
phyloglm19alt2
#
#K12972
phyloglm20 <- phyloglm(class ~ ko.K12972, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 20)
summary(phyloglm20)
#phylogenetic linear regression with Pagel's lambda
phylolm20alt <- phylolm(class ~ ko.K12972, data= traitdata_kinetic, phy= tree_classifier_26tips, method='lambda')
summary(phylolm20alt)
#linear regression without considering phylogeny
lm20 <- lm(class ~ ko.K12972, data= traitdata_kinetic)
summary(lm20)
#Pagel's method for correlated binary trait evolution
dummy4 <- traitdata_kinetic$ko.K12972
names(dummy4) <- traitdata_kinetic$Tipname
phyloglm20alt2 <- fitPagel(tree= tree_classifier_26tips, dummy4, dummy2)
phyloglm20alt2
#
#K02106
phyloglm21 <- phyloglm(class ~ ko.K02106, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 20)
summary(phyloglm21)
#phylogenetic linear regression with Pagel's lambda
phylolm21alt <- phylolm(class ~ ko.K02106, data= traitdata_kinetic, phy= tree_classifier_26tips, method='lambda')
summary(phylolm21alt)
#linear regression without considering phylogeny
lm21 <- lm(class ~ ko.K02106, data= traitdata_kinetic)
summary(lm21)
#Pagel's method for correlated binary trait evolution
dummy5 <- traitdata_kinetic$ko.K02106
names(dummy5) <- traitdata_kinetic$Tipname
phyloglm21alt2 <- fitPagel(tree= tree_classifier_26tips, dummy5, dummy2)
phyloglm21alt2
#
#K03521
phyloglm22 <- phyloglm(class ~ ko.K03521, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm22)
#
#write phylogreg summary data to file
output_phylogreg <- rbind(summary(phyloglm1)$coefficients, summary(phyloglm2)$coefficients, summary(phyloglm3)$coefficients, summary(phyloglm4)$coefficients, summary(phyloglm5)$coefficients, summary(phyloglm6)$coefficients, summary(phyloglm7)$coefficients, summary(phyloglm8)$coefficients, summary(phyloglm9)$coefficients, summary(phyloglm10)$coefficients, summary(phyloglm11)$coefficients, summary(phyloglm12)$coefficients, summary(phyloglm13)$coefficients, summary(phyloglm14)$coefficients, summary(phyloglm15)$coefficients, summary(phyloglm16)$coefficients, summary(phyloglm17)$coefficients, summary(phyloglm18)$coefficients, summary(phyloglm19)$coefficients, summary(phyloglm20)$coefficients, summary(phyloglm21)$coefficients, summary(phyloglm22)$coefficients)
#make new column to signify the term (slope or intercept) and the KO (always paired sets of the KO indicated by the slope term)
term <- rownames(output_phylogreg)
#convert output to a data frame
output_phylogreg <- data.frame(output_phylogreg)
#add the new column term to the output dataframe
output_phylogreg$term <- term
#change order of columns in the dataset
output_phylogreg <- output_phylogreg[,c(5,1:4)]
#inspect dataset
output_phylogreg
write.csv(output_phylogreg, row.names=F, file = "output_PhyLogReg_limokinetic_KO.csv")
#
#write alternative phylogenetic linear regression summary data to file
output_phylinreg <- rbind(summary(phylolm1alt)$coefficients, summary(phylolm19alt)$coefficients, summary(phylolm20alt)$coefficients, summary(phylolm21alt)$coefficients)
output_phylinreg
term_phylinreg <- rownames(output_phylinreg)
#convert output to a data frame
output_phylinreg <- data.frame(output_phylinreg)
#add the new column term to the output dataframe
output_phylinreg$term <- term_phylinreg
#change order of columns in the dataset
output_phylinreg <- output_phylinreg[,c(5,1:4)]
#inspect dataset
output_phylinreg
write.csv(output_phylinreg, row.names=F, file = "output_PhyLinearReg_limokinetic_KO.csv")
#
#write alternative linear regression summary data to file
output_linreg <- rbind(summary(lm1)$coefficients, summary(lm19)$coefficients, summary(lm20)$coefficients, summary(lm21)$coefficients)
output_linreg
term_linreg <- rownames(output_linreg)
#convert output to a data frame
output_linreg <- data.frame(output_linreg)
#add the new column term to the output dataframe
output_linreg$term <- term_linreg
#change order of columns in the dataset
output_linreg <- output_linreg[,c(5,1:4)]
#inspect dataset
output_linreg
#adjust column names
colnames(output_linreg) <- c("term", "Estimate","StdErr","t.value","p.value")
#inspect dataset
output_linreg
write.csv(output_linreg, row.names=F, file = "output_LinearReg_limokinetic_KO.csv")
#
#write alternative phylogenetic method Pagel's correlated binary trait evolution summary data to file
output_PagelBinaryTraits <- data.frame(KO = c("K07666", "K03410", "K12972", "K02106"), pvalue = c(phyloglm1alt2$P[1], phyloglm19alt2$P[1], phyloglm20alt2$P[1], phyloglm21alt2$P[1]))
output_PagelBinaryTraits
write.csv(output_PagelBinaryTraits, file = "output_PagelBinary_limokinetic_KO.csv")
