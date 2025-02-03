#phylogenetic regression of limokinetic motility strategy KEGG orthologs
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
head(traitdata_kinetic)
traitdata_kinetic[,c(4:25)] <- ifelse(traitdata_kinetic[,c(4:25)] > 0,1,0)
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
phyloglm1 <- phyloglm(predictions ~ ko.K07666, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 25)
summary(phyloglm1)
#phylogenetic linear regression with Pagel's lambda
phylolm1alt <- phylolm(predictions ~ ko.K07666, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm1alt)
#Pagel's method for correlated binary trait evolution
dummy1 <- traitdata_kinetic$ko.K07666
dummypred <- traitdata_kinetic$predictions
names(dummy1) <- traitdata_kinetic$Tipname
names(dummypred) <- traitdata_kinetic$Tipname
phyloglm1alt2 <- fitPagel(tree= tree_classifier_26tips, dummy1, dummypred)
phyloglm1alt2
#
#K03782
phyloglm2 <- phyloglm(predictions ~ ko.K03782, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm2)
#phylogenetic linear regression with Pagel's lambda
phylolm2alt <- phylolm(predictions ~ ko.K03782, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm2alt)
#Pagel's method for correlated binary trait evolution
dummy2 <- traitdata_kinetic$ko.K03782
names(dummy2) <- traitdata_kinetic$Tipname
phyloglm2alt2 <- fitPagel(tree= tree_classifier_26tips, dummy2, dummypred)
phyloglm2alt2
#
#K01507
phyloglm3 <- phyloglm(predictions ~ ko.K01507, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm3)
#phylogenetic linear regression with Pagel's lambda
phylolm3alt <- phylolm(predictions ~ ko.K01507, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm3alt)
#Pagel's method for correlated binary trait evolution
dummy3 <- traitdata_kinetic$ko.K01507
names(dummy3) <- traitdata_kinetic$Tipname
phyloglm3alt2 <- fitPagel(tree= tree_classifier_26tips, dummy3, dummypred)
phyloglm3alt2
#
#K03606
phyloglm4 <- phyloglm(predictions ~ ko.K03606, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 11)
summary(phyloglm4)
#
#K16554
phyloglm5 <- phyloglm(predictions ~ ko.K16554, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm5)
#
#K01690
phyloglm6 <- phyloglm(predictions ~ ko.K01690, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm6)
#phylogenetic linear regression with Pagel's lambda
phylolm6alt <- phylolm(predictions ~ ko.K01690, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm6alt)
#Pagel's method for correlated binary trait evolution
dummy6 <- traitdata_kinetic$ko.K01690
names(dummy6) <- traitdata_kinetic$Tipname
phyloglm6alt2 <- fitPagel(tree= tree_classifier_26tips, dummy6, dummypred)
phyloglm6alt2
#
#K02302
phyloglm7 <- phyloglm(predictions ~ ko.K02302, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm7)
#phylogenetic linear regression with Pagel's lambda
phylolm7alt <- phylolm(predictions ~ ko.K02302, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm7alt)
#Pagel's method for correlated binary trait evolution
dummy7 <- traitdata_kinetic$ko.K02302
names(dummy7) <- traitdata_kinetic$Tipname
phyloglm7alt2 <- fitPagel(tree= tree_classifier_26tips, dummy7, dummypred)
phyloglm7alt2
#
#K04047
phyloglm8 <- phyloglm(predictions ~ ko.K04047, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm8)
#phylogenetic linear regression with Pagel's lambda
phylolm8alt <- phylolm(predictions ~ ko.K04047, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm8alt)
#Pagel's method for correlated binary trait evolution
dummy8 <- traitdata_kinetic$ko.K04047
names(dummy8) <- traitdata_kinetic$Tipname
phyloglm8alt2 <- fitPagel(tree= tree_classifier_26tips, dummy8, dummypred)
phyloglm8alt2
#
#K07645
phyloglm9 <- phyloglm(predictions ~ ko.K07645, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm9)
#phylogenetic linear regression with Pagel's lambda
phylolm9alt <- phylolm(predictions ~ ko.K07645, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm9alt)
#Pagel's method for correlated binary trait evolution
dummy9 <- traitdata_kinetic$ko.K07645
names(dummy9) <- traitdata_kinetic$Tipname
phyloglm9alt2 <- fitPagel(tree= tree_classifier_26tips, dummy9, dummypred)
phyloglm9alt2
#
#K20977
phyloglm10 <- phyloglm(predictions ~ ko.K20977, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm10)
#phylogenetic linear regression with Pagel's lambda
phylolm10alt <- phylolm(predictions ~ ko.K20977, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm10alt)
#Pagel's method for correlated binary trait evolution
dummy10 <- traitdata_kinetic$ko.K20977
names(dummy10) <- traitdata_kinetic$Tipname
phyloglm10alt2 <- fitPagel(tree= tree_classifier_26tips, dummy10, dummypred)
phyloglm10alt2
#
#K20978
phyloglm11 <- phyloglm(predictions ~ ko.K20978, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm11)
#phylogenetic linear regression with Pagel's lambda
phylolm11alt <- phylolm(predictions ~ ko.K20978, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm11alt)
#Pagel's method for correlated binary trait evolution
dummy11 <- traitdata_kinetic$ko.K20978
names(dummy11) <- traitdata_kinetic$Tipname
phyloglm11alt2 <- fitPagel(tree= tree_classifier_26tips, dummy11, dummypred)
phyloglm11alt2
#
#K03411
phyloglm12 <- phyloglm(predictions ~ ko.K03411, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm12)
#phylogenetic linear regression with Pagel's lambda
phylolm12alt <- phylolm(predictions ~ ko.K03411, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm12alt)
#Pagel's method for correlated binary trait evolution
dummy12 <- traitdata_kinetic$ko.K03411
names(dummy12) <- traitdata_kinetic$Tipname
phyloglm12alt2 <- fitPagel(tree= tree_classifier_26tips, dummy12, dummypred)
phyloglm12alt2
#
#K20920
phyloglm13 <- phyloglm(predictions ~ ko.K20920, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm13)
#
#K20988
phyloglm14 <- phyloglm(predictions ~ ko.K20988, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm14)
#
#K09797
phyloglm15 <- phyloglm(predictions ~ ko.K09797, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm15)
#phylogenetic linear regression with Pagel's lambda
phylolm15alt <- phylolm(predictions ~ ko.K09797, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm15alt)
#Pagel's method for correlated binary trait evolution
dummy15 <- traitdata_kinetic$ko.K09797
names(dummy15) <- traitdata_kinetic$Tipname
phyloglm15alt2 <- fitPagel(tree= tree_classifier_26tips, dummy15, dummypred)
phyloglm15alt2
#
#K00140
phyloglm16 <- phyloglm(predictions ~ ko.K00140, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm16)
#phylogenetic linear regression with Pagel's lambda
phylolm16alt <- phylolm(predictions ~ ko.K00140, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm16alt)
#Pagel's method for correlated binary trait evolution
dummy16 <- traitdata_kinetic$ko.K00140
names(dummy16) <- traitdata_kinetic$Tipname
phyloglm16alt2 <- fitPagel(tree= tree_classifier_26tips, dummy16, dummypred)
phyloglm16alt2
#
#K00523
phyloglm17 <- phyloglm(predictions ~ ko.K00523, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm17)
#phylogenetic linear regression with Pagel's lambda
phylolm17alt <- phylolm(predictions ~ ko.K00523, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm17alt)
#Pagel's method for correlated binary trait evolution
dummy17 <- traitdata_kinetic$ko.K00523
names(dummy17) <- traitdata_kinetic$Tipname
phyloglm17alt2 <- fitPagel(tree= tree_classifier_26tips, dummy17, dummypred)
phyloglm17alt2
#
#K00641
phyloglm18 <- phyloglm(predictions ~ ko.K00641, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm18)
#phylogenetic linear regression with Pagel's lambda
phylolm18alt <- phylolm(predictions ~ ko.K00641, data= traitdata_kinetic, phy= tree_classifier_26tips, model=c('lambda'))
summary(phylolm18alt)
#Pagel's method for correlated binary trait evolution
dummy18 <- traitdata_kinetic$ko.K00641
names(dummy18) <- traitdata_kinetic$Tipname
phyloglm18alt2 <- fitPagel(tree= tree_classifier_26tips, dummy18, dummypred)
phyloglm18alt2
#
#K03410
phyloglm19 <- phyloglm(predictions ~ ko.K03410, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol =19)
summary(phyloglm19)
#phylogenetic linear regression with Pagel's lambda
phylolm19alt <- phylolm(predictions ~ ko.K03410, data= traitdata_kinetic, phy= tree_classifier_26tips, method='lambda')
summary(phylolm19alt)
#Pagel's method for correlated binary trait evolution
dummy19 <- traitdata_kinetic$ko.K03410
names(dummy19) <- traitdata_kinetic$Tipname
phyloglm19alt2 <- fitPagel(tree= tree_classifier_26tips, dummy19, dummypred)
phyloglm19alt2
#
#K12972
phyloglm20 <- phyloglm(predictions ~ ko.K12972, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 20)
summary(phyloglm20)
#phylogenetic linear regression with Pagel's lambda
phylolm20alt <- phylolm(predictions ~ ko.K12972, data= traitdata_kinetic, phy= tree_classifier_26tips, method='lambda')
summary(phylolm20alt)
#Pagel's method for correlated binary trait evolution
dummy20 <- traitdata_kinetic$ko.K12972
names(dummy20) <- traitdata_kinetic$Tipname
phyloglm20alt2 <- fitPagel(tree= tree_classifier_26tips, dummy20, dummypred)
phyloglm20alt2
#
#K02106
phyloglm21 <- phyloglm(predictions ~ ko.K02106, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 20)
summary(phyloglm21)
#phylogenetic linear regression with Pagel's lambda
phylolm21alt <- phylolm(predictions ~ ko.K02106, data= traitdata_kinetic, phy= tree_classifier_26tips, method='lambda')
summary(phylolm21alt)
#Pagel's method for correlated binary trait evolution
dummy21 <- traitdata_kinetic$ko.K02106
names(dummy21) <- traitdata_kinetic$Tipname
phyloglm21alt2 <- fitPagel(tree= tree_classifier_26tips, dummy21, dummypred)
phyloglm21alt2
#
#K03521
phyloglm22 <- phyloglm(predictions ~ ko.K03521, data= traitdata_kinetic, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm22)
#phylogenetic linear regression with Pagel's lambda
phylolm22alt <- phylolm(predictions ~ ko.K03521, data= traitdata_kinetic, phy= tree_classifier_26tips, method='lambda')
summary(phylolm22alt)
#Pagel's method for correlated binary trait evolution
dummy22 <- traitdata_kinetic$ko.K03521
names(dummy22) <- traitdata_kinetic$Tipname
phyloglm22alt2 <- fitPagel(tree= tree_classifier_26tips, dummy22, dummypred)
phyloglm22alt2
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
output_phylinreg <- rbind(summary(phylolm2alt)$coefficients, summary(phylolm3alt)$coefficients, summary(phylolm6alt)$coefficients, summary(phylolm7alt)$coefficients, summary(phylolm8alt)$coefficients, summary(phylolm9alt)$coefficients, summary(phylolm10alt)$coefficients, summary(phylolm11alt)$coefficients, summary(phylolm12alt)$coefficients, summary(phylolm15alt)$coefficients, summary(phylolm16alt)$coefficients, summary(phylolm17alt)$coefficients, summary(phylolm18alt)$coefficients, summary(phylolm19alt)$coefficients, summary(phylolm20alt)$coefficients, summary(phylolm21alt)$coefficients, summary(phylolm22alt)$coefficients)
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
#write alternative phylogenetic method Pagel's correlated binary trait evolution summary data to file
output_PagelBinaryTraits <- data.frame(KO = c("K07666", "K03782", "K01507", "K01690", "K02302", "K04047", "K07645", "K20977", "K20978", "K03411", "K09797", "K00140", "K00523", "K00641", "K03410", "K12972", "K02106", "K03521"), pvalue = c(phyloglm1alt2$P[1], phyloglm2alt2$P[1], phyloglm3alt2$P[1], phyloglm6alt2$P[1], phyloglm7alt2$P[1], phyloglm8alt2$P[1], phyloglm9alt2$P[1], phyloglm10alt2$P[1], phyloglm11alt2$P[1], phyloglm12alt2$P[1], phyloglm15alt2$P[1], phyloglm16alt2$P[1], phyloglm17alt2$P[1], phyloglm18alt2$P[1], phyloglm19alt2$P[1], phyloglm20alt2$P[1], phyloglm21alt2$P[1], phyloglm22alt2$P[1]))
output_PagelBinaryTraits
write.csv(output_PagelBinaryTraits, file = "output_PagelBinary_limokinetic_KO.csv")