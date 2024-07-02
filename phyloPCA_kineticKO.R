#Phylogenetic PCA of limokinetic dataset
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
#data prep for phylogenetic PCA
#keep only tips with trait data
tree_classifier_26tips <- keep.tip(tree_full, tip = c(traitdata_kinetic$Tipname))
tree_classifier_26tips
#
#plot trimmed tree for visual inspection
plot(tree_classifier_26tips) 
#
#requires tree tip labels as rownames
row.names(traitdata_kinetic) <- traitdata_kinetic$Tipname
#
#visually inspect data
head(traitdata_kinetic)
#
#
#phylogenetic PCA of limokinetic data
#trim traitdata dataset to remove first 2 columns
traitdata_kinetic_ppca <- traitdata_kinetic[,-c(1,2)]
pca.limo <- phyl.pca(tree_classifier_26tips, traitdata_kinetic_ppca, method = "BM", mode = "cov")
#
#save the data points in PCA axis space
pca.limo.points <- pca.limo$S
#extract the real number components of the complex numbers in points
pca.limo.points.real <- Re(pca.limo.points)
pca.limo.load <- pca.limo$L
#extract the real number components of the complex numbers in loadings
pca.limo.load.real <- Re(pca.limo.load)
pca.limo.eigvec <- pca.limo$Evec
pca.limo.eigval <- diag(pca.limo$Eval)
dim(pca.limo.points.real)
pca.limo.eigval.real <- Re(pca.limo.eigval)
pca.limo.eigval.real.pos <- subset(pca.limo.eigval.real, pca.limo.eigval.real > 0)
#confirm order of factors dataframe matches pca points dataframe
head(pca.limo.points[,1:5])
head(traitdata_kinetic_ppca)
tail(pca.limo.points[,1:5])
tail(traitdata_kinetic_ppca)
#variance explained by PC axes
pc1 <- pca.limo.eigval.real.pos[1]/sum(pca.limo.eigval.real.pos)
pc2 <- pca.limo.eigval.real.pos[2]/sum(pca.limo.eigval.real.pos)
pc3 <- pca.limo.eigval.real.pos[3]/sum(pca.limo.eigval.real.pos)
pc4 <- pca.limo.eigval.real.pos[4]/sum(pca.limo.eigval.real.pos)
pc5 <- pca.limo.eigval.real.pos[5]/sum(pca.limo.eigval.real.pos)
pc6 <- pca.limo.eigval.real.pos[6]/sum(pca.limo.eigval.real.pos)
pc7 <- pca.limo.eigval.real.pos[7]/sum(pca.limo.eigval.real.pos)
pc8 <- pca.limo.eigval.real.pos[8]/sum(pca.limo.eigval.real.pos)
pc9 <- pca.limo.eigval.real.pos[9]/sum(pca.limo.eigval.real.pos)
pc10 <- pca.limo.eigval.real.pos[10]/sum(pca.limo.eigval.real.pos)
pc1
pc2
pc3
pc4
pc5
pc6
pc7
pc8
pc9
pc10
#variance explained by sum of first six pPCA axes
sum(pca.limo.eigval.real.pos[1:6])/sum(pca.limo.eigval.real.pos)
#ordered loadings of PC1
#save all loadings on first pPCA axes
pca.limo.load.pc1 <- pca.limo.load.real[,1]
#save the top 100 loadings for the first 6 pPCA axes
#first make a separate dataframe for the first 6 pPCA axes that is ordered
pca.limo.load.pc1.ord <- as.data.frame(pca.limo.load.pc1[order(-abs(pca.limo.load.pc1))])
pca.limo.load.pc1.ord
#
#Regress pPCA axis (trait data) vs. classifier outcome
#make a data frame with each ppca axis and the trait data, since order of data frames was already confirmed above
pca.limo.points.class <- data.frame(pca.limo.points)
pca.limo.points.class$class <- traitdata_kinetic$class
#phylogenetic logistic regression limokinetic vs. ppc1
glm1 <- glm(class ~ PC1, data= pca.limo.points.class, family = "binomial")
summary(glm1)
phyloglm1 <- phyloglm(class ~ PC1, data= pca.limo.points.class, phy= tree_classifier_26tips, method='logistic_IG10', btol = 11)
summary(phyloglm1)
#phylogenetic logistic regression limokinetic vs. ppc2
glm2 <- glm(class ~ PC2, data= pca.limo.points.class, family = "binomial")
summary(glm2)
phyloglm2 <- phyloglm(class ~ PC2, data= pca.limo.points.class, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm2)
#phylogenetic logistic regression limokinetic vs. ppc3
glm3 <- glm(class ~ PC3, data= pca.limo.points.class, family = "binomial")
summary(glm3)
phyloglm3 <- phyloglm(class ~ PC3, data= pca.limo.points.class, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm3)
#phylogenetic logistic regression limokinetic vs. ppc4
glm4 <- glm(class ~ PC4, data= pca.limo.points.class, family = "binomial")
summary(glm4)
phyloglm4 <- phyloglm(class ~ PC4, data= pca.limo.points.class, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm4)
#phylogenetic logistic regression limokinetic vs. ppc5
glm5 <- glm(class ~ PC5, data= pca.limo.points.class, family = "binomial")
summary(glm5)
phyloglm5 <- phyloglm(class ~ PC5, data= pca.limo.points.class, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm5)
#phylogenetic logistic regression limokinetic vs. ppc6
glm6 <- glm(class ~ PC6, data= pca.limo.points.class, family = "binomial")
summary(glm6)
phyloglm6 <- phyloglm(class ~ PC6, data= pca.limo.points.class, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm6)