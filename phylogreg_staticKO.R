#phylogenetic regression of limostatic motility strategy KEGG orthologs
#
#load packages
require(phylolm)
require(phytools)
#
#
#load phylogenetic tree
tree_full <- read.tree("dataset_fulltree.treefile",tree.names=T)
tree_full
#
#
#load limokinetic data
traitdata_static <- read.csv('dataset_limostatic_KO.csv', header = T)
#
#convert traitdata_static kegg data to binary variables
dim(traitdata_static)
traitdata_static[15:20,1:10]
traitdata_static[,c(3:125)] <- ifelse(traitdata_static[,c(3:125)] > 0,1,0)
traitdata_static[15:20,1:10]
#
#
#data prep for phylogenetic regression
#keep only tips with trait data
tree_classifier_26tips <- keep.tip(tree_full, tip = c(traitdata_static$Tipname))
tree_classifier_26tips
#
#plot trimmed tree
plot(tree_classifier_26tips) 
#
#tree tip labels as rownames
row.names(traitdata_static) <- traitdata_static$Tipname
#
#visually inspect data
head(traitdata_static)
#
#
#phylogenetic logistic regression limostatic classifier outcome vs. all features
#
phyloglm1 <- phyloglm(class ~ traitdata_static[,3], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm1)
#
phyloglm2 <- phyloglm(class ~ traitdata_static[,4], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm2)
#
phyloglm3 <- phyloglm(class ~ traitdata_static[,5], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm3)
#
phyloglm4 <- phyloglm(class ~ traitdata_static[,6], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm4)
#
phyloglm5 <- phyloglm(class ~ traitdata_static[,7], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm5)
#
phyloglm6 <- phyloglm(class ~ traitdata_static[,8], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm6)
#
phyloglm7 <- phyloglm(class ~ traitdata_static[,9], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm7)
#
phyloglm8 <- phyloglm(class ~ traitdata_static[,10], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm8)
#
phyloglm9 <- phyloglm(class ~ traitdata_static[,11], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm9)
#
phyloglm10 <- phyloglm(class ~ traitdata_static[,12], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm10)
#
phyloglm11 <- phyloglm(class ~ traitdata_static[,13], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm11)
#
phyloglm12 <- phyloglm(class ~ traitdata_static[,14], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm12)
#
phyloglm13 <- phyloglm(class ~ traitdata_static[,15], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm13)
#
phyloglm14 <- phyloglm(class ~ traitdata_static[,16], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm14)
#
phyloglm15 <- phyloglm(class ~ traitdata_static[,17], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm15)
#
phyloglm16 <- phyloglm(class ~ traitdata_static[,18], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm16)
#
phyloglm17 <- phyloglm(class ~ traitdata_static[,19], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm17)
#
phyloglm18 <- phyloglm(class ~ traitdata_static[,20], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm18)
#
phyloglm19 <- phyloglm(class ~ traitdata_static[,21], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm19)
#
phyloglm20 <- phyloglm(class ~ traitdata_static[,22], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm20)
#
phyloglm21 <- phyloglm(class ~ traitdata_static[,23], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm21)
#
phyloglm22 <- phyloglm(class ~ traitdata_static[,24], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm22)
#
phyloglm23 <- phyloglm(class ~ traitdata_static[,25], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm23)
#
phyloglm24 <- phyloglm(class ~ traitdata_static[,26], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm24)
#
phyloglm25 <- phyloglm(class ~ traitdata_static[,27], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm25)
#
phyloglm26 <- phyloglm(class ~ traitdata_static[,28], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm26)
#
phyloglm27 <- phyloglm(class ~ traitdata_static[,29], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm27)
#
phyloglm28 <- phyloglm(class ~ traitdata_static[,30], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm28)
#
phyloglm29 <- phyloglm(class ~ traitdata_static[,31], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm29)
#
phyloglm30 <- phyloglm(class ~ traitdata_static[,32], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm30)
#
phyloglm31 <- phyloglm(class ~ traitdata_static[,33], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm31)
#
phyloglm32 <- phyloglm(class ~ traitdata_static[,34], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm32)
#
phyloglm33 <- phyloglm(class ~ traitdata_static[,35], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm33)
#
phyloglm34 <- phyloglm(class ~ traitdata_static[,36], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm34)
#
phyloglm35 <- phyloglm(class ~ traitdata_static[,37], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm35)
#
phyloglm36 <- phyloglm(class ~ traitdata_static[,38], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm36)
#
phyloglm37 <- phyloglm(class ~ traitdata_static[,39], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm37)
#
phyloglm38 <- phyloglm(class ~ traitdata_static[,40], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm38)
#
phyloglm39 <- phyloglm(class ~ traitdata_static[,41], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm39)
#
phyloglm40 <- phyloglm(class ~ traitdata_static[,42], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm40)
#
phyloglm41 <- phyloglm(class ~ traitdata_static[,43], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm41)
#
phyloglm42 <- phyloglm(class ~ traitdata_static[,44], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm42)
#
phyloglm43 <- phyloglm(class ~ traitdata_static[,45], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm43)
#
phyloglm44 <- phyloglm(class ~ traitdata_static[,46], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm44)
#
phyloglm45 <- phyloglm(class ~ traitdata_static[,47], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm45)
#
phyloglm46 <- phyloglm(class ~ traitdata_static[,48], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm46)
#
phyloglm47 <- phyloglm(class ~ traitdata_static[,49], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm47)
#
phyloglm48 <- phyloglm(class ~ traitdata_static[,50], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm48)
#
phyloglm49 <- phyloglm(class ~ traitdata_static[,51], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm49)
#
phyloglm50 <- phyloglm(class ~ traitdata_static[,52], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm50)
#
phyloglm51 <- phyloglm(class ~ traitdata_static[,53], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm51)
#
phyloglm52 <- phyloglm(class ~ traitdata_static[,54], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm52)
#
phyloglm53 <- phyloglm(class ~ traitdata_static[,55], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm53)
#
phyloglm54 <- phyloglm(class ~ traitdata_static[,56], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm54)
#
phyloglm55 <- phyloglm(class ~ traitdata_static[,57], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm55)
#
phyloglm56 <- phyloglm(class ~ traitdata_static[,58], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm56)
#
phyloglm57 <- phyloglm(class ~ traitdata_static[,59], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm57)
#
phyloglm58 <- phyloglm(class ~ traitdata_static[,60], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm58)
#
phyloglm59 <- phyloglm(class ~ traitdata_static[,61], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm59)
#
phyloglm60 <- phyloglm(class ~ traitdata_static[,62], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm60)
#
phyloglm61 <- phyloglm(class ~ traitdata_static[,63], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm61)
#
phyloglm62 <- phyloglm(class ~ traitdata_static[,64], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm62)
#
phyloglm63 <- phyloglm(class ~ traitdata_static[,65], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm63)
#
phyloglm64 <- phyloglm(class ~ traitdata_static[,66], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm64)
#
phyloglm65 <- phyloglm(class ~ traitdata_static[,67], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm65)
#
phyloglm66 <- phyloglm(class ~ traitdata_static[,68], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm66)
#
phyloglm67 <- phyloglm(class ~ traitdata_static[,69], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm67)
#
phyloglm68 <- phyloglm(class ~ traitdata_static[,70], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm68)
#
phyloglm69 <- phyloglm(class ~ traitdata_static[,71], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm69)
#
phyloglm70 <- phyloglm(class ~ traitdata_static[,72], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm70)
#
phyloglm71 <- phyloglm(class ~ traitdata_static[,73], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm71)
#
phyloglm72 <- phyloglm(class ~ traitdata_static[,74], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm72)
#
phyloglm73 <- phyloglm(class ~ traitdata_static[,75], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm73)
#
phyloglm74 <- phyloglm(class ~ traitdata_static[,76], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm74)
#
phyloglm75 <- phyloglm(class ~ traitdata_static[,77], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm75)
#
phyloglm76 <- phyloglm(class ~ traitdata_static[,78], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm76)
#
phyloglm77 <- phyloglm(class ~ traitdata_static[,79], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm77)
#
phyloglm78 <- phyloglm(class ~ traitdata_static[,80], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm78)
#
phyloglm79 <- phyloglm(class ~ traitdata_static[,81], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm79)
#
phyloglm80 <- phyloglm(class ~ traitdata_static[,82], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm80)
#
phyloglm81 <- phyloglm(class ~ traitdata_static[,83], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm81)
#
phyloglm82 <- phyloglm(class ~ traitdata_static[,84], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm82)
#
phyloglm83 <- phyloglm(class ~ traitdata_static[,85], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm83)
#
phyloglm84 <- phyloglm(class ~ traitdata_static[,86], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm84)
#
phyloglm85 <- phyloglm(class ~ traitdata_static[,87], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm85)
#
phyloglm86 <- phyloglm(class ~ traitdata_static[,88], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm86)
#
phyloglm87 <- phyloglm(class ~ traitdata_static[,89], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm87)
#
phyloglm88 <- phyloglm(class ~ traitdata_static[,90], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm88)
#
phyloglm89 <- phyloglm(class ~ traitdata_static[,91], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm89)
#
phyloglm90 <- phyloglm(class ~ traitdata_static[,92], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm90)
#
phyloglm91 <- phyloglm(class ~ traitdata_static[,93], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm91)
#
phyloglm92 <- phyloglm(class ~ traitdata_static[,94], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm92)
#
phyloglm93 <- phyloglm(class ~ traitdata_static[,95], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm93)
#
phyloglm94 <- phyloglm(class ~ traitdata_static[,96], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm94)
#
phyloglm95 <- phyloglm(class ~ traitdata_static[,97], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm95)
#
phyloglm96 <- phyloglm(class ~ traitdata_static[,98], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm96)
#
phyloglm97 <- phyloglm(class ~ traitdata_static[,99], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm97)
#
phyloglm98 <- phyloglm(class ~ traitdata_static[,100], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm98)
#
phyloglm99 <- phyloglm(class ~ traitdata_static[,101], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm99)
#
phyloglm100 <- phyloglm(class ~ traitdata_static[,102], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm100)
#
phyloglm101 <- phyloglm(class ~ traitdata_static[,103], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm101)
#
phyloglm102 <- phyloglm(class ~ traitdata_static[,104], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm102)
#
phyloglm103 <- phyloglm(class ~ traitdata_static[,105], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm103)
#
phyloglm104 <- phyloglm(class ~ traitdata_static[,106], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm104)
#
phyloglm105 <- phyloglm(class ~ traitdata_static[,107], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm105)
#
phyloglm106 <- phyloglm(class ~ traitdata_static[,108], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm106)
#
phyloglm107 <- phyloglm(class ~ traitdata_static[,109], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm107)
#
phyloglm108 <- phyloglm(class ~ traitdata_static[,110], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm108)
#
phyloglm109 <- phyloglm(class ~ traitdata_static[,111], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm109)
#
phyloglm110 <- phyloglm(class ~ traitdata_static[,112], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm110)
#
phyloglm111 <- phyloglm(class ~ traitdata_static[,113], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm111)
#
phyloglm112 <- phyloglm(class ~ traitdata_static[,114], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE', btol = 21)
summary(phyloglm112)
#
phyloglm113 <- phyloglm(class ~ traitdata_static[,115], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm113)
#
phyloglm114 <- phyloglm(class ~ traitdata_static[,116], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm114)
#
phyloglm115 <- phyloglm(class ~ traitdata_static[,117], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm115)
#
phyloglm116 <- phyloglm(class ~ traitdata_static[,118], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm116)
#
phyloglm117 <- phyloglm(class ~ traitdata_static[,119], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm117)
#
phyloglm118 <- phyloglm(class ~ traitdata_static[,120], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm118)
#
phyloglm119 <- phyloglm(class ~ traitdata_static[,121], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm119)
#
phyloglm120 <- phyloglm(class ~ traitdata_static[,122], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm120)
#
phyloglm121 <- phyloglm(class ~ traitdata_static[,123], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm121)
#
phyloglm122 <- phyloglm(class ~ traitdata_static[,124], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm122)
#
phyloglm123 <- phyloglm(class ~ traitdata_static[,125], data= traitdata_static, phy= tree_classifier_26tips, method='logistic_MPLE')
summary(phyloglm123)
#
#
#output all phyloglm coefficients as a table
plogreg.all <- rbind(summary(phyloglm1)$coefficients, summary(phyloglm2)$coefficients, summary(phyloglm3)$coefficients, summary(phyloglm4)$coefficients, summary(phyloglm5)$coefficients, summary(phyloglm6)$coefficients, summary(phyloglm7)$coefficients, summary(phyloglm8)$coefficients, summary(phyloglm9)$coefficients, summary(phyloglm10)$coefficients, summary(phyloglm11)$coefficients, summary(phyloglm12)$coefficients, summary(phyloglm13)$coefficients, summary(phyloglm14)$coefficients, summary(phyloglm15)$coefficients, summary(phyloglm16)$coefficients, summary(phyloglm17)$coefficients, summary(phyloglm18)$coefficients, summary(phyloglm19)$coefficients, summary(phyloglm20)$coefficients, summary(phyloglm21)$coefficients, summary(phyloglm22)$coefficients, summary(phyloglm23)$coefficients, summary(phyloglm24)$coefficients, summary(phyloglm25)$coefficients, summary(phyloglm26)$coefficients, summary(phyloglm27)$coefficients, summary(phyloglm28)$coefficients, summary(phyloglm29)$coefficients, summary(phyloglm30)$coefficients, summary(phyloglm31)$coefficients, summary(phyloglm32)$coefficients, summary(phyloglm33)$coefficients, summary(phyloglm34)$coefficients, summary(phyloglm35)$coefficients, summary(phyloglm36)$coefficients, summary(phyloglm37)$coefficients, summary(phyloglm38)$coefficients, summary(phyloglm39)$coefficients, summary(phyloglm40)$coefficients, summary(phyloglm41)$coefficients, summary(phyloglm42)$coefficients, summary(phyloglm43)$coefficients, summary(phyloglm44)$coefficients, summary(phyloglm45)$coefficients, summary(phyloglm46)$coefficients, summary(phyloglm47)$coefficients, summary(phyloglm48)$coefficients, summary(phyloglm49)$coefficients, summary(phyloglm50)$coefficients, summary(phyloglm51)$coefficients, summary(phyloglm52)$coefficients, summary(phyloglm53)$coefficients, summary(phyloglm54)$coefficients, summary(phyloglm55)$coefficients, summary(phyloglm56)$coefficients, summary(phyloglm57)$coefficients, summary(phyloglm58)$coefficients, summary(phyloglm59)$coefficients, summary(phyloglm60)$coefficients, summary(phyloglm61)$coefficients, summary(phyloglm62)$coefficients, summary(phyloglm63)$coefficients, summary(phyloglm64)$coefficients, summary(phyloglm65)$coefficients, summary(phyloglm66)$coefficients, summary(phyloglm67)$coefficients, summary(phyloglm68)$coefficients, summary(phyloglm69)$coefficients, summary(phyloglm70)$coefficients, summary(phyloglm71)$coefficients, summary(phyloglm72)$coefficients, summary(phyloglm73)$coefficients, summary(phyloglm74)$coefficients, summary(phyloglm75)$coefficients, summary(phyloglm76)$coefficients, summary(phyloglm77)$coefficients, summary(phyloglm78)$coefficients, summary(phyloglm79)$coefficients, summary(phyloglm80)$coefficients, summary(phyloglm81)$coefficients, summary(phyloglm82)$coefficients, summary(phyloglm83)$coefficients, summary(phyloglm84)$coefficients, summary(phyloglm85)$coefficients, summary(phyloglm86)$coefficients, summary(phyloglm87)$coefficients, summary(phyloglm88)$coefficients, summary(phyloglm89)$coefficients, summary(phyloglm90)$coefficients, summary(phyloglm91)$coefficients, summary(phyloglm92)$coefficients, summary(phyloglm93)$coefficients, summary(phyloglm94)$coefficients, summary(phyloglm95)$coefficients, summary(phyloglm96)$coefficients, summary(phyloglm97)$coefficients, summary(phyloglm98)$coefficients, summary(phyloglm99)$coefficients, summary(phyloglm100)$coefficients, summary(phyloglm101)$coefficients, summary(phyloglm102)$coefficients, summary(phyloglm103)$coefficients, summary(phyloglm104)$coefficients, summary(phyloglm105)$coefficients, summary(phyloglm106)$coefficients, summary(phyloglm107)$coefficients, summary(phyloglm108)$coefficients, summary(phyloglm109)$coefficients, summary(phyloglm110)$coefficients, summary(phyloglm111)$coefficients, summary(phyloglm112)$coefficients, summary(phyloglm113)$coefficients, summary(phyloglm114)$coefficients, summary(phyloglm115)$coefficients, summary(phyloglm116)$coefficients, summary(phyloglm117)$coefficients, summary(phyloglm118)$coefficients, summary(phyloglm119)$coefficients, summary(phyloglm120)$coefficients, summary(phyloglm121)$coefficients, summary(phyloglm122)$coefficients, summary(phyloglm123)$coefficients)
#check data frame
head(plogreg.all)
#make new column to signify the term (slope or intercept) and the KO (always paired sets of the KO indicated by the slope term)
term <- rownames(plogreg.all)
#convert output to a data frame
plogreg.all.df <- data.frame(plogreg.all)
#add the new column term to the output dataframe
plogreg.all.df$term <- term
#change order of columns in the dataset
plogreg.all.df <- plogreg.all.df[,c(5,1:4)]
#inspect dataset
plogreg.all.df
#replace the nondescriptive generic name of each slope term (the bracketed position name in the original dataset) with the KO name from the original dataset
plogreg.all.df$term[seq(2,246, by=2)] <-c(colnames(traitdata_static[3:125])) 
head(plogreg.all.df)
write.csv(plogreg.all.df, row.names=F, "output_PhyLogReg_limostatic_KO.csv")