#phylogenetic regression motility strategies
#
#load packages
require(tidyverse)
require(broom)
#
#
#load limokinetic data
traitdata_kinetic <- read.csv('dataset_limokinetic_KO.csv', header = T)
#convert traitdata_kinetic kegg data to binary variables
dim(traitdata_kinetic)
traitdata_kinetic[,c(3:24)] <- ifelse(traitdata_kinetic[,c(3:24)] > 0,1,0)
head(traitdata_kinetic)
#
#
#data prep for regression, convert variable Tipname to rownames
row.names(traitdata_kinetic) <- traitdata_kinetic$Tipname
#remove Tipname as a variable 
traitdata_kinetic <- traitdata_kinetic[,-c(1)]
#
#visually inspect data
head(traitdata_kinetic)
#
#log reg
output <- traitdata_kinetic %>%
	names() %>%
	paste("class ~", .) %>%
	map_df(~glm(as.formula(.x), data= traitdata_kinetic, family = "binomial") %>% tidy())
#examine output
output
#change colnames in output
colnames(output) <- c("term", "Estimate", "StdErr", "z.value", "p.value")
#remove the first row from output, which represents a model only considering the class variable (no slope fitted) so only an intercept and makes the paired data following it confusing to interpret
output_trim <- output[-1,]
output_trim
write.csv(output_trim, row.names=F, file = "output_logreg_limokinetic_KO.csv")