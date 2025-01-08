#logistic regression for orthologs in limostatic motility strategy classifier
#
#load packages
require(tidyverse)
require(broom)
#
#
#load limostatic data
traitdata_static <- read.csv('dataset_limostatic_KO.csv', header = T)
#convert traitdata_static kegg data to binary variables
dim(traitdata_static)
traitdata_static[15:20,1:10]
traitdata_static[,c(3:125)] <- ifelse(traitdata_static[,c(3:125)] > 0,1,0)
traitdata_static[15:20,1:10]
#
#
#data prep for regression, convert variable Tipname to rownames
row.names(traitdata_static) <- traitdata_static$Tipname
#remove tipname as a variable now that it is a row name
traitdata_static <- traitdata_static[,-c(1)]
#
#visually inspect data
traitdata_static[1:5,1:10]
#
#log reg
output <- traitdata_static %>%
	names() %>%
	paste("class ~", .) %>%
	map_df(~glm(as.formula(.x), data= traitdata_static, family = "binomial") %>% tidy())
#examine output
output
#change colnames in output
colnames(output) <- c("term", "Estimate", "StdErr", "z.value", "p.value")
#remove the first row from output, which represents a model only considering the class variable (no slope fitted) so only an intercept and makes the paired data following it confusing to interpret
output_trim <- output[-1,]
output_trim
#	
write.csv(output_trim, row.names=F, file = "output_logreg_limostatic_KO.csv")
