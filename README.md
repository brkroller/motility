# Description
Data and scripts for phylogenetic comparative analyses of motility strategies associated with the manuscript: 2024 Keegstra et al. Risk-reward trade-off in motility endurance generates dichotomy in search strategies among copiotrophic marine bacteria.

# System requirements
- R version 4.4.1 (2024-06-14) “Race for Your Life” with packages “tidyverse – 2.0.0”, “broom – 1.0.6”, “phylolm – 2.6.2”, & “phytools – 2.3.0” installed (along with their dependencies).
 on Mac OSX/Windows/Linux operating systems.
- Software tested on R version 4.4.1 installed on macOS Sequoia 15.2
- No required non-standard hardware.

# Installation guide
- Download and install R (https://cran.r-project.org). Open R and install the following packages with the following generic command.
- install.packages(“name of package”, dependencies = T)
- Package names and versions “tidyverse” – v2.0.0, “broom” – v1.0.6”, “phylolm” – v2.6.2, & “phytools” – v2.3.0.
- Download the 9 files in this repository.
- Typical install time on a “normal” desktop computer: minutes.

# Demo (see instructions for use below)
# Instructions for use
- Overview: Each R script corresponds to 1 or more data file(s) in csv format, which together regenerate the analyses in the paper. Scripts and data files are named after their analysis method and dataset. The R scripts reference the data file names as provided in this repository, so if you change the names of the files you will need to change the corresponding file names in the scripts.  Please note that in the R scripts, the load data section will require the user to either: 1) store the R scripts and data files in the same directory, 2) adjust the filename in the script to include the path to the stored data file or 3) reset the working directory in the R console to the path where the data file is stored.

- File names and brief descriptions.
1) dataset_fulltree.treefile = Phylogenetic tree dataset used as input for the scripts. 
2) dataset_limokinetic_KO.csv = KEGG ortholog dataset containing 22 orthologs for 26 genomes which are used by the limokinetic classifier. Each row represents a genome, with the first column representing the name of the genome in the phylogenetic tree. The second column represents the input for the Bayesian classifier (whether the genome is identified as either 1 =the classification indicated in the filename or 0 =the other classification). The third column represents the output from the Bayesian classifier (whether the genome is identified as either 1 =the classification indicated in the filename or 0 =the other classification).  All remaining columns represent the number genes from each orthologous group present in the genome.
3) dataset_limostatic_KO.csv = KEGG ortholog dataset containing 121 orthologs for 26 genomes which are used by the limostatic classifier. Each row represents a genome, with the first column representing the name of the genome in the phylogenetic tree. The second column represents the input for the Bayesian classifier (whether the genome is identified as either 1 =the classification indicated in the filename or 0 =the other classification). The third column represents the output from the Bayesian classifier (whether the genome is identified as either 1 =the classification indicated in the filename or 0 =the other classification).
4) logreg_kinetic_KO.R = R script for performing logistic regression of limokinetic classification outcome against each orthologous group presence/absence without accounting for phylogenetic relationships among genomes.
5) logreg_static_KO.R = R script for performing logistic regression of limostatic classification outcome against each orthologous group presence/absence without accounting for phylogenetic relationships among genomes.
6) phylopca_kineticKO.R = R script for performing phylogenetic PCA of orthologous groups used by the limokinetic classifier and assessing if orthogonal combinations of features were associated with the classifier outcome after accounting for phylogenetic relationships among genomes.
7) phylopca_staticKO.R = R script for performing phylogenetic PCA of orthologous groups used by the limostatic classifier and assessing if orthogonal combinations of features were associated with the classifier outcome after accounting for phylogenetic relationships among genomes.
8) phylogreg_kineticKO.R = R script for performing logistic regression of limokinetic classification outcome against each orthologous group presence/absence while accounting for phylogenetic relationships among genomes.
9) phylogreg_staticKO.R = R script for performing logistic regression of limostatic classification outcome against each orthologous group presence/absence while accounting for phylogenetic relationships among genomes.

- Open the R script file, execute the code to perform the desired analysis.
- Expected output: 6 total comma separated values files containing summary data for logistic regression and phylogenetic logistic regression analysis. Output for phylogenetic PCA is not summarized as 1 output file, but as individual saved objects in the R software. See comments in scripts for descriptions of each saved object.
- Expected run time on a “normal” desktop computer: minutes.
