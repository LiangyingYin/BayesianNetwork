# BN-GWAS
BN-GWAS is a framework to estimate the gene-phenotype network and to quantify the direct and indirect causal effects of genes on the studied disease/trait. In addition, joint causal effects (the causal effects of each gene when intervention is performed on a group of genes) are computed under the same framework, leveraging GWAS-predicted expression profiles and raw expression data from a reference sample (e.g., GTEx).The framework is based on a Bayesian network (BN) approach, which produces a directed graph showing the causal relationships among genes and between genes and the phenotype. It aims to identify the overall causal network and quantify the effects of genes on the studied trait, including both direct and indirect effects. The inferred causal network can be used to predict the consequences of external interventions on various genes, in isolation or combination, on target traits. 

# Source code
# System Requirement
The R source codes are expected to work under Mac, Linux and Windows operation systems.The codes have been tested on R version 3.4.4 and 3.6.1 and are expected to work on newer R versions as well. 
# Installation guide
To successfully run the codes, you need to install all dependency R packages listed in the Rcode files, e.g., pcalg, glassoFast, graph,coop, ect. Typical install time depends on how many packages you need to install in R. The installation time for all required R pacakges should be less than 1 hour on a desktop computer. To speed up the installation process, you are recommended to install miniconda(https://docs.conda.io/en/latest/miniconda.html) and mamba on your desktop(https://github.com/mamba-org/mamba). Please be kindly noted that traditional way of package installation also works.

# Getting started 
# Simulation
The code is available under the Simulation folder. 

# Applications
The code is available under the Application folder. Due to privacy issues, the original data from UK Biobank is not provided here.

In order to use BN-GWAS to identify gene-phenotype network for trait of interest, you need to separately infer the gene-gene, gene-outcome network first, then merge them into a complete gene-phenotype network.

# Causalgraphs
This folder contains the inferred causal graphs for some studied traits.

# License
This project is licensed under GNU GPL v3.

# Authors
Liangying Yin(The Chinese University of Hong Kong)

Yaning Feng(The Chinese University of Hong Kong)

Hon-Cheong So(The Chinese University of Hong Kong)
