# BN-GWAS
BN-GWAS is a framework to estimate the gene-phenotype network and to quantify the direct and indirect causal effects of genes on the studied disease/trait. In addition, joint causal effects (the causal effects of each gene when intervention is performed on a group of genes) are computed under the same framework, leveraging GWAS-predicted expression profiles and raw expression data from a reference sample (e.g., GTEx).The framework is based on a Bayesian network (BN) approach, which produces a directed graph showing the causal relationships among genes and between genes and the phenotype. It aims to identify the overall causal network and quantify the effects of genes on the studied trait, including both direct and indirect effects. The inferred causal network can be used to predict the consequences of external interventions on various genes, in isolation or combination, on target traits. 


# Getting started 
# Simulation
The code is available under the Simulation folder

# Applications
The code is available under the Application folder

In order to use BN-GWAS to identify gene-phenotype network for trait of interest, you need to separately infer the gene-gene, gene-outcome network first, then merge them into a complete gene-phenotype network

# Causalgraphs
This folder demonstrates the inferred causal graphs for some studied traits

# License
This project is licensed under GNU GPL v3.

# Authors
Liangying Yin(The Chinese University of Hong Kong)

Yaning Feng(The Chinese University of Hong Kong)

Hon-Cheong So(The Chinese University of Hong Kong)
