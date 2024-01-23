# SLE-Gene_Expression-Analysis
Publicly available expression profiling by array data was used to examine differential gene expression in samples with systemic lupus erythematosus and controls.

Data was derived from GSE17755 and methods were modified from https://doi.org/10.1186/s12967-023-03943-9 . 

This script uses a python implementation of limma with a Bayesian adjustment and also performs principle component analysis. A significance threshold is defined as log2(FC)|> 1.5 and p-value < 0.05. Results can be viewed in a volcano plot and there is a counting function to determine the most common GO terms associated with significant genes. 
