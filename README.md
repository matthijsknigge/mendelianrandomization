# Mendelian Randomization

Mendelian Randomization (MR) refers to the random segregation and assortment of genes from parents to offspring that occur during gamete formation and provides a method of using genetic variants in observational settings to make casual inferences regarding the relationship between exposure and outcomes. The basic principle utilized in the MR framework is that if genetic variants either alter the level of or mirror the biological effects of a modifiable exposure that itself alters disease risk, then these genetic variants should be related to disease risk. MR studies aim to provide evidence for or against a causal relationship between a modifiable exposure variable and a disease. Genetic variants are used because these are less susceptible for confounding because of it is subjected to Mendel’s first law, the law of of segregation. These genetic variants segregate independently and randomly from environmental factors, and it can be assumed that genetic variants segregate independently from other traits.

This package provides functionality for the following operations:
  * Calculation the standard deviation from the effectsize or log odd score when it is not present.
  
  * Clumping, for pruning SNPs that are in linkage disequilibrium (LD). ALso is provided a method for finding proxy SNPs for replacing SNPs that are in LD.
  
  * Cochran's Q test, for determining SNPs that overfit the model, or the ones that introduce pleiotropic effects.
  
  * Mendelian Randomization Egger method (MR-egger method) for estimating causality, testing causality, and for testing the overall pleiotropy within the data set.
  
  * When your data set misses allelic information, this can be queried by using a reference file.
  
  * Forest plotting the MR analysis, for seeing the overall weight a SNP brings into the study.
  
  * Funnel plotting the methods used with the MR analysis to detect study bias.
  
  * Get the chromosome number and position of SNPs.
  
  * Harmonization of the data set. Align SNPs, remove problematic SNPs, for example palindromic SNPs, mismatch SNPs, and SNPs that have a wrong reference.
  
  * Inverse-Variance Weighted method for averaging the estimate ratios.
  
  * A highly fasionable way of plotting MR results.
  
  * The functionality to pre-process data, which test the data set on missing alleles, missing beta's, selects for genome-wide significance, removes duplicates,  and removes alleles from which it is not possible to measure direction.
  
  * qq-plot for the p-value distribution of a chosen method.
  
  * A normal qq-plot for plotting the theoretical quantiles against the normal quantiles.
  
  * Remove a certain region within a chromosome.
  
  * Perform wald-ratio for obtaining a causal estimate based on the exposure regression on genotype and the outcome regression of genotype.
  
  * Test data for trying the package.

# Installing Mendelian Randomization
The package is hosted on bitbucket, and this allows for a smooth installation, and updates are easy to install. Before installen Mendelian Randomization, make sure you have installed <code>devtools<code>:











