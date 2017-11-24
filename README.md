# Mendelian Randomization

---

Mendelian Randomization (MR) is the process that refers to the random segregation and assortment of genes from ancestors to offspring that takes place during gamete formation and gives a method of using genetic variants to make casual inferences regarding the relationship between exposure and outcomes. The basic principle utilized in the MR pipeline, is that if a genetic variant either alters the level of or mimics the biological effects of a exposure that itself alters disease risk, then these genetic variants should be related to disease risk.The goal of MR studies is to provide evidence for or against a causal relationship between a exposure and a disease. Genetic variants are used because these are less susceptible for confounding because of it is subjected to Mendel’s first law, the law of of segregation. These genetic variants segregate independently and randomly from environmental factors, and it can be assumed that genetic variants segregate independently from other traits.

This package provides functionality for the following operations:

  * Calculate the standard deviation from the effectsize or log odd score when it is not present. `mr.calculate.se()`
  
  * Clumping, for pruning SNPs that are in linkage disequilibrium (LD). ALso is provided a method for finding proxy SNPs for replacing SNPs that are in LD. `mr.clump()`
  
  * Cochran's Q test, for determining SNPs that overfit the model, or the ones that introduce pleiotropic effects. `mr.cochran.Q.test()`
  
  * Mendelian Randomization Egger method (MR-egger method) for estimating causality, testing causality, and for testing the overall pleiotropy within the data set. `mr.egger.method()`
  
  * When your data set misses allelic information, this can be queried by using a reference file. `mr.find.missing.allelic.information()`
  
  * Forest plotting the MR analysis, for seeing the overall weight a SNP brings into the study. `mr.forest.plot()`
  
  * Funnel plotting the methods used with the MR analysis to detect study bias. `mr.funnel.plot()`
  
  * Get the chromosome number and position of SNPs. `mr.get.chr.pos()`
  
  * Harmonization of the data set. Align SNPs, remove problematic SNPs, for example palindromic SNPs, mismatch SNPs, and SNPs that have a wrong reference. `mr.harmonize()`
  
  * Inverse-Variance Weighted method for averaging the estimate ratios. `mr.inverse.variance.weighted.method()`
  
  * A highly fasionable way of plotting MR results. `mr.plot()`
  
  * The functionality to pre-process data, which test the data set on missing alleles, missing beta's, selects for genome-wide significance, removes duplicates,  and removes alleles from which it is not possible to measure direction. `mr.pre.process()`
  
  * qq-plot for the p-value distribution of a chosen method. `mr.qq.p.distribution()`
  
  * A normal qq-plot for plotting the theoretical quantiles against the normal quantiles. `mr.qq.plot()`
  
  * Remove a certain region within a chromosome. `mr.remove.region()`
  
  * Perform wald-ratio for obtaining a causal estimate based on the exposure regression on genotype and the outcome regression of genotype. `mr.wald.ratio()`
  
  * Test data for trying the package.

# Installing Mendelian Randomization

---

The package is hosted on bitbucket, and this allows for a smooth installation, and updates are easy to install. Before installing Mendelian Randomization, make sure you have installed `devtools`:

```
install.packages("devtools")
```

And then you are ready to install the `mendelianRandomization` package:

```
devtools::install_bitbucket("matthijsknigge/mendelianRandomization")
```

Other libraries that are needed in this package:

```
install.packages("stringr")
install.packages("readr")
install.packages("ggplot2")
install.packages("ggExtra")
install.packages("gridExtra")
install.packages("latex2exp")
```
This package needs R version 3.2.0 or greater.

# Tutorial

---

The package also contains test data for doing a basic Mendelian Randomization analysis. The first step is to read the data. For this analysis we want to infer causality between an exposure and outcome. In this setup the exposure is Celiac Disease, and the outcome is High-density-lipoproteïne (HDL).

```
# the exposure
exposure <- data("celiac")

# the outcome
outcome  <- data("hdl")
```

# Tutorial

---

The package also contains test data for doing a basic Mendelian Randomization analysis. The first step is to read the data. For this analysis we want to infer causality between an exposure and outcome. In this setup the exposure is Inflammatory bowel disease, and the outcome is Celiac Disease.

```{r}
# the exposure
data("celiac")
outcome <- celiac

# the outcome
data("Inflammatory.bowel.disease")
exposure  <- Inflammatory.bowel.disease
```


Lets check out the data.
```
head(outcome)
```
Here we se a column with SNP identifiers, the effect allele, the effectsize, the pvalue, and the standard deviation.
|SNP        |effect_allele |       Z_OR|         P|        se|
|:----------|:-------------|----------:|---------:|---------:|
|rs61733845 |T             |  0.0353671| 0.4249000| 0.0443226|
|rs1320571  |A             |  0.0188218| 0.6590000| 0.0426513|
|rs9729550  |A             |  0.1004835| 0.0000025| 0.0213295|
|rs1815606  |G             |  0.0677437| 0.0007151| 0.0200204|
|rs7515488  |T             | -0.1028082| 0.0001195| 0.0267232|
|rs11260562 |A             | -0.0393647| 0.3344000| 0.0407802|





