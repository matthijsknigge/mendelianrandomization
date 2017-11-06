# Mendelian Randomization

---

Mendelian Randomization (MR) refers to the random segregation and assortment of genes from parents to offspring that occur during gamete formation and provides a method of using genetic variants in observational settings to make casual inferences regarding the relationship between exposure and outcomes. The basic principle utilized in the MR framework is that if genetic variants either alter the level of or mirror the biological effects of a modifiable exposure that itself alters disease risk, then these genetic variants should be related to disease risk. MR studies aim to provide evidence for or against a causal relationship between a modifiable exposure variable and a disease. Genetic variants are used because these are less susceptible for confounding because of it is subjected to Mendel’s first law, the law of of segregation. These genetic variants segregate independently and randomly from environmental factors, and it can be assumed that genetic variants segregate independently from other traits.

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

# Tutorial

---

The package also contains test data for doing a basic Mendelian Randomization analysis. The first step is to read the data. For this analysis we want to infer causality between an exposure and outcome. In this setup the exposure is Celiac Disease, and the outcome is High-density-lipoproteïne (HDL).

```
# the exposure
exposure <- data("celiac")

# the outcome
outcome  <- data("hdl")
```

Lets check out the data.
```
head(exposure)
```
|SNPs       | effect_allele | Z_OR       | P       | se          |
|-----------|---------------|------------|---------|-------------|
|rs61733845 | T             | 0.03536714 |4.249e-01| 0.04432255  |
|rs1320571  | A             | 0.01882175 |6.590e-01| 0.04265126  |
|rs9729550  | A             | 0.10048354 |2.465e-06| 0.02132954  |
|rs1815606  | G             | 0.06774365 |7.151e-04| 0.02002045  |
|rs7515488  | T             | -0.10280822|1.195e-04| 0.02672322  |
|rs11260562 | A             | -0.03936472|3.344e-01| 0.04078024  |


```
head(outcome)
```

| SNP_hg19  |  rsid    |A1 |A2 | Beta |    SE  | N     |P.value   | Freq.A1.ESP.EUR |
|-----------|----------|---|---|------|--------|-------|----------|-----------------|
|15:58723675 |rs1800588|  T|  C| 0.1180589| 0.002991| 316391|       0  |       0.24468|
|16:56985139 |rs9989419|  G | A| 0.1306400| 0.002602| 316391|       0 |        0.60004|
|16:56988044|  rs173539|  T | C| 0.2301930| 0.002842| 292803|       0    |     0.32579|
|16:56989590 | rs247616|  T | C| 0.2424130| 0.002758| 316391|       0   |      0.30817|
|16:56993324 |rs3764261|  A | C| 0.2387770| 0.002775| 310132|       0  |       0.31250|
|16:56995236 |rs1800775|  A | C| 0.1901360| 0.002690| 285413|       0 |        0.50581|




