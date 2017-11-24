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

```
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

|SNP        |effect_allele |Z_OR       |P         |se        |
|-----------|--------------|-----------|----------|----------|
|rs61733845 |T             |  0.0353671| 0.4249000| 0.0443226|
|rs1320571  |A             |  0.0188218| 0.6590000| 0.0426513|
|rs9729550  |A             |  0.1004835| 0.0000025| 0.0213295|
|rs1815606  |G             |  0.0677437| 0.0007151| 0.0200204|
|rs7515488  |T             | -0.1028082| 0.0001195| 0.0267232|
|rs11260562 |A             | -0.0393647| 0.3344000| 0.0407802|


```
head(exposure)
```
Here there is a column with the SNP identifiers, effect_allele, the effectsize, the standard deviation of the genetic effect, and the p value.

|SNP        |effect_allele |      beta|        se|  pval|
|------|-----------|--------------|----------|----------|------|
|13665 |rs1003342  |NA            |        NA|        NA| 0e+00|
|13666 |rs10051722 |A             | 0.0616269| 0.0107204| 0e+00|
|13667 |rs10061469 |A             | 0.0518248| 0.0105946| 1e-06|
|13668 |rs10065637 |G             | 0.0686809| 0.0128937| 1e-07|
|13669 |rs10142466 |NA            |        NA|        NA| 0e+00|
|13671 |rs10486483 |A             | 0.0602257| 0.0102041| 2e-07|

We need both files at least to contain the SNP id, beta, se, pval, effect allele. And the exposure must containt both alleles to infer in what the direction the effect takes place. Since the other allele is missing for the exposure, we have to query it. The alleles are queried with `mr.find.missing.alleles`.

```
exposure <- mr.find.missing.allelic.information(data = exposure, 
                                              thousand.G = "/path/to/reference.bim")
head(exposure)
```
|      |SNP        |effect_allele |      beta|        se|  pval|other_allele |
|------|-----------|--------------|----------|----------|------|-------------|
|13665 |rs1003342  |NA            |        NA|        NA| 0e+00|             |
|13666 |rs10051722 |A             | 0.0616269| 0.0107204| 0e+00|C            |
|13667 |rs10061469 |A             | 0.0518248| 0.0105946| 1e-06|             |
|13668 |rs10065637 |G             | 0.0686809| 0.0128937| 1e-07|             |
|13669 |rs10142466 |NA            |        NA|        NA| 0e+00|             |
|13671 |rs10486483 |A             | 0.0602257| 0.0102041| 2e-07|G            |

Now that we have what we need, we can start pre-processing the data. But first let's check out how many SNPs we have for exposure and outcome.

```
# exposure amount of SNPs
length(exposure$SNP)
> 196

# outcome amount of SNPs
length(outcome$rsid)
> 97434
```
The next step is to pre-process the exposure, and only the exposure because we are interested in the effect of the exposure on the outcome. We select for genome-wide significance, thus the p-value of the exposure must > 5*10^-8. Also the SNPs without effectsize will be removed, the SNPs from which it is impossible to measure the strand will be removed, and the duplicates will be removed. Besides this, also all negative effectsizes are flipped, and their alleles are also flipped. This is done because we want to measure the positive effect of the exposure on the outcome. This is done with `mr.pre.process`.

```
exposure <- mr.pre.process(B = exposure$beta, 
                           B.se = exposure$se, 
                           pval = exposure$pval, 
                           effect_allele = exposure$effect_allele, 
                           other_allele = exposure$other_allele, 
                           SNP = exposure$SNP)

# amount of SNP left
length(exposure$SNP)
> 52

head(exposure)

```
The product of this is a data frame containing only genome-wide significant SNPs, all signs of effect sizes that are negative are flipped because we are interested in the positive effect of the exposure on the outcome, and ambiguous alleles have been removed from the data set.

|SNP        |      beta|        se|effect_allele |other_allele |pval  |
|---|-----------|----------|----------|--------------|-------------|------|
|2  |rs10051722 | 0.0616269| 0.0107204|A             |C            |9e-09 |
|8  |rs10761659 | 0.1538111| 0.0100349|G             |A            |5e-53 |
|26 |rs11641184 | 0.0772895| 0.0102041|A             |C            |1e-14 |
|38 |rs12627970 | 0.1093163| 0.0127551|G             |A            |2e-18 |
|39 |rs12718244 | 0.0761709| 0.0102041|A             |G            |3e-14 |
|40 |rs12720356 | 0.1483323| 0.0204082|C             |A            |4e-16 |

We do not want our data set to suffer from pleiotropic effects which can be reintroduced by linkage disequilibrium. So we will clump our data set, and the SNPs that are in LD with each other will be removed and will be replaced by a suitable proxy. The clump methods has various other parameters that can be used, but here we use the default, this is done with `mr.clump`. Within the table itself nothing changes, it is the same as the previous table shown. This function parameter `verbose` is set to `TRUE` because we want to use suitable proxy SNPs instead of just removing the SNPs. Hence that also an vector of the SNPs of the outcome is given as parameter, this is because a proxy SNP must be present in both files.

```
exposure <- mr.clump(data = exposure, 
                     refdat = "/path/to/reference", 
                     verbose = T, 
                     SNPs.of.opposite.file = outcome$SNP)
length(exposure$SNP)
> 44
```

The next step is to combine both data sets, with the `mr.harmonize` function. This functions aligns the outcome on the exposure, and takes the intercept between the exposure and outcome.

```
h <- mr.harmonize(By = outcome$Z_OR, Bx = exposure$beta, 
                  By.se = outcome$se, Bx.se = exposure$se, 
                  outcome.pval = outcome$P, exposure.pval = exposure$pval, 
                  outcome.effect_allele = outcome$effect_allele, 
                  exposure.effect_allele = exposure$effect_allele, 
                  exposure.other_allele = exposure$other_allele, 
                  outcome.SNP = outcome$SNP, 
                  exposure.SNP = exposure$SNP)
```

Now it is time to estimate the causalty bewteen regression of the outcome on the genotype and the regression from the exposure on the genotype. Here we use the `mr.wald.method` to calculate the causal estimate between genetic variants. This method return the causal estimate, the standard error of the causal estimate, and the p-value.

```
```
