# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'



data("celiac")
data("hdl")


hdl <- mr.pre.process(B = hdl$Beta, B.se = hdl$SE, pval = hdl$P.value, effect_allele = hdl$A1, other_allele = hdl$A2, SNP = hdl$rsid)

h <- mr.harmonize(By = celiac$Z_OR, Bx = hdl$beta, By.se = celiac$se, Bx.se = hdl$se, outcome.pval = celiac$P, exposure.pval = hdl$pval, outcome.effect_allele = celiac$GiantMinorAllele,
                  exposure.effect_allele = hdl$effect_allele, exposure.other_allele = hdl$other_allele, outcome.SNP = celiac$rs, exposure.SNP = hdl$SNP)

