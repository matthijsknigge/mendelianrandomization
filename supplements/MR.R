# libraries
library(mendelianRandomization)
library(data.table)

# outcome
# o <- fread("~/Desktop/eo_baso_sum_build37_171771_20161212.txt", header = T)
o <- fread("~/Bioinformatics/Src/R/cluster_results/GLGC_exome/data/GLGC_exome_TC_release.txt", header = T)

# exposure
e <- fread("~/Bioinformatics/Src/R/cluster_results/outcome_without_HLL_region/outcome/celiac_removed/celiac_2011.txt", header = T, sep = " ")

# calc se
e$se <- mr.calculate.se(B = e$Z_OR, pval = e$P)
o$se <- mr.calculate.se(B = o$Z_OR, pval = o$P)

# pre process
e <- mr.pre.process(B = e$Z_OR, B.se = e$se, pval = e$P, effect_allele = e$effect_allele, other_allele = e$other_allele, SNP = e$SNP)

# harmonize
h <- mr.harmonize(By = o$Beta, Bx = e$beta, By.se = o$SE, Bx.se = e$se, outcome.pval = o$`P-value`, exposure.pval = e$pval, outcome.effect_allele = o$A1, exposure.effect_allele = e$effect_allele, exposure.other_allele = e$other_allele,
                  outcome.SNP = o$rsid, exposure.SNP = e$SNP)

# clump
h <- mr.clump(data = h, refdat = "~/Bioinformatics/1000G/EUR_1000G_phase3_vcf", verbose = T)


# wald
wald <- mr.wald.ratio(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)

# ivw
ivw <- mr.inverse.variance.weighted.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)



# egger
egger <- mr.egger.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)

# combine
h$iv <- wald$iv; h$iv.se <- wald$iv.se; h$iv.p <- wald$iv.p
#

h$cochran.Q <- mr.cochran.Q.test(data = h, pval = .05)$cochran.Q

# after Q
egger.Q <- mr.egger.method(By = h[which(h$cochran.Q == 1), ]$By, Bx = h[which(h$cochran.Q == 1), ]$Bx, By.se = h[which(h$cochran.Q == 1), ]$By.se, Bx.se = h[which(h$cochran.Q == 1), ]$Bx.se)

# ivw
ivw.Q <- mr.inverse.variance.weighted.method(By = h[which(h$cochran.Q == 1), ]$By, Bx = h[which(h$cochran.Q == 1), ]$Bx, By.se = h[which(h$cochran.Q == 1), ]$By.se, Bx.se = h[which(h$cochran.Q == 1), ]$Bx.se)

# plot

p <- mr.plot(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se, iv = h$iv, iv.se = h$iv.se, ivw = ivw$ivw, egger = egger$egger, egger.i = egger$egger.i, ivw.Q = ivw.Q$ivw,
             egger.Q = egger.Q$egger, egger.i.Q = egger.Q$egger.i, egger.p.fdr = egger$egger.p, ivw.p.fdr = ivw$ivw.p, egger.i.p = egger$egger.i.p, exposure.name =  "Celiac 2011",
             outcome.name = "eo_baso_sum_build37_171771_20161212")
ggdraw(p)

ggsave(filename = "~/Desktop/test.png", plot = ggdraw(p))


ivw




