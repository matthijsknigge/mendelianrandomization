# libraries
library(mendelianRandomization)
library(data.table)

# outcome
# o <- fread("~/Desktop/eo_baso_sum_build37_171771_20161212.txt", header = T)
e <- fread("~/Bioinformatics/Src/R/cluster_results/GLGC_exome/data/GLGC_exome_TC_release.txt", header = T)

# exposure
o <- fread("~/Bioinformatics/Src/R/cluster_results/outcome_without_HLL_region/outcome/celiac_removed/celiac_2011.txt", header = T, sep = " ")

# calc se
e$se <- mr.calculate.se(B = e$Z_OR, pval = e$P)
o$se <- mr.calculate.se(B = o$Z_OR, pval = o$P)

# pre process
e <- mr.pre.process(B = e$Beta, B.se = e$SE, pval = e$`P-value`, effect_allele = e$A1, other_allele = e$A2, SNP = e$rsid)

# harmonize
h <- mr.harmonize(By = o$Z_OR, Bx = e$beta, By.se = o$se, Bx.se = e$se, outcome.pval = o$P, exposure.pval = e$pval, outcome.effect_allele = o$effect_allele, exposure.effect_allele = e$effect_allele, exposure.other_allele = e$other_allele,
                  outcome.SNP = o$SNP, exposure.SNP = e$SNP)

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

p <- mr.plot(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se, iv = h$iv, iv.se = h$iv.se, ivw = ivw$ivw, egger = egger$egger, egger.i = egger$egger.i, ivw.Q = ivw.Q$ivw, chochran.Q = h$cochran.Q,
             egger.Q = egger.Q$egger, egger.i.Q = egger.Q$egger.i, egger.p.fdr = egger$egger.p, ivw.p.fdr = ivw$ivw.p, egger.i.p = egger$egger.i.p, exposure.name =  "Celiac 2011",
             outcome.name = "GLGC_exome_TC_release")
ggdraw(p)

ggsave(filename = "~/Desktop/test.png", plot = ggdraw(p))


ivw

ivw.fdr[length(ivw.fdr-1)]

ivw.fdr <- p.adjust(p = c(m.2011$beta.ivw.p, ivw$ivw.p), method = "fdr", n = length(m.2011$beta.ivw.p)+1)
egger.fdr <- p.adjust(p = c(m.2011$beta.egger.p, egger$egger.p), method = "fdr", n = length(m.2011$beta.egger.p)+1)


m.2011 <- read.table("~/Bioinformatics/Src/R/cluster_results/outcome_without_HLL_region/final_tables/master.table.methods.2011.txt", header = T)
m.2010 <- read.table("~/Bioinformatics/Src/R/cluster_results/outcome_without_HLL_region/final_tables/master.table.methods.2010.txt", header = T)

m.2011 <- m.2011[which(m.2011$beta.ivw.p.fdr < .05 | m.2011$beta.egger.p.fdr < .05), ]
m.2010 <- m.2010[which(m.2010$beta.ivw.p.fdr < .05 | m.2010$beta.egger.p.fdr < .05), ]



