# libraries
library(data.table)
library(mendelianRandomization)

# data
celiac.2011.exposures <- fread("~/Bioinformatics/mr.results/all.combined.corrected/celiac_2011~exposures.txt")
celiac.2011.exposures.s <- fread("~/Bioinformatics/mr.results/all.combined.corrected/celiac_2011~exposures_SNPs.txt")


m <- celiac.2011.exposures[which(celiac.2011.exposures$phenotype == "I10-I15 Hypertensive diseases"), ]
s <- celiac.2011.exposures.s[which(celiac.2011.exposures.s$phenotype == "I10-I15 Hypertensive diseases"), ]

paths <- mr.Pascal(SNPs = s$SNP, pval = s$pval, Pascal.root = "~/Bioinformatics/PASCAL", cochran.Q = s$cochran.Q, p.threshold = .2)



m1 <- celiac.2011.exposures[which(celiac.2011.exposures$phenotype == "E00-E07 Disorders of thyroid gland"), ]
s1 <- celiac.2011.exposures.s[which(celiac.2011.exposures.s$phenotype == "E00-E07 Disorders of thyroid gland"), ]

paths1 <- mr.Pascal(SNPs = s1$SNP, pval = s1$pval, Pascal.root = "~/Bioinformatics/PASCAL", cochran.Q = s1$cochran.Q, p.threshold = .2, supress.output = F, up = 50000, down = 50000)

paths2 <- mr.Pascal(SNPs = s1$SNP, pval = s1$pval, Pascal.root = "~/Bioinformatics/PASCAL", cochran.Q = s1$cochran.Q, p.threshold = .05, supress.output = F, up = 60000, down = 60000)

# pascal default = 50,000 bp

# 250000 plink default for clumping is 250.000 bp
