# library
library(data.table)
library(mendelianRandomization)
library(optparse)

# options 0S4MD5VC
option_list = list(
  make_option(c("-t", "--trait"),
              help=" "),
  make_option(c("-y", "--year"),
              help=" "),
  make_option(c("-p", "--port"),
              help=" "),
  make_option(c("-a", "--address"),
              help=" "))

opt_parser    <- OptionParser(option_list=option_list, description="")
opt           <- parse_args(opt_parser)


# globals
trait_name <- opt$t
year <- opt$y
port <- opt$p
location <- opt$a
path.to.data <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/uk.biobank.gwas.2/files/"
path.to.outcome <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/celiac.without.HLA/celiac_"
path.to.ref <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/1000G/EUR_1000G_phase3_vcf.bim"
path.to.1000G <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/1000G/EUR_1000G_phase3_vcf"
path.to.output.directory <- "/groups/umcg-wijmenga/tmp04/umcg-mknigge/uk.biobank.gwas.2/bidirectionalMR/"
forward <- "celiac~exposure/"
reverse <- "exposure~celiac/"

# read data
celiac <- fread(paste0(path.to.outcome, year, ".txt"), header = T)
trait <- fread(paste0(path.to.data, trait_name), header = T)

# calc se
celiac$se <- mr.calculate.se(B = celiac$Z_OR, pval = celiac$P)$se

# check if effectsize or odds ratio
if(length(which(trait$beta < 0)) == 0){
  trait$beta <- log(trait$beta)
  trait$se <- mr.calculate.se(B = trait$beta, pval = trait$pval)$se
}



# ------------------------------- FORWARD ----------------------------------------------
outcome <- celiac
exposure <- trait

# pre process
exposure <- mr.pre.process(B = exposure$beta, B.se = exposure$se, pval = exposure$pval, effect_allele = exposure$other_allele, other_allele = exposure$effect_allele, SNP = exposure$rsid)

# clump
exposure <- mr.clump(data = exposure, refdat = path.to.1000G, verbose = T, SNPs.of.opposite.file = outcome$SNP)

# harmonize
h <- mr.harmonize(By = outcome$Z_OR, Bx = exposure$beta, By.se = outcome$se, Bx.se = exposure$se, outcome.pval = outcome$P, exposure.pval = exposure$pval, outcome.effect_allele = outcome$effect_allele,
                  exposure.effect_allele = exposure$effect_allele, exposure.other_allele = exposure$other_allele, outcome.SNP = outcome$SNP, exposure.SNP = exposure$SNP)
# do wald
iv <- mr.wald.ratio(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)

# add
h$iv <- iv$iv; h$iv.se <- iv$iv.se; h$iv.p <-  iv$iv.p

# cochran's Q
h$cochran.Q <- mr.cochran.Q.test(data = h, pval = .05)$cochran.Q

# ivw
ivw <- mr.inverse.variance.weighted.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)
ivw.Q <- mr.inverse.variance.weighted.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se, subset = h$cochran.Q)

# egger
egger <- mr.egger.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)
egger.Q <- mr.egger.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se, subset = h$cochran.Q)

# bind
if(is.data.frame(h) && nrow(h) != 0){

  # bind frame
  forward.method <- data.frame(ivw = ivw$ivw, ivw.se = ivw$ivw.se, ivw.p = ivw$ivw.p,
                               egger = egger$egger, egger.se = egger$egger.se, egger.p = egger$egger.p,
                               egger.i = egger$egger.i, egger.i.se = egger$egger.i.se, egger.i.p = egger$egger.i.p,
                               ivw.Q = ivw.Q$ivw, ivw.Q.se = ivw.Q$ivw.se, ivw.Q.p = ivw.Q$ivw.p,
                               egger.Q = egger.Q$egger, egger.Q.se = egger.Q$egger.se, egger.Q.p = egger.Q$egger.p,
                               egger.Q.i = egger.Q$egger.i, egger.Q.i.se = egger.Q$egger.i.se, egger.Q.i.p = egger.Q$egger.i.p,
                               phenotype = gsub(".txt", "", trait_name), nSNP = length(h$By))
  # bind SNPs
  h$phenotype <- gsub(".txt", "", trait_name)
  forward.SNPs <- h

  # write output
  write.table(x = forward.SNPs, file = paste0(path.to.output.directory, forward, year, "/", gsub(".txt", "", trait_name), "_SNPs.txt"), quote = T, row.names = F)
  write.table(x = forward.method, file = paste0(path.to.output.directory, forward, year, "/", gsub(".txt", "", trait_name), "_method.txt"), quote = T, row.names = F)
}

# ------------------------------- REVERSE ----------------------------------------------
outcome <- trait
exposure <- celiac

# pre process
exposure <- mr.pre.process(B = exposure$Z_OR, B.se = exposure$se, pval = exposure$P, effect_allele = exposure$effect_allele, other_allele = exposure$other_allele, SNP = exposure$SNP)

# clump
exposure <- mr.clump(data = exposure, refdat = path.to.1000G, verbose = T, SNPs.of.opposite.file = outcome$rsid)

# harmonize
h <- mr.harmonize(By = outcome$beta, Bx = exposure$beta, By.se = outcome$se, Bx.se = exposure$se, outcome.pval = outcome$pval, exposure.pval = exposure$pval, outcome.effect_allele = outcome$other_allele,
                  exposure.effect_allele = exposure$effect_allele, exposure.other_allele = exposure$other_allele, outcome.SNP = outcome$rsid, exposure.SNP = exposure$SNP)

# do wald
iv <- mr.wald.ratio(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)

# add
h$iv <- iv$iv; h$iv.se <- iv$iv.se; h$iv.p <-  iv$iv.p

# cochran's Q
h$cochran.Q <- mr.cochran.Q.test(data = h, pval = .05)$cochran.Q

# ivw
ivw <- mr.inverse.variance.weighted.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)
ivw.Q <- mr.inverse.variance.weighted.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se, subset = h$cochran.Q)

# egger
egger <- mr.egger.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se)
egger.Q <- mr.egger.method(By = h$By, Bx = h$Bx, By.se = h$By.se, Bx.se = h$Bx.se, subset = h$cochran.Q)

# bind
if(is.data.frame(h) && nrow(h) != 0){

  # bind frame
  reverse.method <- data.frame(ivw = ivw$ivw, ivw.se = ivw$ivw.se, ivw.p = ivw$ivw.p,
                               egger = egger$egger, egger.se = egger$egger.se, egger.p = egger$egger.p,
                               egger.i = egger$egger.i, egger.i.se = egger$egger.i.se, egger.i.p = egger$egger.i.p,
                               ivw.Q = ivw.Q$ivw, ivw.Q.se = ivw.Q$ivw.se, ivw.Q.p = ivw.Q$ivw.p,
                               egger.Q = egger.Q$egger, egger.Q.se = egger.Q$egger.se, egger.Q.p = egger.Q$egger.p,
                               egger.Q.i = egger.Q$egger.i, egger.Q.i.se = egger.Q$egger.i.se, egger.Q.i.p = egger.Q$egger.i.p,
                               phenotype = gsub(".txt", "", trait_name), nSNP = length(h$By))
  # bind SNPs
  h$phenotype <- gsub(".txt", "", trait_name)
  reverse.SNPs <- h

  # write output
  write.table(x = reverse.SNPs, file = paste0(path.to.output.directory, reverse, year, "/", gsub(".txt", "", trait_name), "_SNPs.txt"), quote = T, row.names = F)
  write.table(x = reverse.method, file = paste0(path.to.output.directory, reverse, year, "/", gsub(".txt", "", trait_name), "_method.txt"), quote = T, row.names = F)
}

# talk to socket, letting know that MR analysis is done
con <- socketConnection(host = location, port = port, blocking = F)
# close connection
close(con)




