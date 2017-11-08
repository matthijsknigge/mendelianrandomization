# libraries and functions
source("MRfunctions.R")
require(MRInstruments)
library(data.table)
library(R.utils)
library(optparse)


# options
option_list = list(
  make_option(c("-f", "--file"), type="character",
              help="Input file containing paths to summary statistics."))
opt_parser    <- OptionParser(option_list=option_list, description="")
opt           <- parse_args(opt_parser)

# global
temdir.clump <- "/groups/umcg-wijmenga/tmp03/users/umcg-mknigge/hemoglobin/temp_clump_2011"
refdata <- "/groups/umcg-wijmenga/tmp03/users/umcg-mknigge/hemoglobin/1000G/EUR_1000G_phase3_vcf"
plink_bin_location <- "/groups/umcg-wijmenga/tmp03/users/umcg-mknigge/hemoglobin/PLINK/"
output_directory_data  <- "/groups/umcg-wijmenga/tmp03/users/umcg-mknigge/hemoglobin/results_2011/"
bim_file <- "/groups/umcg-wijmenga/tmp03/users/umcg-mknigge/hemoglobin/1000G/EUR_1000G_phase3_vcf.bim"
outcome <-  "/groups/umcg-wijmenga/tmp03/users/umcg-mknigge/hemoglobin/outcome/celiac_2011.txt"
data <- "/groups/umcg-wijmenga/tmp03/users/umcg-mknigge/hemoglobin/data/"

# outcome data
outcome <- fread(outcome, header = T); colnames(outcome) <- c("SNP", "effect_allele", "beta", "pval")

# calculate se for outcome
outcome$se <-sqrt(outcome$beta^2/(qchisq(outcome$pval, 1, lower.tail = FALSE)))

# global frames
all.methods <<- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c("phenotype", "beta.ivw", "beta.ivw.se", "beta.ivw.p", "beta.egger", "beta.egger.se", "beta.egger.p", "beta.egger.i", "beta.egger.i.se", "beta.egger.i.p"))
all.calcuations <<- setNames(data.frame(matrix(ncol = 9, nrow = 0)), c("phenotype", "SNP" ,"beta.exposure","se.exposure", "beta.outcome", "se.outcome", "beta.iv", "beta.iv.se", "beta.iv.p"))

# exposure data
exposure <- fread(paste0(data, opt$file), header = T); 
exposure$VARIANT <- NULL; exposure$CHR <- NULL; exposure$BP <- NULL; exposure$ALT_MINOR <- NULL; exposure$DIRECTION <- NULL; exposure$MLOG10P <- NULL; exposure$ALT_FREQ <- NULL; exposure$MA_FREQ <- NULL

colnames(exposure) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval")

# save outcome
celiac <- outcome

# pre-process
exposure <- mr.pre.process(exposure)

# clump data
exposure <- mr.clump.data(data = exposure, snp.exp = exposure$SNP, snp.out = celiac$SNP, refdata = refdata,
                          temdir.clump = temdir.clump, plink_bin_location = plink_bin_location)

# perform calculations
calculations <- mr.do.calculation(celiac, exposure, gsub(".tsv", "", opt$file))

plot(mr.calculations$beta.exposure, mr.calculations$beta.outcome)


# save calculations global
all.methods <<- rbind(all.methods, calculations$mr.methods)
all.calcuations <<- rbind(all.calcuations, calculations$mr.calculations)

# # save data
write.table(file = paste0(output_directory_data, gsub(".tsv", "", opt$file), "_methods.txt"), all.methods, quote = T, row.names = F)
write.table(file = paste0(output_directory_data, gsub(".tsv", "", opt$file), "_calculations.txt"), all.calcuations, quote = T, row.names = F)









