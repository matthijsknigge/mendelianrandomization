#' Allelic Scoring
#' @author Matthijs Knigge
#'
#' @param SNP Vector of SNPs
#' @param B Vector of genetic effects
#' @param alleles Vector of effect alleles associated with genetic effects
#' @param blife absolute path to the bfile
#'
#' @keywords se
#' @export
#' @examples
#' mr.allelic.scoring()
#'
#' @return List with the following elements:
#'         se: standard error of genetic effect
mr.allelic.scoring <- function(SNP, B, alleles, bfile){
  # check if tempt exists
  if(!dir.exists(paste0(system.file(package="mendelianRandomization", "executables"), "/PLINK/temp_clump"))){
    dir.create(file.path(system.file(package="mendelianRandomization", "executables"), "/PLINK/temp_clump"))
  }
  # plink_bin
  plink_bin <- paste0(system.file(package="mendelianRandomization", "executables"), "/PLINK/")
  # temp folder for clumping
  tempdir <-   paste0(system.file(package="mendelianRandomization", "executables"), "/PLINK/temp_clump")
  # Make textfile
  fn <- tempfile(tmpdir = tempdir)
  write.table(data.frame(SNP=SNP, B=B, alleles=alleles), file=fn, row.names=F, col.names=T, quote=F)
  # call plink
  fun2 <- paste0(plink_bin, "./plink", " --bfile ", bfile, " --score ", fn, " 1 3 2 ", " header ", " --out ", fn)
  # execute function
  system(fun2)
  # read output
  profile <- fread(paste0(fn, ".profile"))
  # remove
  unlink(paste0(fn, "*"))
  # change phenotypes 1 to 0
  profile[which(profile$PHENO == 1), ]$PHENO <- 0
  # change phenotypes 2 to 1
  profile[which(profile$PHENO == 2), ]$PHENO <- 1
  # logistic regression
  logit <- glm(profile$PHENO ~ profile$SCORE, family = binomial(link='logit'))
  # summary
  smmry <- summary(logit)
  print(smmry)
  # get coeff
  estimate <- smmry$coefficients[1,1]
  # get pval
  pval <- smmry$coefficients[1,4]
  # return
  return(list(estimate = estimate, pval = pval))
}
