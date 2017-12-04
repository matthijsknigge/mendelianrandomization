#' PLINK show-tag for finding proxy SNPs
#' @author Matthijs Knigge
#'
#'
#' @param data this must be the data structure given from mr.harmonize
#' @param refdat path to reference data for tagging (EUR_1000G_phase3_vcf)
#' @param originaldat path to the original data, this wil act as a reference for choosing which proxy SNPs to use. Should contain SNP, beta, se
#' @param popdat path to the population data.
#' @param tag_kb window size, default = 1000
#' @param tag_r2 clumping r2 cutoff, default = .1
#' @keywords show tag
#' @export
#' @examples
#' mr.show.tag(SNP = SNP, refdat = '/path/to/refdata/EUR_1000G_phase3_vcf', originaldat = 'path/to/original/gwas.txt', popdat = 'path/to/population/data.bim')
#'
#' @return List with the following: vector of SNPs that was given. vector replacement with corresponding replacement or indication that there is no replacement
mr.show.tag <-function(data, refdat, originaldat, popdat, tag_r2 = .8, tag_kb = 1000){
  require(R.utils); require(data.table); require(stringr); require(readr); library(progress)
  # read popdat
  popdat <- fread(popdat); colnames(popdat) <- c("chr_name", "SNP", "start", "pos", "A1", "A2")
  # check if there are any SNPs that are not in popdat
  int <- Reduce(intersect, list(data$SNP, popdat$SNP))
  # SNPs in popdat
  SNPs.in.pop <- data[data$SNP %in% int, ]
  # SNPs not in popdat
  SNPs.not.in.pop <- data[!data$SNP %in% int, ]
  # if everything is matched, return
  if(length(SNPs.not.in.pop$SNP) == 0){
    return(list(data = data))
  }
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
  write.table(data.frame(SNP=SNPs.not.in.pop$SNP), file=fn, row.names=F, col.names=T, quote=F)
  # call plink
  fun2 <- paste0(plink_bin, "./plink", " --bfile ", refdat, " --show-tags ", fn, " --tag-r2 ", tag_r2, " --tag-kb ", tag_kb, " --out ", fn, " --list-all")
  # perform function
  system(fun2)
  # read out from PLINK
  tags.list <- fread(paste0(fn, ".tags.list"))
  # unlink temp
  # unlink(paste0(fn, "*"))
  # read original GWAS
  originaldat <- fread(originaldat)
  # # function for lapply
  parse.tag <- function(tags){
    # update
    pb$tick()
    # parse all tags into list
    SNP.tag <- str_extract_all(string = tags, pattern = regex("rs\\d*"), simplify = TRUE)
    # intercept with population and original GWAS
    int <- Reduce(intersect, list(SNP.tag, popdat$SNP, originaldat$SNP))
    # return NA if there are no matches
    if(length(int) == 0){
      return(NA)
    }
    # if there are matches, return first element
    if(length(int) != 0){
      return(int[1])
    }
  }
  # filter data
  print("Parsing tags")
  pb <- progress_bar$new(
    format = "[:bar] :percent in :elapsed",
    total = length(tags.list$TAGS), clear = FALSE, width = 60)
  tags.list$TAGS <- lapply(tags.list$TAGS, FUN = function(x)  parse.tag(x))
  # remove NA
  tags.list <- tags.list[!which(is.na(tags.list$TAGS)), ]
  # remove from SNPs.not.in.pop that is not in tags.list
  SNPs.not.in.pop <- SNPs.not.in.pop[SNPs.not.in.pop$SNP %in% tags.list$SNP]
  # filter original dat
  originaldat <- originaldat[originaldat$SNP %in% tags.list$TAGS, ]
  # swap SNP, beta, se, and pval
  SNPs.not.in.pop[which(SNPs.not.in.pop$SNP == tags.list$SNP), c("SNP", "Bx", "Bx.se", "pval")] <- originaldat[which(originaldat$SNP == tags.list$TAGS), c("SNP", "beta", "se", "pval")]
  # bind result together
  data <- rbind(SNPs.in.pop, SNPs.not.in.pop)
  # annotate the proxy SNPs
  data$proxy <- 0
  # 1 means that it is proxy
  data[data$SNP %in% SNPs.not.in.pop$SNP]$proxy <- 1
  # return data
  return(list(data = data))

}


