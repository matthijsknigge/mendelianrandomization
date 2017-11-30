#' PLINK show-tag for finding proxy SNPs
#' @author Matthijs Knigge
#'
#'
#' @param SNP vector containing SNPs to be proxied
#' @param refdat path to reference data for tagging (EUR_1000G_phase3_vcf)
#' @param originaldat path to the original data, this wil act as a reference for choosing which proxy SNPs to use.
#' @param popdat path to the population data.
#' @param tag_kb window size, default = 1000
#' @param tag_r2 clumping r2 cutoff, default = .1
#' @keywords show tag
#' @export
#' @examples
#' mr.show.tag(SNP = SNP, refdat = '/path/to/refdata/EUR_1000G_phase3_vcf')
#'
#' @return List with the following: vector of SNPs that was given. vector replacement with corresponding replacement or indication that there is no replacement
mr.show.tag <-function(SNP, refdat, originaldat, popdat, tag_r2 = .8, tag_kb = 1000){
  require(R.utils); require(data.table); require(stringr); require(readr)
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
  write.table(data.frame(SNP=SNP), file=fn, row.names=F, col.names=T, quote=F)
  # call plink
  fun2 <- paste0(plink_bin, "./plink", " --bfile ", refdat, " --show-tags ", fn, " --tag-r2 ", tag_r2, " --tag-kb ", tag_kb, " --out ", fn, " --list-all")
  # perform function
  system(fun2)
  # read out from PLINK
  tags.list <- fread(paste0(fn, ".tags.list"))
  # unlink temp
  unlink(paste0(fn, "*"))
  # filter data
  tags.list$TAGS <- lapply(tags.list$TAGS, function(x) if(x != "NONE"){
    str_extract_all(string = x, pattern = regex("rs\\d*"), simplify = TRUE)[1]
  })
  # return data
  return(list(SNP = tags.list$SNP, replacement = tags.list$TAGS))

}


