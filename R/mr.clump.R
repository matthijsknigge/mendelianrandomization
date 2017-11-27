#' Perform clumping
#' @author Matthijs Knigge
#'
#'
#' @param data data.frame to be clumped. Must contain columns: pval, SNP
#' @param verbose boolean TRUE for clumping and using proxy SNP, FALSE for just removing SNPs in LD. Default: False
#' @param refdat path to reference data for clumping (EUR_1000G_phase3_vcf)
#' @param clump_kb window size, default = 1000
#' @param clump_r2 clumping r2 cutoff, default = .1
#' @param clump_p1 clumping sig level for index SNPs, default = 5*10^-8
#' @param clump_p2 clumping sig level for index SNPs, default = 5*10^-8
#' @param SNPs.of.opposite.file vector of SNPs of the opposite file for reference
#' @keywords clump
#' @export
#' @examples
#' mr.clump(data = data.frame, refdat = '/path/to/refdata/EUR_1000G_phase3_vcf', verbose = TRUE)
#'
#' @return data.frame with removed snps
mr.clump <- function(data, refdat, clump_kb = 1000, clump_r2 = .1, clump_p1 = 5*10^-8, clump_p2 = 5*10^-8, verbose = FALSE, SNPs.of.opposite.file){
  require(R.utils)
  # check size
  if(length(data$SNP) == 1){
    return(data)
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
  write.table(data.frame(SNP=data$SNP, P=data$pval), file=fn, row.names=F, col.names=T, quote=F)
  # call plink
  if(verbose == TRUE){
    fun2 <- paste0(plink_bin, "./plink", " --bfile ", refdat, " --clump ", fn, " --clump-p1 ", clump_p1, " --clump-p2 ", clump_p2, " --clump-r2 ", clump_r2, " --clump-kb ", clump_kb, " --clump-verbose", " --out ", fn)
  }
  if(verbose == FALSE){
    fun2 <- paste0(plink_bin, "./plink", " --bfile ", refdat, " --clump ", fn, " --clump-p1 ", clump_p1, " --clump-p2 ", clump_p2, " --clump-r2 ", clump_r2, " --clump-kb ", clump_kb, " --out ", fn)
  }
  # perform function
  system(fun2)

  if(file.exists(paste(fn, ".clumped", sep="")) == TRUE){
    # clump verbose
    if(verbose == TRUE){
      SNP <- mr.parse.verbose(file = paste(fn, ".clumped", sep=""), SNPs.list.1 = data$SNP, SNPs.list.2 = SNPs.of.opposite.file)
    }
    # clump
    if(verbose == FALSE){
      SNP <- read.table(paste(fn, ".clumped", sep=""), header=T)$SNP
    }
    # remove temp clump files
    unlink(paste0(fn, "*"))
    # remove SNPs in LD
    data <- data[data$SNP %in% SNP, ]
    # clear work space
    rm(snp.list, envir = .GlobalEnv); rm(proxy.snp, envir = .GlobalEnv)
    # return result
    return(data)
  } else {
    message("Nothing within Linkage Disequilibrium")
    return(data)
  }
}


#' Parse verbose, get proxy snps
#' @author Matthijs Knigge
#'
#' @param file file conatianing snps and pval to clump
#' @param SNPs list of snps to compare with proxy determination
#'
#' @return list of SNPs
mr.parse.verbose <- function(file, SNPs.list.1, SNPs.list.2){
  require(stringr); require(readr)
  # read entire file as single character
  file <- read_file(file)
  # pattern for finding chunks
  pattern <- regex("(SNP)([.\\S\\s]*?)(--+)"); message("chunking for pattern.....")
  # split verbose file in chunks: (SNP) [any character between] (--+)
  matches <- str_extract_all(string = file, pattern = pattern, simplify = TRUE)
  # store SNPs found
  snp.list <<- c()
  # store proxy SNPs found
  proxy.snp <<- data.frame()
  # for every chunk found
  message("evaluating chunks.....")
  for(match in matches){
    # check if the chunk is a proxy chunk
    if(grepl("INDEX", match) == TRUE){
      # pattern for proxy chunk
      pattern <- regex("(KB)([.\\S\\s]*?)(--+)")
      # get all proxy SNPs and there RSQ
      proxy.chunk <- str_extract_all(string = match, pattern = pattern, simplify = TRUE)
      # pattern for SNP chunk
      pattern <- regex("(rs\\d*)(.*)")
      # extract SNP and RSQ from proxy chunk
      snp.chunk <- str_extract_all(string = proxy.chunk, pattern = pattern, simplify = TRUE)
      # split snp chunk into frame
      for(snp in snp.chunk){
        SNPs <- unlist(strsplit(snp, "\\s+"))
        proxy.snp <<- rbind(proxy.snp, data.frame("SNP" = SNPs[1], "RSQ" = as.numeric(SNPs[3])))
      }
      # order proxy SNP data.frame
      proxy.snp <<- proxy.snp[order(proxy.snp$RSQ, decreasing = T), ]
      # remove all below 0.8
      if(length(which(proxy.snp$RSQ < .8)) > 0){
        proxy.snp <<- proxy.snp[-which(proxy.snp$RSQ < .8), ]
      }
      # determine wich SNPs to use as proxy
      int <- Reduce(intersect, list(SNPs.list.1, SNPs.list.2, proxy.snp$SNP))
      # apply overlap
      proxy.snp <<- proxy.snp[proxy.snp$SNP %in% int, ]
      # get max RSQ
      snp.list <<- c(snp.list, as.character(proxy.snp[which.max(proxy.snp$RSQ), ]$SNP))
      # empty the frame for next iteration
      proxy.snp <<- data.frame()
    }
    # if it is a chunk without proxy
    if(grepl("INDEX", match) == FALSE){
      # pattern for SNP
      pattern <- regex("(rs\\d*)")
      # get the SNP
      snp <- str_extract_all(string = match, pattern = pattern, simplify = TRUE)
      # match SNP with outcome and exposure
      int <- Reduce(intersect, list(SNPs.list.1, SNPs.list.2, snp))
      # apply overlap
      snp.list <<- c(snp.list, as.character(int))
    }
  }
  message("Done clumping verbose")
  # report
  for(snp.found in snp.list){
    message(paste0("found ", snp.found))
  }
  # return clumped verbose data
  return(snp.list)
}
