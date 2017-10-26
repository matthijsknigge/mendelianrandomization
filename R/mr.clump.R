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
#' @keywords clump
#' @export
#' @examples
#' mr.clump(data = data.frame, refdat = '/path/to/refdata/EUR_1000G_phase3_vcf', verbose = TRUE)
#'
#' @return data.frame with removed snps
mr.clump <- function(data, refdat, clump_kb = 1000, clump_r2 = .1, clump_p1 = 5*10^-8, clump_p2 = 5*10^-8, verbose = FALSE){
  # plink_bin
  plink_bin <- "exe/"
  # temp folder for clumping
  tempdir <-   "exe/temp_clump"
  # Make textfile
  fn <- tempfile(tmpdir = tempdir)
  write.table(data.frame(SNP=data$SNP, P=data$pval), file=fn, row=F, col=T, qu=F)
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
      SNP <- mr.parse.verbose(file = paste(fn, ".clumped", sep=""), SNPs = data$SNP)
    }
    # clump
    if(verbose == FALSE){
      SNP <- read.table(paste(fn, ".clumped", sep=""), header=T)$SNP
    }
    unlink(paste0(tempdir, fn, ".*"))

    # intersect between data sets
    int <- Reduce(intersect, list(data$SNP, SNP))
    # remove SNPs in LD
    data <- data[data$SNP %in% int, ]
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
mr.parse.verbose <- function(file, SNPs){
  snp <<- c()
  previous.line <<- ""
  track.snp <<- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("SNP", "r_2"))
  RSQ <<- 0
  count <<- 0
  snp.proxy.added <<- ""
  previous.RSQ <<- 0
  proxy.added <<- FALSE
  clump = readLines(file)
  for (i in 1:length(clump)){
    # add to count
    count <<- count + 1
    #------------------- remove trailing spaces and characters---------------
    clump[i] <- substr(clump[i], 70, nchar(clump[i]))
    if (clump[i] != ""){
      x <- unlist(strsplit(clump[i], "\\s+"))
      ws <- which(x == "")
      if(length(ws) != 0){
        x <- x[-which(x == "")]
        #------------------- remove trailing spaces and characters---------------
      }
      print(x)
      current.line <<- x
      if(grepl("rs", current.line[1]) == TRUE & previous.line[1] == "SNP"){
        snp <<- c(snp, current.line[1])
      }
      if(i == length(clump)){
        if(grepl("rs", current.line[1]) == TRUE & previous.line[1] == "SNP"){
          snp <<- c(snp, current.line[1])
        }
      }

      if(current.line[1] == "KB" & grepl("rs", previous.line[1]) == TRUE){
        snp <<- snp[-which(snp == previous.line[1])]
      }

      if(grepl("rs", previous.line[1]) == TRUE & grepl("rs", current.line[1]) == TRUE){

        if(current.line[1] %in% SNPs == TRUE){
          if(as.numeric(current.line[3]) > 0.8){
            # message(paste0("--------------------------------------", current.line[1], " ", current.line[3]))

            # add to frame
            track.snp <<- rbind(track.snp, data.frame("SNP" = current.line[1],
                                                      "r_2" = as.numeric(current.line[3])))
          }
        }
      }

      if(grepl("rs", previous.line[1]) == TRUE & current.line[1] == "SNP"){
        if(length(track.snp$SNP) > 0){
          message(paste0("--------------------------------------", track.snp[which.max(track.snp$r_2),]$SNP, " ", track.snp[which.max(track.snp$r_2),]$r_2))
          snp <<- c(snp, as.character(track.snp[which.max(track.snp$r_2),]$SNP))
          track.snp$SNP <<- NULL; track.snp$r_2 <<- NULL
        }
      }
      if(grepl("rs", previous.line[1]) == TRUE & current.line[1] == "KB"){

      }
      previous.line <<- x
    }
  }
  return(snp)
}
