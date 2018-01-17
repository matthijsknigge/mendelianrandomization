#' Pascal (Pathway scoring algorithm)
#' @author Matthijs Knigge
#'
#' @description Pascal is a tool that allows gene and pathway-level analysis of GWAS association results without the need to access genotypic data. This is an executable that can not be exported. This must be downloaded and installed by the user. Pascal is available for download here: https://www2.unil.ch/cbg/index.php?title=Pascal. Individual SNPs are mapped to genes and their association p-values are combined into gene scores, and genes are grouped into pathways and their gene scores are combined into pathway scores.
#'
#' @param SNPs vector of SNPs.
#' @param pval vector p values associated with the SNPs.
#' @param Pascal.root path to the root folder of Pascal
#' @param cochran.Q if in the data set a subset of SNPs is removed, this vector can be used to look for pathways that distinguish both subsets of SNPs. Default is NULL
#'
#' @keywords Pascal
#' @export
#' @examples
#' mr.Pascal()
#'
#' @return List with the following elements:
#'         pathways: significant pathways found by Pascal on SNPs that are not removed by Cochran's Q test
#'         pathways.Q: significant pathways found by Pascal on SNPs that are removed by Cochran's Q test
#'
mr.Pascal <- function(SNPs, pval, Pascal.root, cochran.Q = NULL){
  require(data.table)
  # check cochran's Q
  if(sum(cochran.Q) == 0){
    cochran.Q <- NULL
  }
  # output folder from pascal
  Pascal.bin <- Pascal.root
  # Pascal output
  Pascal.out <- paste0(Pascal.root, "/output/")
  # tempdir for Pascal
  if(!dir.exists(paste0(Pascal.root, "/temp"))){
    # check if tempdir exists, if not create
    dir.create(file.path(paste0(Pascal.root, "/temp")))
  }
  # Pascal tempdir
  tempdir <- paste0(Pascal.root, "/temp")
  # get pathways for SNPs if Cochran's Q test was not applied
  if(is.null(cochran.Q)){
    # call Pascal and get pathways
    pathways <- mr.perform.Pascal(Pascal.bin = Pascal.bin, Pascal.out = Pascal.out, tempdir = tempdir, SNPs = SNPs, pval = pval)$pathways
    # force equal lengths of vectors
    paths <- list(pathways = pathways, pathways.Q = NA)
    # get max length
    n <- max(lengths(paths))
    # introduce NA
    paths <- lapply(paths, `length<-`, n)
    # return pathways
    return(list(pathways = paths$pathways, pathways.Q = paths$pathways.Q))
  }
  # get pathways for SNPs if cochran's Q test was applied
  if(!is.null(cochran.Q)){
    # call Pascal subset of SNPs that was not removed by Cochran's Q test
    pathways <- mr.perform.Pascal(Pascal.bin = Pascal.bin, Pascal.out = Pascal.out, tempdir = tempdir, SNPs = SNPs[which(cochran.Q == 1)], pval = pval[which(cochran.Q == 1)])$pathways
    # call Pascal subset of SNPs that was not removed by Cochran's Q test
    pathways.Q <- mr.perform.Pascal(Pascal.bin = Pascal.bin, Pascal.out = Pascal.out, tempdir = tempdir, SNPs = SNPs[which(cochran.Q == 0)], pval = pval[which(cochran.Q == 0)])$pathways
    # force equal lengths of vectors
    paths <- list(pathways = pathways, pathways.Q = pathways.Q)
    # get max length
    n <- max(lengths(paths))
    # introduce NA
    paths <- lapply(paths, `length<-`, n)
    # return pathways
    return(list(pathways = paths$pathways, pathways.Q = paths$pathways.Q))
  }
}


mr.perform.Pascal <- function(Pascal.bin, Pascal.out, tempdir, SNPs, pval){
  # Make textfile
  fn <- tempfile(tmpdir = tempdir)
  # output table
  write.table(data.frame(SNP=SNPs, snpPvalCol=pval), file=fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep = "\t")
  # function for executing Pascal
  fun2 <- paste0(Pascal.bin, "./Pascal", " --pval ", fn, " --genescoring=sum", " --runpathway=on")
  # perform function
  print("Performing pathway analysis......"); system(fun2, ignore.stdout = T, ignore.stderr = T); print("done")
  # get tmp file name
  tmp <- unlist(strsplit(x = fn, split = "/")); tmp <- tmp[length(tmp)]
  # paths to output data from Pascal
  pathway.path     <- paste0(Pascal.out, tmp, ".PathwaySet--msigBIOCARTA_KEGG_REACTOME--sum.txt")
  fusion.path      <- paste0(Pascal.out, tmp, ".sum.fusion.genescores.txt")
  gene.scores.path <- paste0(Pascal.out, tmp, ".sum.genescores.txt")
  SNPs.errors.path <- paste0(Pascal.out, tmp, ".sum.numSnpError.txt")
  # read data
  pathway <- fread(pathway.path)
  # get significant path ways
  pathways <- pathway[which(pathway$chi2Pvalue < .05), ]$Name
  # unlink tmp files
  unlink(fn); unlink(pathway.path); unlink(fusion.path); unlink(gene.scores.path); unlink(SNPs.errors.path);
  # return significant pathways
  return(list(pathways = pathways))
}
