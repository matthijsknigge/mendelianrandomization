#' Pascal (Pathway scoring algorithm)
#' @author Matthijs Knigge
#'
#' @description Pascal is a tool that allows gene and pathway-level analysis of GWAS association results without the need to access genotypic data. This is an executable that can not be exported. This must be downloaded and installed by the user. Pascal is available for download here: https://www2.unil.ch/cbg/index.php?title=Pascal. Individual SNPs are mapped to genes and their association p-values are combined into gene scores; second, genes are grouped into pathways and their gene scores are combined into pathway scores.
#'
#' @param SNPs vector of SNPs. Default is NULL.
#' @param pval vector p values associated with the SNPs. Default is NULL.
#' @param path.to.ref Give the path to where the gnu-zipped tped-files are stored.
#'
#' @keywords Pascal
#' @export
#' @examples
#' mr.Pascal()
#'
#' @return List with the following elements:
#'
mr.Pascal <- function(SNPs = NULL, pval = NULL, path.to.ref){
  require(data.table)
  # output folder from pascal
  Pascal.bin <- paste0(system.file(package="mendelianRandomization", "executables", "PASCAL"), "/")
  # Pascal output
  Pascal.out <- paste0(system.file(package="mendelianRandomization", "executables", "PASCAL"), "/output/")
  # tempdir for Pascal
  if(!dir.exists(paste0(system.file(package="mendelianRandomization", "executables"), "/PASCAL/temp"))){
    # check if tempdir exists, if not create
    dir.create(file.path(system.file(package="mendelianRandomization", "executables"), "/PASCAL/temp"))
  }
  # Pascal tempdir
  tempdir <- file.path(system.file(package="mendelianRandomization", "executables"), "PASCAL/temp")
  # Make textfile
  fn <- tempfile(tmpdir = tempdir)
  write.table(data.frame(SNP=SNPs, snpPvalCol=pval), file=fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep = "\t")
  # function for executing Pascal
  fun2 <- paste0(Pascal.bin, "./Pascal", " --pval ", fn, " --genescoring=sum", " --runpathway=on")
  # perform function
  system(fun2)
  # get tmp file name
  tmp <- unlist(strsplit(x = fn, split = "/")); tmp <- tmp[length(tmp)]
  # paths to data from Pascal
  pathway.path     <- paste0(Pascal.out, tmp, ".PathwaySet--msigBIOCARTA_KEGG_REACTOME--sum.txt")
  fusion.path      <- paste0(Pascal.out, tmp, ".sum.fusion.genescores.txt")
  gene.scores.path <- paste0(Pascal.out, tmp, ".sum.genescores.txt")
  SNPs.errors.path <- paste0(Pascal.out, tmp, ".sum.numSnpError.txt")
  # read data
  pathway     <<- fread(pathway.path)
  fusion      <<- fread(fusion.path)
  gene.scores <<- fread(gene.scores.path)
  SNPs.errors <<- fread(SNPs.errors.path)

  # unlink tmp files
  unlink(fn); unlink(pathway.path); unlink(fusion.path); unlink(gene.scores.path); unlink(SNPs.errors.path);
}
