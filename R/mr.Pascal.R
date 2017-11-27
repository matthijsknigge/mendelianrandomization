#' Pascal (Pathway scoring algorithm)
#' @author Matthijs Knigge
#'
#' @description Pascal is a tool that allows gene and pathway-level analysis of GWAS association results without the need to access genotypic data. This is an executable that can not be exported. This must be downloaded and installed by the user. Pascal is available for download here: https://www2.unil.ch/cbg/index.php?title=Pascal. Individual SNPs are mapped to genes and their association p-values are combined into gene scores; second, genes are grouped into pathways and their gene scores are combined into pathway scores.
#'
#' @param SNPs vector of SNPs. Default is NULL.
#' @param pval vector p values associated with the SNPs. Default is NULL.
#' @param custom path to the root folder of Pascal
#'
#' @keywords Pascal
#' @export
#' @examples
#' mr.Pascal()
#'
#' @return List with the following elements:
#'
mr.Pascal <- function(SNPs, pval, Pascal.root){
  require(data.table)
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
  # Make textfile
  fn <- tempfile(tmpdir = tempdir)
  write.table(data.frame(SNP=SNPs, snpPvalCol=pval), file=fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep = "\t")
  # function for executing Pascal
  fun2 <- paste0(Pascal.bin, "/./Pascal", " --pval ", fn, " --genescoring=sum", " --runpathway=on")
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
  pathway <- fread(pathway.path)
  # get significant path ways
  pathways <- pathway[which(pathway$chi2Pvalue < .05), ]$Name
  # unlink tmp files
  unlink(fn); unlink(pathway.path); unlink(fusion.path); unlink(gene.scores.path); unlink(SNPs.errors.path);
  # return significant pathways
  return(list(pathways = pathways))
}
