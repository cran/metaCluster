#' @title Calculation of GC content
#'
#' @description This function will calculate GC content from each sequence or contigs of a FASTA file.
#'
#' @param fasta_file
#'
#' @return GC content
#'
#' @examples
#'
#' @export

GC.content <- function(fasta_file){
  requireNamespace("seqinr")
  x <- fasta_file
  tt<-function(x){
    res<-GC(x)
    val=round(res,4)
    return(val)
  }

  f_res<-lapply(x,tt)
  s=data.frame(f_res)

  rownames(s) <- c("GC-content")

  w=t(s)
  return(w)
}
