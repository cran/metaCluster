#' @title Oligonucleotide Frequency
#'
#' @description This function will calculate oligonucleotide frequency of each sequence or contig from a FASTA file.
#'
#' @param fasta_file
#'
#' @param f
#'
#' @return Oligonucleotide frequency
#'
#' @examples
#'
#' @export
  oligo.freq <- function(fasta_file,f){
      requireNamespace("Biostrings")

  x<- readDNAStringSet(fasta_file)
  y <- oligonucleotideFrequency(x,width = f)
  z <- data.frame(y)
  rownames(z) <- names(x)

  return(z)
}
