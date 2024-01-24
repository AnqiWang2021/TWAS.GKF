#' Match alleles between summary statistics and SNP information.
#'
#' This function performs allele flip by comparing input alleles with reference alleles.
#' It determines which alleles should be kept based on matching or complementarity with the reference.
#' The function also determines if alleles need to be flipped to match the reference.
#'
#' @param a1 Vector of first alleles from the dataset.
#' @param a2 Vector of alternative alleles from the dataset.
#' @param ref1 Vector of first reference alleles.
#' @param ref2 Vector of alternative reference alleles.
#' @return A list with two elements:
#' \itemize{
#'   \item \code{keep}: Logical vector indicating if the allele pair should be kept.
#'   \item \code{flip}: Logical vector indicating if the allele pair should be flipped to match the reference.
#' }
#' @export
allele.qc= function(a1,a2,ref1,ref2){

  ref = ref1
  flip = ref

  flip[ref=="T"]=="C"
  flip[ref=="A"]=="C"
  flip[ref=="C"]=="G"
  flip[ref=="A"]=="G"
  flip[ref=="C"]=="T"
  flip[ref=="C"]=="A"
  flip[ref=="G"]=="C"
  flip[ref=="G"]=="A"
  flip[ref=="T"]=="G"
  flip[ref=="T"]=="A"
  flip[ref=="G"]=="T"
  flip[ref=="A"]=="T"
  flip1 = flip

  ref = ref2
  flip = ref
  flip[ref=="T"]=="C"
  flip[ref=="A"]=="C"
  flip[ref=="C"]=="G"
  flip[ref=="A"]=="G"
  flip[ref=="C"]=="T"
  flip[ref=="C"]=="A"
  flip[ref=="G"]=="C"
  flip[ref=="G"]=="A"
  flip[ref=="T"]=="G"
  flip[ref=="T"]=="A"
  flip[ref=="G"]=="T"
  flip[ref=="A"]=="T"
  flip2 = flip

  snp = list()
  snp[["keep"]] = !(a1=="T"&a2=="C")|(a1=="A" & a2=="C")|(a1=="A"&a2=="G")|(a1=="C"&a2=="G")|(a1=="C"&a2=="T")|(a1=="C"& a2=="A")|(a1=="G"&a2=="A")|(a1=="G"&a2=="C")|(a1=="T"& a2=="G")|(a1=="T" & a2=="A")|(a1=="G"& a2=="T")|(a1=="A"& a2=="T")
  snp[["keep"]][a1!="T"&a1!="A"&a1!="C"&a1!="G"]=F
  snp[["keep"]][a2!="T"&a2!="A"&a2!="C"&a2!="G"]=F
  snp[["flip"]] = (a1==ref2&a2 ==ref1)|(a1==flip2&a2 ==flip1)

  return(snp)

}
