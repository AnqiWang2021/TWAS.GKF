#' Perform Transcriptome-Wide Association Study (TWAS)
#'
#' This function conducts a TWAS using genomic data to assess its association with a phenotype.
#' It applies elastic net regularization to build a predictive model.
#' Cross-validation is used to evaluate the model, and the function outputs the weights and
#' predictions for the best model based on performance metrics.
#'
#' @param genomic_data Matrix of genomic data where rows are samples and columns are genetic features.
#' @param phenotype_data Numeric vector of phenotype data where each element corresponds to a sample.
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{wgt.matrix}}{Matrix of effect weights from the best model, zeroed if performance criteria are not met.}
#'     \item{\code{pred.wgt}}{Matrix of predicted values based on the best model, zeroed if performance criteria are not met.}
#'   }
#' @importFrom glmnet cv.glmnet
#' @importFrom stats cor.test
#' @export
#'
twas <- function(genomic_data, phenotype_data){
  ##### Elastic Net
  enet = cv.glmnet(x = genomic_data,y = phenotype_data,nfolds=5,keep=T,alpha = 0.5,intercept=F)
  ##### Cross validation analysis
  eff.wgt = as.numeric(coefficients(enet,s=enet$lambda.min))[-1]
  pred.wgt = as.matrix(enet$fit.preval[,which(enet$lambda==enet$lambda.min)])
  ##### Performance
  performance = cor.test(phenotype_data,pred.wgt)
  ##### Fitted_value of performance
  fit_r = as.numeric(performance$estimate)
  fit_p = as.numeric(performance$p.value)

  fit_r = ifelse(is.na(fit_r),0,fit_r)
  fit_p = ifelse(is.na(fit_p),1,fit_p)

  ##### Validate the performance
  if (fit_r>0.1 & fit_p<0.05)
  {
    eff.wgt = as.matrix(eff.wgt)
    pred.wgt = as.matrix(pred.wgt)
  }else
  {
    eff.wgt = as.matrix(rep(0,length(eff.wgt)))
    pred.wgt = as.matrix(rep(0,nrow(pred.wgt)))
  }

  return(list(wgt.matrix = eff.wgt, pred.wgt = pred.wgt))
}

#' Transcriptome-Wide Association Study with GhostKnockoff Filter
#'
#' This function performs a TWAS using the GhostKnockoff Filter (GKF) to control the false discovery rate.
#' It calculates the variance and effect size estimate of gene expression based on genetic variants,
#' and then determines the significance of each gene's association with the phenotype using z-scores and knockoff statistics.
#'
#' @param weight A numeric vector of SNP weights for predicting gene expression.
#' @param LD The linkage disequilibrium matrix for the SNPs.
#' @param sumstats A data frame or matrix containing summary statistics from GWAS, with the columns "beta" and "se", including beta coefficients and standard errors.
#' @param gene_matrix A matrix of gene expression data.
#' @param n.study A vactor with the number of studies or samples.
#' @param M The number of knockoff variables to be generated, default is 5.
#' @param FDR_level The desired false discovery rate level, default is 0.05.
#' @return A list containing:
#'   \itemize{
#'     \item \code{var}: Variance of the gene expression regulated by the SNPs.
#'     \item \code{est}: Effect size estimate.
#'     \item \code{z}: Z-score for the association.
#'     \item \code{qvalue}: qß-values from the knockoff filter.
#'     \item \code{gene_select}: Indices of the selected genes passing the FDR threshold.
#'   }
#' @importFrom stats cor
#' @importFrom corpcor cor.shrink
#' @export
TWAS_GKF <- function(weight, LD, sumstats, gene_matrix, n.study, M = 5, FDR_level = 0.05){
    #### Effect size variance
    var_GReX = as.numeric(t(weight)%*%as.matrix(LD)%*%weight)
    #### Effect size estimate
    effect_size = sum(weight*sumstats$beta*diag(as.matrix(LD))/var_GReX)
    #### Z-score
    zscore = sum(weight*sqrt(diag(as.matrix(LD)))*sumstats$beta/sqrt(var_GReX)/sumstats0$se,na.rm=T)

    cor = matrix(as.numeric(corpcor::cor.shrink(gene_matrix,verbose=F)),nrow = ncol(gene_matrix))
    prelim = GhostKnockoff.prelim(cor,M=M)
    filter = GhostKnockoff.fit(as.matrix(zscore),n.study,fit.prelim=prelim)

    ff = GhostKnockoff.filter(filter$T_0,filter$T_k)
    gene_select = which(ff$q <= FDR_level)
    return(var = var_GReX, est = effect_size, z = zscore, qvalue = ff$q, gene_select = gene_select)
}
