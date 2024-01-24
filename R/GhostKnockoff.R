#' Preliminary Computations for Ghost Knockoff Variables
#'
#' This function performs preliminary computations needed to generate Ghost Knockoff variables.
#' It takes a covariance matrix of genes, the number of knockoff variables to generate,
#' the method for solving the semidefinite programming problem, and an optional size limit
#' for the problem. It returns matrices necessary for the subsequent steps in the Ghost
#' Knockoff generation process.
#'
#' @param cor.G A covariance matrix representing the relationships between genes.
#' @param M The number of knockoff variables to generate, default is 5.
#' @param method The method used for solving the optimization problem, either 'sdp' or 'asdp', default is 'asdp'.
#' @param max.size The maximum size for the semidefinite programming problem, default is 500.
#' @return A list containing the P.each matrix, the V.left matrix (left null space), the Normal_50Studies matrix,
#' and permute.index, which are used in the knockoff variable generation process.
#' @importFrom stats rnorm
#' @import Matrix
#' @export
GhostKnockoff.prelim<-function(cor.G, M=5, method='asdp', max.size=500){
  #method = 'asdp'
  temp.index<-1:nrow(cor.G)
  n.G<-nrow(cor.G)
  #Permutation test for constant variants in the reference panel
  permute.index<-rep(0,n.G)
  permute.index[-temp.index]<-1

  Normal_50Studies<-matrix(rnorm(n.G*M*50),n.G*M,50)
  P.each<-matrix(0,n.G,n.G)
  if(length(temp.index)!=0){
    Sigma<-cor.G[temp.index,temp.index,drop=F]
    SigmaInv<-solve(Sigma,verbose=F)#solve(Sigma)#invcov.shrink(Sigma,verbose=F)

    if(method=='sdp'){temp.s<-create.solve_sdp_M(Sigma,M=M)}
    if(method=='asdp'){temp.s<-create.solve_asdp_M(Sigma,M=M,max.size=max.size)}
    s<-temp.s
    diag_s<-diag(s,length(s))

    if(sum(s)==0){
      V.left<-matrix(0,length(temp.index)*M,length(temp.index)*M)
    }else{
      Sigma_k<-2*diag_s - s*t(s*SigmaInv)
      V.each<-Matrix(forceSymmetric(Sigma_k-diag_s))

      #random part of knockoff
      V<-matrix(1,M,M)%x%V.each
      diag(V)<-diag(V)+rep(s,M)
      V.left<-try(t(chol(V)),silent=T)
      if(class(V.left)=="try-error"){
        svd.fit<-svd(V)
        u<-svd.fit$u
        svd.fit$d[is.na(svd.fit$d)]<-0
        cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
        n.svd<-which(cump>=0.999)[1]
        if(is.na(n.svd)){n.svd<-nrow(V)}
        svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
        V.left<-t(sqrt(svd.fit$d[svd.index])*t(u[,svd.index,drop=F]))
      }
    }
    P.each[temp.index,temp.index]<-diag(1,length(s))-s*SigmaInv
    V.index<-rep(temp.index,M)+rep(0:(M-1),each=length(temp.index))*n.G
    Normal_50Studies[V.index,]<-as.matrix(V.left%*%matrix(rnorm(ncol(V.left)*50),ncol(V.left),50))

    #Permutation test for tightly linked variants
    permute.index[temp.index[s==0]]<-1
  }
  permute.V.index<-rep(permute.index,M)
  P.each[permute.index==1,]<-0
  Normal_50Studies[permute.V.index==1,]<-matrix(rnorm(sum(permute.index)*M*50),sum(permute.index)*M,50)

  return(list(P.each=as.matrix(P.each), V.left=V.left, Normal_50Studies=as.matrix(Normal_50Studies), permute.index=permute.index, M=M))
}

#' Fit Ghost Knockoff Model
#'
#' @param Zscore_0 Matrix of Z-scores for the original variables.
#' @param n.study Vector of the number of studies for each variable.
#' @param fit.prelim List of preliminary computations from `GhostKnockoff.prelim`.
#' @param gamma Adjustment parameter for the knockoff generation, default is 1.
#' @param weight.study Vector of weights for each study; if NULL, weights are computed as the square root
#' of the number of studies divided by the square root of the sum of all studies.
#' @return A list containing the weighted Z-scores for the original variables (`GK.Zscore_0`),
#' the knockoff variables (`GK.Zscore_k`), and the test statistics for the original (`T_0`)
#' and knockoff variables (`T_k`).
#' @import stats
#' @export
GhostKnockoff.fit<-function(Zscore_0, n.study, fit.prelim, gamma=1, weight.study=NULL){
  if(length(weight.study)==0){weight.study<-sqrt(n.study)/sqrt(sum(n.study))}
  #to account for study specific variants
  W.missing<-apply(!is.na(Zscore_0),2,as.numeric)
  W<-t(weight.study*t(W.missing))

  Zscore_0[is.na(Zscore_0)]<-0
  M<-fit.prelim$M
  n.G<-nrow(Zscore_0)
  P.each<-fit.prelim$P.each
  Normal_50Studies<-fit.prelim$Normal_50Studies

  for (i in 1:ncol(Zscore_0)){
    Normal_k<-matrix(Normal_50Studies[,i],nrow=n.G)
    Zscore_i0<-Zscore_0[,i,drop=F]
    Zscore_ik<-as.vector(P.each%*%Zscore_0[,i,drop=F])+gamma*Normal_k
    if(i==1){
      GK.Zscore_0<-W[,i]*Zscore_i0
      GK.Zscore_k<-W[,i]*Zscore_ik
    }else{
      GK.Zscore_0<-GK.Zscore_0+W[,i]*Zscore_i0
      GK.Zscore_k<-GK.Zscore_k+W[,i]*Zscore_ik
    }
  }

  T_0<-(GK.Zscore_0)^2
  T_k<-(GK.Zscore_k)^2

  return(list(GK.Zscore_0=GK.Zscore_0,GK.Zscore_k=GK.Zscore_k,T_0=T_0,T_k=T_k))
}

#' Calculate Study-Specific Correlation Matrix
#'
#' @import stats
#' @export
GhostKnockoff.GetCorStudy<-function(Zscore_0, fit.prelim){
  temp<-Zscore_0
  temp[is.na(temp)]<-0
  PZscore_0<-fit.prelim$P.each%*%temp

  cor.Z<-Zscore_0-PZscore_0
  cor.Z[abs(temp)>1.96]<-NA
  cor.study<-cor(cor.Z,use='pairwise.complete.obs')

  return(cor.study)
}

#' Meta Analysis for Preliminary Ghost Knockoff Computations
#' @export
GhostKnockoff.prelim.Meta<-function(cor.study, n.study){

  sigma <- cor.study
  k=nrow(sigma)
  A <- diag(sqrt(n.study))
  B <- diag(1,k)

  w <- Variable(k)
  objective <- Minimize(quad_form(w,cor.study))
  constraint1 <- t(as.matrix(sqrt(n.study))) %*% w -1 == 0
  constraint2 <- B %*% w >=0
  problem <- Problem(objective, constraints = list(constraint1, constraint2))
  result <- solve(problem)
  w_opt <- as.vector(result$getValue(w))
  w_opt<-w_opt/sqrt(result$value)

  N.ratio<-sum(w_opt^2)/sum(w_opt*t(w_opt*cor.study))
  gamma<-sqrt(1+1/N.ratio-N.ratio)

  return(list(w_opt=w_opt,gamma=gamma))
}

#' @import stats
#' @export
GhostKnockoff.filter<-function (T_0,T_k){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  M<-ncol(T_k);Rej.Bound<-10000

  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0

  which.max.alt<-function(x){
    temp.index<-which(x==max(x))
    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
  }
  kappa<-apply(T.temp,1,which.max.alt)-1

  Get.OtherMedian<-function(x){median(x[-which.max(x)])}
  tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)


  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));index_bound<-max(which(tau[b]>0))
  for(i in 1:length(b)){
    temp.index<-i:min(length(b),Rej.Bound,index_bound)
    if(length(temp.index)==0){next}
    q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  q[q>1]<-1
  return(list(kappa=kappa,tau=tau,q=q))
}

#' @importFrom Rdsdp dsdp
#' @importFrom Matrix Diagonal cov2cor
#' @export
create.solve_sdp_M <- function(Sigma, M=1, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma, check.attributes = F))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]

  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }

  # Convert problem for SCS

  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)

  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),]
  Cs = c((M+1)/M*G) ##change from 2 to (M+1)/M

  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p #not sure if it should be changed - may be not as it is the dimention of the linear part.

  # Objective
  b = rep(1,p)

  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")

  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }

  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1

  # Compensate for numerical errors (feasibility)
  if(verbose) cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-8
  while ((psd==0) & (s_eps<=0.1)) {
    if (is_posdef((M+1)/M*G-diag(s*(1-s_eps),length(s)),tol=1e-9)) { ##change from 2 to (M+1)/M
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  s[s<0]=0
  if(verbose) cat("done. \n")

  # Verify that the solution is correct
  if (all(s==0)) {
    warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }

  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}

#' Approximate Solution to Semidefinite Programming for Knockoff Variable Creation
#' @importFrom Rdsdp dsdp
#' @importFrom Matrix Diagonal cov2cor
#' @importFrom gtools binsearch
#' @importFrom RSpectra eigs
#' @export
create.solve_asdp_M <- function(Sigma, M=1, max.size=500, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma, check.attributes = F))

  if(ncol(Sigma) <= max.size) return(create.solve_sdp_M(Sigma, M=M, gaptol=gaptol, maxit=maxit, verbose=verbose))

  # Approximate the covariance matrix as block diagonal
  if(verbose) cat(sprintf("Dividing the problem into subproblems of size <= %s ... ", max.size))
  cluster_sol = divide.sdp(Sigma, max.size=max.size)
  n.blocks = max(cluster_sol$clusters)
  if(verbose) cat("done. \n")

  # Solve the smaller SDPs corresponding to each block
  if(verbose) cat(sprintf("Solving %s smaller SDPs ... \n", n.blocks))
  s_asdp_list = list()
  if(verbose) pb <- utils::txtProgressBar(min = 0, max = n.blocks, style = 3)
  for(k in 1:n.blocks) {
    s_asdp_list[[k]] = create.solve_sdp_M(as.matrix(cluster_sol$subSigma[[k]]), M=M, gaptol=gaptol, maxit=maxit)
    if(verbose) utils::setTxtProgressBar(pb, k)
  }
  if(verbose) cat("\n")

  # Assemble the solutions into one vector of length p
  p = dim(Sigma)[1]
  idx_count = rep(1, n.blocks)
  s_asdp = rep(0,p)
  for( j in 1:p ){
    cluster_j = cluster_sol$clusters[j]
    s_asdp[j] = s_asdp_list[[cluster_j]][idx_count[cluster_j]]
    idx_count[cluster_j] = idx_count[cluster_j]+1
  }

  # Maximize the shrinkage factor
  if(verbose) cat(sprintf("Combinining the solutions of the %s smaller SDPs ... ", n.blocks))
  tol = 1e-9
  maxitr=1000
  gamma_range = c(seq(0,0.1,len=11)[-11],seq(0.1,1,len=10)) # change from 100 to 20 to make it accurate near 0 and scalable.
  #options(warn=-1)
  gamma_opt = gtools::binsearch( function(i) {
    G = (M+1)/M*Sigma - gamma_range[i]*diag(s_asdp)
    lambda_min = suppressWarnings(RSpectra::eigs(G, 1, which = "SR", opts = list(retvec = FALSE, maxitr=maxitr, tol=tol))$values)
    if (length(lambda_min)==0) {
      #lambda_min = 1  # Not converged
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(G)$values)
    }
    lambda_min
  }, range=c(1,length(gamma_range)) )
  s_asdp_scaled = gamma_range[min(gamma_opt$where)]*s_asdp
  options(warn=0)
  if(verbose) cat("done. \n")

  if(verbose) cat("Verifying that the solution is correct ... ")
  # Verify that the solution is correct
  if (!is_posdef((M+1)/M*Sigma-diag(s_asdp_scaled,length(s_asdp_scaled)))) {
    warning('In creation of approximate SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
    s_asdp_scaled = 0*s_asdp_scaled
  }
  if(verbose) cat("done. \n")

  # Return result
  s_asdp_scaled
}

#' Divide Covariance Matrix for Semidefinite Programming
#' @importFrom stats hclust cutree dist
#' @importFrom Matrix cov2cor
#' @export
divide.sdp <- function(Sigma, max.size) {
  # Convert the covariance matrix into a dissimilarity matrix
  # Add a small perturbation to stabilize the clustering in the case of highly symmetrical matrices
  p = ncol(Sigma)
  Eps = matrix(rnorm(p*p),p)*1e-6
  dissimilarity = 1 - abs(cov2cor(Sigma)+Eps)
  distance = as.dist(dissimilarity)

  # Hierarchical clustering
  fit = hclust(distance, method="single")
  # Cut tree into clusters of size smaller than a threshold
  n.blocks.min = 1
  n.blocks.max = ncol(Sigma)
  for(it in 1:100) {
    n.blocks = ceiling((n.blocks.min+n.blocks.max)/2)
    clusters = cutree(fit, k=n.blocks)
    size = max(table(clusters))
    if(size <= max.size) {
      n.blocks.max = n.blocks
    }
    if(size >= max.size) {
      n.blocks.min = n.blocks
    }
    if(n.blocks.min == n.blocks.max) {
      break
    }
  }

  # Merge small clusters
  clusters.new = merge.clusters(clusters, max.size)
  while(sum(clusters.new != clusters)>0) {
    clusters = clusters.new
    clusters.new = merge.clusters(clusters, max.size)
  }
  clusters = clusters.new

  # Create covariance submatrices for each cluster
  subSigma = vector("list", max(clusters))
  for( k in 1:length(subSigma) ) {
    indices_k = clusters==k
    subSigma[[k]] = Sigma[indices_k,indices_k]
  }

  # Return the cluster assignments and the cluster covariance submatrices
  structure(list(clusters=clusters, subSigma=subSigma), class='knockoff.clusteredCovariance')
}

#' Merge Clusters Based on Size
#' @export
merge.clusters <- function(clusters, max.size) {
  cluster.sizes = table(clusters)
  clusters.new = rep(0, length(clusters))
  g = 1
  g.size = 0
  for(k in 1:max(clusters)) {
    if(g.size + cluster.sizes[k] > max.size) {
      g = g + 1
      g.size = 0
    }
    clusters.new[clusters==k] = g
    g.size = g.size + cluster.sizes[k]
  }
  return(clusters.new)
}

#' Check for Positive Definiteness of a Matrix
#' @importFrom RSpectra eigs
#' @export
is_posdef = function(A, tol=1e-9) {
  p = nrow(matrix(A))

  if (p<500) {
    lambda_min = min(eigen(A)$values)
  }
  else {
    oldw <- getOption("warn")
    #options(warn = -1)
    lambda_min = suppressWarnings(RSpectra::eigs(A, 1, which="SM", opts=list(retvec = FALSE, maxitr=100, tol))$values)
    options(warn = oldw)
    if( length(lambda_min)==0 ) {
      # RSpectra::eigs did not converge. Using eigen instead."
      lambda_min = min(eigen(A)$values)
    }
  }
  return (lambda_min>tol*10)
}

