#' Get the solution that is a (partial) minimal (or/and) maximal l1 norm solution for some lambda, when the lasso penalty is included in the penalty matrix.
#'
#' This function gets the solution that is a (partial) minimal (or/and) maximal l1 norm solution for some lambda, when the lasso penalty is included in the penalty matrix.
#'
#' @param object ESPgenlasso object.
#' @param min.indices beta coef indices that would expected to have the minimum l1 norm.
#' @param max.indices beta coef indices that would expected to have the maximum l1 norm.
#' @param lambda.seq lambda seq.
#' @param tol tolerance.
#' @return solutions for each lambda seq.
#' @export
#'
spec.l1.solution <- function(object, lambda.seq, min.indices = c(), max.indices = c(), tol = 1e-10)
{
  contain.lasso <- .check.lasso(object$D)
  spec.beta.list <- list()
  if(contain.lasso$check && (length(min.indices)!=0 || length(max.indices)!=0) && (length(object$W)!=0)){
    spec.beta.indices <- 0
    for(lambda in lambda.seq){
      spec.beta.indices <- spec.beta.indices + 1
      object.beta <- get.beta(object, lambda , tol)
      object.u <- get.u(object,lambda)
      sgn.vec <- sign(object.u)

      char.param <- make.abc(object, lambda, tol)
      A <- char.param$A
      B <- char.param$B
      C <- char.param$C

      beta.sign.vec <- sgn.vec[contain.lasso$indices]
      d <- as.numeric(beta.sign.vec)[min.indices]%*%object.beta[min.indices]-as.numeric(beta.sign.vec)[max.indices]%*%object.beta[max.indices]
      a <- (as.numeric(beta.sign.vec)[min.indices]%*%C[min.indices,])-(as.numeric(beta.sign.vec)[max.indices]%*%C[max.indices,])

      cons <- optiSolve::lincon(A, d=rep(0, nrow(A)), dir=c(rep(">=",nrow(A))), val=B,
                                use=rep(TRUE,nrow(A)),name=c(1:nrow(A)))
      loss <- optiSolve::linfun(a, d=as.numeric(d), id=1:length(a), name="lin.fun")
      op <- optiSolve::cop(loss,lc=cons)
      op.sol <- optiSolve::solvecop(op,solver = "alabama",maxit = 500,itmax = 500,ilack.max = 500, quiet=TRUE)

      add_coef <- C%*%op.sol$x
      spec.beta <- object.beta + add_coef
      spec.beta.list[[spec.beta.indices]] <- spec.beta
    }
  } else{
    if(contain.lasso$check){
      if(length(object$W)==0){
        warning("Please set genlasso.option to False")
      } else{
        warning("Please set min.indices or max.indices.")
      }
    } else{
      warning("Not implemented in this penalty matrix.\nImplemented only for the penalty matrix including lasso.")
    }
  }
  return(spec.beta.list)
}

