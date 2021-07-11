#' Make the parameters in order to find a specific solution for a fixed lambda.
#'
#' In order to find a specific solution, we need some parameters A,B and C which are explained in the paper. This function makes the parameters for a fixed lambda.
#'
#' @param object ESPgenlasso object.
#' @param lambda a specific lambda
#' @param tol tolerance.
#' @return parameters A, B and C for the characterization.
#' @export
#'
make.abc <- function(object,lambda,tol = 1e-10)
{
  object.beta <- get.beta(object, lambda, tol)
  object.u <- get.u(object,lambda)
  active.set <- which(abs(object$D%*%object.beta) > tol)

  sgn.vec <- sign(object.u)
  if(length(active.set)==0){
    Fs <- ((object$D%*%object$W)*sgn.vec)
  } else{
    Fs <- ((object$D%*%object$W)*sgn.vec)[-active.set,]
  }
  E <- matrix(1,nrow=1,ncol=length(active.set))%*%((object$D%*%object$W)*as.numeric(sign(object$D%*%object.beta)))[active.set,]
  Gs <- matrix(1,nrow=1,ncol=dim(Fs)[1])%*%Fs
  E.Gs <- E+Gs
  dim.n <- dim(E.Gs)[1]
  dim.p <- dim(E.Gs)[2]
  E.Gs1 <- rbind(E.Gs,matrix(rep(0,(dim.p-dim.n)*dim.p),ncol=dim.p))
  Hs <- matrix(.null_ftn(t(E.Gs1),tol),nrow=dim(E.Gs1)[2])
  M <- as.numeric(sign(object$D%*%object.beta)[active.set,])*((object$D%*%object$W)[active.set,])
  N <- (-1*as.numeric(sign((object$D%*%object.beta))[active.set,]))*((object$D%*%object.beta)[active.set,])

  A <- rbind(Fs%*%Hs,M%*%Hs)
  B <- matrix(c(rep(0,dim(A)[1]-length(N)),N),ncol=1)
  C <- object$W%*%Hs

  return(list(A = A, B = B, C = C))
}
