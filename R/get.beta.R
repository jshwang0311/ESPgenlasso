#' Get the primal solution of the genlasso problem at a fixed lambda.
#'
#' This function gets the primal solution of the genlasso problem at a fixed lambda.
#'
#' @param object ESPgenlasso object.
#' @param lambda a specific lambda.
#' @param tol tolerance.
#' @return the primal solution at lamdba.
#' @examples
#' y <- matrix(c(1000), nrow = 1)
#' X <- matrix(c(1,1,0),nrow = 1)
#' D <- matrix(c(1,0,0,0,1,1,0,0,-1), nrow = 3)
#' object <- ESPgenlasso(y,X,D,genlasso.option = FALSE)
#' get.beta(object, 700)
#' @export
#'
get.beta <- function(object, lambda, tol = 1e-10)
{
  object.index <- which(object$lambda == lambda)
  if(length(object.index) == 0){
    object.index <- which(object$lambda < lambda)[1] - 1
    if(object.index == 0){
      object.beta <- object$beta[, 1]
      object.beta[which(abs(object.beta) <= tol)] <- 0
    } else{
      object.beta <- -(object$beta[, object.index]-object$beta[, (object.index+1)])/
        (object$lambda[object.index] - object$lambda[object.index+1]) * (object$lambda[object.index]-lambda) + object$beta[, object.index]
      object.beta[which(abs(object.beta) <= tol)] <- 0
    }
  } else{
    object.beta <- object$beta[, object.index[1]]
    object.beta[which(abs(object.beta) <= tol)] <- 0
  }
  return(object.beta)
}
