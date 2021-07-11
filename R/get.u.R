#' Get the dual solution of the genlasso problem at a fixed lambda.
#'
#' This function gets the dual solution of the genlasso problem at a fixed lambda.
#'
#' @param object ESPgenlasso object.
#' @param lambda a specific lambda.
#' @return the dual solution at lambda.
#' @examples
#' y <- matrix(c(1000), nrow = 1)
#' X <- matrix(c(1,1,0),nrow = 1)
#' D <- matrix(c(1,0,0,0,1,1,0,0,-1), nrow = 3)
#' object <- ESPgenlasso(y,X,D,genlasso.option = FALSE)
#' get.u(object, 700)
#' @export
#'
get.u <- function(object, lambda)
{
  object.index <- which(object$lambda == lambda)
  if(length(object.index) == 0){
    object.index <- which(object$lambda < lambda)[1] - 1
    if(object.index == 0){
      object.u <- object$u[, 1]
    } else{
      object.u <- -(object$u[, object.index]-object$u[, (object.index+1)])/
        (object$lambda[object.index] - object$lambda[object.index+1]) * (object$lambda[object.index]-lambda) + object$u[, object.index]
    }
  } else{
    object.u <- object$u[, object.index[1]]
  }
  return(object.u)
}
