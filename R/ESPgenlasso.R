#' Get Exact Solution Path for genlasso problem.
#'
#' This function solves the genlasso problem for all lambda. So, we can get the solution path object which contains dual solution u, break point lambda, primal solution beta and so on.
#'
#' @param y a numeric response vector.
#' @param X a matrix of predictor variables.
#' @param D a penalty matrix.
#' @param genlasso.option whether to use the genlasso package.
#' @param thres.lambda minimum threshold lambda.
#' @param iter maximum threshold iteration number.
#' @param tol tolerance.
#' @return Exact Solution Path of genlasso object (ESPgenlasso object).
#' @examples
#' y <- matrix(c(1000), nrow = 1)
#' X <- matrix(c(1,1,0),nrow = 1)
#' D <- matrix(c(1,0,0,0,1,1,0,0,-1), nrow = 3)
#' ESPgenlasso(y,X,D,genlasso.option = FALSE)
#' @export
#'
ESPgenlasso <- function(y, X, D, genlasso.option=F, thres.lambda = 0.1, iter = 100000, tol=1e-12)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(D)[1]
  lbd <- c()

  lbd[1] <- .Machine$double.xmax
  inv.x <- .ginv_ftn(t(X) %*% (X),tol) %*% t(X)
  til.d <- D%*%inv.x
  til.y <- X%*%inv.x%*%y
  if(n<p){
    X1 <- rbind(X,matrix(rep(0,(p-n)*p),ncol=p))
    w.mat <- matrix(.null_ftn(t(X1),tol),nrow=dim(X1)[2])
  }
  if(n>=p){
    w.mat <- matrix(.null_ftn(t(X),tol),nrow=dim(X)[2])
  }
  null.size <- dim(w.mat)[2]
  if(null.size==0){
    pack.genlasso.rslt <- genlasso::genlasso(y, X, D)
    out <- list(u = pack.genlasso.rslt$u,delta = NULL,delta_slope = NULL, lambda=pack.genlasso.rslt$lambda, hit=pack.genlasso.rslt$hit,
                beta=pack.genlasso.rslt$beta, fit=pack.genlasso.rslt$fit,
                X = X, y = y, D=D, W=NULL)
    class(out) <- c("ESPgenlasso", "list")
    return(out)
  } else if(genlasso.option){
    pack.genlasso.rslt <- genlasso::genlasso(y, X, D,svd=T,eps=1e-8)
    out <- list(u = pack.genlasso.rslt$u,delta = NULL,delta_slope = NULL, lambda=pack.genlasso.rslt$lambda, hit=pack.genlasso.rslt$hit,
                beta=pack.genlasso.rslt$beta, fit=pack.genlasso.rslt$fit,
                X = X, y = y, D=D, W=NULL)
    class(out) <- c("ESPgenlasso", "list")
    return(out)
  } else{
    u.mat <- matrix(nrow=iter,ncol=m)
    delta.mat <- matrix(nrow = iter,ncol=null.size)
    delta.mat.slope <- matrix(nrow = iter,ncol=null.size)
    bd <- list()
    sgn <- list()
    solv.temp <- list()
    hit <- c()
    lbd <- c()
    bd.set <- c()
    sgn.set <- c()
    delta <- c()
    for(k in c(1:iter)){
      if(k==1){
        bd[[k]] <- bd.set
        sgn[[k]] <- sgn.set
        inv.mat <- .gen_mat_ftn(D,inv.x,w.mat,bd.set,tol)
        solv.mat <- (inv.mat%*%rbind((til.d%*%til.y),as.matrix(rep(0,null.size))))
        solv.temp[[k]] <- solv.mat

        u.mat[k,] <- solv.mat[1:m,]
        delta.mat[k,] <- solv.mat[c((dim(solv.mat)[1]-null.size + 1):(dim(solv.mat)[1])),]
        delta <- delta.mat[k,]
        lbd[k] <- max(abs(u.mat[k,]))
        bd.set <- which(abs(abs(u.mat[k,])-lbd[k])<tol)
        sgn.set <- sign(u.mat[k,bd.set])
        r.delta <- delta
        hit <- T
      } else{
        bd[[k]] <- bd.set
        sgn[[k]] <- sgn.set
        bf.u <- u.mat[(k-1),]
        bf.lbd <- lbd[k-1]
        int.size <- m-length(bd.set)
        int.set <- c(1:m)[-bd.set]

        if(length(bd.set)==m){
          delta <- r.delta
          delta.mat[k,] <- delta
          if(lbd[k-1]<=thres.lambda){
            break
          }
          slope <- 0
          delta.slope <- 0
          delta.mat.slope[k-1,] <- rep(0,null.size)
          lv.rslt <- .lv_ftn(bf.u,bf.lbd,til.y,m,D,til.d,w.mat,slope,delta,delta.slope,bd.set,sgn.set)
          if(lv.rslt$lbd==0){
            lbd[k] <- 0
            break
          } else{
            lbd[k] <- lv.rslt$lbd
            u.mat[k,-bd.set] <- bf.u[-bd.set] + (lbd[k] - bf.lbd)*slope
            u.mat[k,bd.set] <- lbd[k]*sgn.set

            bf.u <- u.mat[(k),]
            bf.lbd <- lbd[k]
            int.size <- m-length(bd.set)

            sgn.set <- sgn.set[-which(bd.set==lv.rslt$coord)]
            bd.set <- bd.set[-which(bd.set==lv.rslt$coord)]
            hit <- c(hit,F)
            r.delta <- delta + (lbd[k]-bf.lbd)*delta.slope
          }
        } else{
          delta <- r.delta
          delta.mat[k,] <- delta
          if(lbd[k-1]<=thres.lambda){
            break
          }
          temp.mat <- cbind((matrix(til.d[-bd.set,],ncol=dim(til.d)[2])%*%t(matrix(til.d[-bd.set, ], ncol=dim(til.d)[2]))),
                            (D[-bd.set,]%*%w.mat))
          temp.mat <- rbind(temp.mat, cbind(t((D[-bd.set,]%*%w.mat)), matrix(rep(0,(null.size*null.size)),nrow=null.size)))

          inv.mat <- .svdsolve(temp.mat,tol)

          mid.mat1 <- -1*matrix(til.d[-bd.set,],ncol=dim(til.d)[2])%*%(t(matrix(til.d[bd.set,],ncol=dim(til.d)[2]))%*%matrix(sgn.set,ncol=1))
          mid.mat2 <- -1*t(matrix(D[bd.set,],ncol=dim(D)[2])%*%w.mat)%*%matrix(sgn.set,ncol=1)
          mid.mat <- rbind(mid.mat1,mid.mat2)

          rslt.mat <- inv.mat%*%mid.mat


          slope <- rslt.mat[c(1:int.size),]
          delta.slope <- rslt.mat[c((int.size+1):dim(rslt.mat)[1]),]
          delta.mat.slope[k-1,] <- delta.slope

          if(hit[length(hit)]==F){
            lv.coord.slope <- slope[which(int.set==lv.rslt$coord)]
            lv.coord.sign <- sign(bf.u[lv.rslt$coord])
            if(lv.coord.slope*lv.coord.sign<1){
              break
            }
          }
          hit.rslt <- .hit_ftn(bf.u,bf.lbd,slope,bd.set,m)
          lv.rslt <- .lv_ftn(bf.u,bf.lbd,til.y,m,D,til.d,w.mat,slope,delta,delta.slope,bd.set,sgn.set)

          if(hit.rslt$lbd>lv.rslt$lbd){
            lbd[k] <- hit.rslt$lbd
            u.mat[k,-bd.set] <- bf.u[-bd.set] + (lbd[k] - bf.lbd)*slope
            u.mat[k,bd.set] <- lbd[k]*sgn.set

            sgn.set <- c(sgn.set,sign(u.mat[k, hit.rslt$coord[1]]))
            bd.set <- c(bd.set, hit.rslt$coord[1])

            u.mat[k, bd.set] <- lbd[k] * sgn.set
            r.delta <- delta + (lbd[k] - bf.lbd) * delta.slope
            hit <- c(hit,T)
          } else{
            if(lv.rslt$lbd==0){
              lbd[k] <- 0
              break
            }
            lbd[k] <- lv.rslt$lbd
            u.mat[k,-bd.set] <- bf.u[-bd.set] + (lbd[k] - bf.lbd)*slope
            u.mat[k,bd.set] <- lbd[k]*sgn.set

            bf.u <- u.mat[(k),]
            bf.lbd <- lbd[k]

            sgn.set <- sgn.set[-which(bd.set==lv.rslt$coord)]
            bd.set <- bd.set[-which(bd.set==lv.rslt$coord)]
            hit <- c(hit,F)
            r.delta <- delta + (lbd[k]-lbd[k-1])*delta.slope
          }
        }
      }
    }
    u.mat <- matrix(u.mat[!is.na(u.mat[,1]),],ncol=m)
    delta.mat <- matrix(delta.mat[!is.na(delta.mat[,1]),],ncol=null.size)
    delta.mat.slope <- matrix(delta.mat.slope[!is.na(delta.mat.slope[,1]),],ncol=null.size)

    delta.mat[which(abs(delta.mat)<tol)] <- 0
    delta.mat.slope[which(abs(delta.mat.slope)<tol)] <- 0
    rownames(u.mat) <- apply(u.mat,1,function(x) max(abs(x)))
    beta <- .ginv_ftn(t(X) %*% X, tol) %*% (matrix(rep(t(X) %*% y, times=length(lbd)), dim(X)[2], length(lbd)) - t(D) %*% t(u.mat)) - w.mat %*% (t(delta.mat)[,-1])
    beta <- cbind(beta , .ginv_ftn(t(X) %*% X, tol) %*% (t(X) %*% y) )
    colnames(beta)[dim(beta)[2]] <- "0"
    fit <- X %*% beta
    out <- list(u = t(u.mat), delta = (delta.mat), delta_slope = delta.mat.slope,
                lambda = lbd, hit = hit, beta = beta, fit = fit,
                X = X, y = y, D=D, W=w.mat)
    class(out) <- c("ESPgenlasso", "list")
    return(out)
  }
}
