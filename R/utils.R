.ginv_ftn <- function(X, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  if(any(nz)){
    rslt <- s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
  }
  else{
    rslt <- X
  }
  return(rslt)
}

.null_ftn <- function(M,tol=sqrt(.Machine$double.eps))
{
  tmp <- svd(M)
  ind <- which(abs(tmp$d) <tol)
  rslt <- tmp$u[,ind]
  return(rslt)
}

.check.lasso <- function(D)
{
  row.sum <- rowSums(abs(D))
  row.max <- apply(D, 1, max)
  row.argmax <- apply(D, 1, which.max)
  D.indices <- which((row.sum - row.max) == 0)
  sort.result <- sort(row.argmax[D.indices], index.return=T)
  lasso.coef.indices <- sort.result$x
  D.lasso.indices <- D.indices[sort.result$ix]
  if(length(lasso.coef.indices) == dim(D)[2]){
    check.seq <- c(1:(dim(D)[2]))
    if(sum(lasso.coef.indices-check.seq) == 0){
      return(list(check = T, indices = D.lasso.indices))
    } else{
      return(list(check = F, indices = D.lasso.indices))
    }
  } else{
    return(list(check = F, indices = D.lasso.indices))
  }
}

.hit_ftn <- function(bf.u,bf.lbd,slope,bd.set,m)
{
  if(length(bd.set)==0){
    int.set <- c(1:m)
  }
  else{
    int.set <- c(1:m)[-bd.set]
  }
  over.part <- bf.u[int.set] - slope*bf.lbd

  under.part1 <- 1-slope
  under.part2 <- -1-slope

  rslt1 <- over.part/under.part1
  rslt2 <- over.part/under.part2

  hit.rslt <- rep(0,length(int.set))

  for(i in 1:length(rslt1))
  {
    if(rslt1[i]>0 && rslt2[i]>0)
    {
      if(rslt1[i] > rslt2[i]){
        rslt1[i] <- -1
      }
      else{
        rslt2[i] <- -1
      }
    }
  }

  hit.rslt[which(rslt1>0)] <- rslt1[which(rslt1>0)]
  hit.rslt[which(rslt2>0)] <- rslt2[which(rslt2>0)]
  hit.rslt[which(hit.rslt>bf.lbd)] <- bf.lbd

  return(list(lbd=max(hit.rslt),coord = int.set[which.max(hit.rslt)]))
}

.lv_ftn <- function(bf.u,bf.lbd,til.y,m,l.d,til.d,w.mat,slope,delta,delta.slope,bd.set,sgn.set)
{
  int.size <- m - length(bd.set)
  null.size <- dim(w.mat)[2]
  tol=sqrt(.Machine$double.eps)
  if(int.size>0){
    f.part <- matrix(til.d[bd.set,],ncol=dim(til.d)[2])%*%(
      matrix(til.y,ncol=1) - t(matrix(til.d[-bd.set,],ncol=dim(til.d)[2]))%*%matrix(bf.u[-bd.set],ncol=1) +
        bf.lbd*t(matrix(til.d[-bd.set,],ncol=dim(til.d)[2]))%*%matrix(slope,ncol=1)) - ((l.d[bd.set,]%*%w.mat))%*%(matrix(delta,ncol=1) -
                                                                                                                     bf.lbd*matrix(delta.slope,ncol=1))
    s.part <- matrix(til.d[bd.set,],ncol=dim(til.d)[2])%*%(
      t(matrix(til.d[bd.set,],ncol=dim(til.d)[2]))%*%matrix(sgn.set,ncol=1) +
        t(matrix(til.d[-bd.set,],ncol=dim(til.d)[2]))%*%matrix(slope,ncol=1)) + ((l.d[bd.set,]%*%w.mat))%*%matrix(delta.slope,ncol=1)

    f.part <- f.part * sgn.set
    s.part <- s.part * sgn.set
  }
  else{
    d.bd <- matrix(l.d[bd.set,],ncol=dim(l.d)[2])
    til.d.bd <- matrix(til.d[bd.set,],ncol=dim(til.d)[2])
    f.part <- til.d.bd%*%til.y - d.bd%*%w.mat%*%delta
    f.part <- f.part*sgn.set
    s.part <- til.d.bd%*%t(til.d.bd)%*%sgn.set
    s.part <- s.part*sgn.set
  }

  temp.lbd <- c(0)
  temp.coord <- c(1)
  f.part[which(abs(f.part)<tol)] <- 0
  s.part[which(abs(s.part)<tol)] <- 0

  for(i in 1:length(bd.set))
  {
    if(f.part[i]<0){
      if(s.part[i]<0){
        if(temp.lbd < f.part[i]/s.part[i]){
          temp.lbd <- f.part[i]/s.part[i]
          temp.coord <- i
        }
      }
    }
  }
  if(temp.lbd>bf.lbd){
    temp.lbd <- bf.lbd
  }
  return(list(lbd=temp.lbd,coord=bd.set[temp.coord]))
}

.mult_element_ftn <- function(A.mat, B.mat){
  for(i in c(1:dim(A.mat)[1]))
  {
    B.mat[i,] <- B.mat[i,]*A.mat[i,]
  }
  return(B.mat)
}


.gen_mat_ftn <- function(l.d,inv.x,w.mat,bd.set,tol)
{
  til.d <- l.d%*%inv.x
  null.size <- dim(w.mat)[2]
  if(length(bd.set)>0){
    temp.mat <- cbind((matrix(til.d[-bd.set,],ncol=dim(til.d)[2])%*%t(matrix(til.d[-bd.set,],ncol=dim(til.d)[2]))),(l.d[-bd.set,]%*%w.mat))
    temp.mat <- rbind(temp.mat,cbind(t((l.d[-bd.set,]%*%w.mat)),matrix(rep(0,(null.size*null.size)),nrow=null.size)))
  }
  else{
    temp.mat <- cbind((til.d%*%t(til.d)),(l.d%*%w.mat))
    temp.mat <- rbind(temp.mat,cbind(t((l.d%*%w.mat)),matrix(rep(0,(null.size*null.size)),nrow=null.size)))
  }
  if(dim(temp.mat)[1]>10000){
    inv.mat <- .svdsolve_ftn(temp.mat,5000,tol,T)$x
  }
  else{
    inv.mat <- tryCatch({
      inv.mat <- .svdsolve_ftn(temp.mat,5000,tol,F)$x
    }, error = function(err) {
      err$message = paste(err$message, "\n(svd solve using package)", sep = "")
      warning(err)
      inv.mat <- .svdsolve_ftn(temp.mat,min(Matrix::rankMatrix(temp.mat)[1],5000),tol,T)$x
      return(inv.mat)
    })
  }
  return(inv.mat)
}

.svdsolve_ftn <- function(A,nu,rtol = 1e-5,pack.use=F,propack=F) {
  if(pack.use){
    if(propack){
      s= svd::propack.svd(A, verbose=FALSE)
    }
    else{
      maxit = 5
      s = suppressWarnings(irlba::irlba(A, nu, tol=1e-5, maxit = maxit, verbose=FALSE))
      if(s$iter==maxit){
        s = suppressWarnings(irlba::irlba(A, nu, tol=1e-5, maxit = maxit, verbose=FALSE, fastpath=FALSE))
        if(s$iter==maxit){
          for(i in c(1:(nu))){
            if(s$iter==maxit){
              s= suppressWarnings(irlba::irlba(A, nu-i, tol=1e-5, maxit = maxit, verbose=FALSE))
            }
            else{
              break
            }
          }
        }
      }
    }
  }
  else{
    s = svd(A)
  }
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v%*%(di*(t(s$u))),q=sum(ii)))
}

.svdsolve <- function(A,tol){
  if(dim(A)[1]>10000){
    inv.mat <- .svdsolve_ftn(A,5000,tol,T,F)$x
  }
  else{
    inv.mat <- tryCatch({
      inv.mat <- .svdsolve_ftn(A,5000,tol,F,F)$x
    }, error = function(err) {
      inv.mat <- tryCatch({
        inv.mat <- .svdsolve_ftn(A, (Matrix::rankMatrix(A)[1]), tol, T, F)$x
      },error = function(err) {
        inv.mat <- .svdsolve_ftn(A, (Matrix::rankMatrix(A)[1]), tol, T, T)$x
        return(inv.mat)
      }
      )
      return(inv.mat)
    })
  }
  return(inv.mat)
}
