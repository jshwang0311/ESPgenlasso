# ESPgenlasso (Exact Solution Path of generalized lasso problem for rank-deficient design matrices)
## Description
This R package provides a solution path of generalized lasso problem for rank-deficient design matrices, which is an exact entrie solution path. And also when the Lasso penalty is contained in the penalty matrix, the specific solutions, which are the (partial) minimal (or/and) maximal $l_1$ norm solution and the (partial) $l_{\infty}$ norm solution, are obtained. Details are described in the following paper:
- Jaesung Hwang, Joongho Won, Yongdai Kim, Characterization of the Solutions Set of the Generalized LASSO Problems for Non-Full Rank Cases.

## Install
In R terminal,

```R
install.packages("devtools")
devtools::install_github("jshwang0311/ESPgenlasso")
```

## Example
Import library, generate simulated data and make a penalty matrix using the following commands.
```R
library(ESPgenlasso)
library(mvtnorm)

n <- 100
p <- 500

X <- rmvnorm(n,sigma = diag(1,p))
beta.coef <- rep(0,p)
beta.coef[c((p/10+1):(2*p/10))] <- 1
beta.coef[c((2*p/10+1):(4*p/10))] <- 2
beta.coef <- matrix(beta.coef,ncol=1)
err.val <- matrix(rmvnorm(1,sigma = diag(1,n)),ncol=1)
y <- X%*%beta.coef + err.val

D <- matrix(0,ncol=p,nrow=200)
for(i in c(1:200)){
  D[i,(299+i)] <- -1
  D[i,(i+300)] <- 1
}
eps <- 0.1
temp.D <- diag(eps,p)
D <- rbind(D,temp.D)
```

Find the exact solution path with the following command.
```R
new.prop.rslt <- ESPgenlasso::ESPgenlasso(y,X,D)
```

To find a partial minimal and maximal $l_1$ norm solution for the given lambdas, enter the following command.
```R
# Set the indices of beta coef that would expected to have the maximum l1 norm
max.indices <- which(beta.coef==2)[c(1:100)]
# Set the indices of beta coef that would expected to have the minimum l1 norm
min.indices <- which(beta.coef==0)[c(1:50)]
lambda.seq <- c(100,30,20,10,5,3,1)
spec.l1 <- ESPgenlasso::spec.l1.solution(new.prop.rslt, lambda.seq, min.indices, max.indices)
```

To get a partial maximal l infinity norm solution for the given lambdas, use the following command.
```R
# Set beta coef indices that would expected to have the maximum l infinity norm.
max.indices <- which(beta.coef==2)[c(1:50)]
lambda.seq <- c(10,5,3,2,1)
spec.linf <- ESPgenlasso::spec.linf.solution(new.prop.rslt, lambda.seq, max.indices, parallel = T, numWorkers = 8)
```
