###Adjusted mvrnorm so it handles p>n too

mvrnorm2 <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE) 
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  if (missing(EISPACK)) 
    EISPACK <- getOption("mvnorm_use_EISPACK", FALSE)
  eS <- eigen(Sigma, symmetric = TRUE, EISPACK = EISPACK)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    #addition
    if(p>n){
      X <- scale(X, TRUE, FALSE)
      X <- X %*% svd(X, nu = 0, nv = ncol(X))$v #set nv = ncol(X)
      X <- scale(X, FALSE, TRUE)
      print(90)
    }
    if (!p>n){
      X <- scale(X, TRUE, FALSE)
      X <- X %*% svd(X, nu = 0)$v
      X <- scale(X, FALSE, TRUE)
      print(200)
    }
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
}