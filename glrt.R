lor <- function(p,q) {
  log(p)-log(1-p)+log(1-q)-log(q)
}
p2q <- function(p,theta) {
  p/(exp(theta)+(1-exp(theta))*p)
}

n.its <- 1000
p <- .1
q <- .1
theta0 <- lor(q,p)
n <- 1000
us <- rbinom(n.its,n,p)
ts <- rbinom(n.its,n,q)
vs <- log(us)-log(n-us)
ws <- log(ts)-log(n-ts)
thetas <- ws-vs

