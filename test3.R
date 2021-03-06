# Charles Zheng
# FDR for imputed data
# fixed hotspot model

rdirichlet <- function(v) {
  ans <- v*0
  for (i in 1:dim(v)[1]) {
    ans[i,] <- rgamma(dim(v)[2], v[i,])
  }
  ans <- apply(ans, 2, function(v) {v / sum(v)})
  ans
}

rmultinom2 <- function(v,n) {
  ans <- v*0
  for (i in 1:dim(v)[2]) {
    ans[,i] <- rmultinom(1, size=n, prob=v[,i])
  }
  ans
}

# ttest for binomial data
bttest <- function(y1,y2,n) {
  vy1 <- y1 * (n-y1)/n
  vy2 <- y2 * (n-y2)/n
  st <- (y1-y2)/(sqrt(vy1+vy2))
  2 * dt(abs(st), 2*n-2)
}

bhq <- function(p.values, q.thres) {
  m.obs <- sum(!is.na(p.values))
  ind.obs <- which(!is.na(p.values))
  p.complete <- p.values[ind.obs]
  o.pc <- order(p.complete)
  ps.sorted <- sort(p.complete)
  q.values <- ps.sorted * m.obs / (1:m.obs)
  no.reject <- max(c(0,which(q.values < q.thres)))
  hs.rejected <- ind.obs[o.pc[1:no.reject]]
  if (no.reject==0) { hs.rejected=numeric() }
  hs.accepted <- ind.obs[o.pc[(no.reject+1):m.obs]]
  list(hs.rejected=hs.rejected, hs.accepted=hs.accepted, q=q.values)
}

indep <- function(v) {
  p <- v[1]+v[2]
  q <- v[1]+v[3]
  c(p*q, p*(1-q), (1-p)*q, (1-p)*(1-q))
}

condit <- function(v) {
  c(v[c(1,2)]/(v[1]+v[2]),v[c(3,4)]/(v[3]+v[4]))
}

# if p1, p2 ~ F, tries to transforms p2 to make it subuniform with prob 1-alpha
finvert <- function(p1, p2, alpha=1) {
  p1 <- sort(p1)
  mm <- length(p1)
  # find threshold values
  tvals <- c(qbeta(1-alpha/mm, 1:mm, mm:1),1)
  breaks <- c(0, p1, 1)
  p3 <- p2
  for (i in 1:(mm+1)) {
    filt <- (p2 < breaks[i+1]) & (p2 >= breaks[i])
    p3[filt] <- tvals[i]
  }
  p3
}

# MODEL GENERATION PARAMETERS

b.lengths <- rep(10,5) # length of blocks
m <- sum(b.lengths) # total length of chromosome
b.breaks1 <- cumsum(c(0,b.lengths[-length(b.lengths)]))+1
b.breaks2 <- cumsum(b.lengths)
k <- 5 # number of clusters
cl.p <- as.vector(rdirichlet(t(t(rep(1,k))))) # cluster probabilities
clusters <- matrix(rbeta(k*m,.3,.3),k,m)




0













# number of pairs with non-independence
m1 <- 800

# accuracy of imputation rule
rule.accuracy <- 100

# closeness of conditional probs for case and controls
# case.control.closeness <- 100

# model parameters
control.probs <- rdirichlet(matrix(1,4,m))
control.probs[,(m1+1):m] <- apply(control.probs[,(m1+1):m], 2, indep)
case.probs <- control.probs
mix.probs <- rdirichlet(matrix(1,2,m)) %x% matrix(1,2,1)
case.probs <- control.probs * mix.probs
if (rule.accuracy > 0) {
  impute.probs <- rdirichlet(control.probs*rule.accuracy)
}
if (rule.accuracy==0) {
  impute.probs <-  rdirichlet(matrix(1,4,m))
}
# conditional impute probs by x
cimpute.probs <- impute.probs
cimpute.probs[c(1,2),] <- apply(impute.probs[c(1,2),],2,function(v) {v/sum(v)})
cimpute.probs[c(3,4),] <- apply(impute.probs[c(3,4),],2,function(v) {v/sum(v)})



# DATA

# same sample size for each group
n <- 200

# number of observed SNPS
m.obs <- 4000

# indices of observed SNPS
ind.obs <- sort(sample(m, m.obs, F))

case.complete <- rmultinom2(case.probs,n)
control.complete <- rmultinom2(control.probs,n)

# observed matrix of pair counts
case.observed <- case.complete
control.observed <- control.complete
case.observed[,-ind.obs] <- NA
control.observed[,-ind.obs] <- NA

# observed marginal counts
case.x <- apply(case.complete[1:2,],2,sum)
control.x <- apply(control.complete[1:2,],2,sum)
case.y <- apply(case.observed[c(1,3),],2,sum)
control.y <- apply(control.observed[c(1,3),],2,sum)



# PROCEDURE

# parameters for procedure
# FDR threshold for complete data
q.thres <- 0.1

# impute missing sites
case.imputed <- matrix(0,4,m)
control.imputed <- matrix(0,4,m)
case.imputed[c(1,2),] <- t(case.x*t(cimpute.probs[c(1,2),]))
case.imputed[c(3,4),] <- t((n-case.x)*t(cimpute.probs[c(3,4),]))
control.imputed[c(1,2),] <- t(control.x*t(cimpute.probs[c(1,2),]))
control.imputed[c(3,4),] <- t((n-control.x)*t(cimpute.probs[c(3,4),]))

case.y.imputed <- apply(case.imputed[c(1,3),],2,sum)
control.y.imputed <- apply(control.imputed[c(1,3),],2,sum)

# get p-values for complete data
p.values <- numeric(m)+NA
for (i in ind.obs) {
  temp <- prop.test(matrix(c(case.y[i], control.y[i], n-case.y[i],n-control.y[i]),2,2),alternative="two.sided")$p.value
  if (is.na(temp)) { temp <- 1}
  p.values[i] <- temp
}


# get p-values for imputed data
p.values.imputed <- numeric(m)
for (i in 1:m) {
  p.values.imputed[i] <-prop.test(matrix(c(case.y.imputed[i], control.y.imputed[i], n-case.y.imputed[i],n-control.y.imputed[i]),2,2),alternative="two.sided")$p.value
}
p.values.imputed[is.na(p.values.imputed)] <- 1

# carry out BH(q) procedure on complete data
res <- bhq(p.values,q.thres)
hs.rejected <- res$hs.rejected
sum(hs.rejected > m1)/length(hs.rejected)

# for fun, carry out BH(q) procedure on imputed data
res2 <- bhq(p.values.imputed, q.thres)
hs.r.imputed <- res2$hs.rejected
sum(hs.r.imputed > m1)/length(hs.r.imputed)

# find empirical distribution of imputed values for nulls
res1b <- bhq(p.values,3*q.thres)
hs.accepted <- res1b$hs.accepted
ps.imputed <- sort(p.values.imputed[hs.accepted])
#plot(grid,ps.imputed,type="l"); lines(0:1,0:1)

grid <- 1:length(ps.imputed)/(length(ps.imputed)+1)

ps.i2 <- grid
ps.i2[ps.imputed < grid] <- ps.imputed[ps.imputed < grid]

# fix the imputed distribution of p-values and carry out BHq on the imputed

p.i.fixed <- finvert(ps.i2, p.values.imputed, 0.1)

res3 <- bhq(p.i.fixed[-ind.obs], 0.6)
hs.r.i2 <- res3$hs.rejected
sum(hs.r.i2 > m1)/length(hs.r.i2)
length(hs.r.i2)


#(hs.r.imputed < m1) + 0
