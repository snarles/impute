# Charles Zheng
# p-values for imputed data


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

# if p1, p2 ~ F, tries to transforms p2 to make it subuniform with prob 1-alpha
finvert <- function(p1, p2, alpha=1) {
  p1 <- sort(p1)
  mm <- length(p1)
  # find threshold values
  tvals <- c(qbeta(1-alpha/mm, 1:mm, mm:1),1)
  breaks <- c(0, p1, 1)
  p3 <- p2
  for (i in 1:(mm+1)) {
    filt <- (p2 <= breaks[i+1]) & (p2 >= breaks[i])
    p3[filt] <- tvals[i]
  }
  p3
}

# obtains p-values giving tau thresholding
ftau <- function(tau, z) {
  pnorm(z-tau) + pnorm(z+tau)-1
}

# converts counts to p values
counts2p <- function(tau, n, y1,y2) {
  l <- log(y1) - log(n-y1) + log(n-y2) - log(y2)
  l[is.na(l)] <- 0
  #v <- y1 * (n - y1)/n *(1/y1 + 1/(n-y1))^2 + y2*(n-y2)/n*(1/y2+1/(n-y2))^2
  v <- 1/y1 + 1/(n-y1) + 1/y2 + 1/(n-y2)
  z <- (abs(l) - tau)/sqrt(v)
  z[z <0] <- 0
  p <- 2-2*pnorm(z)
  p
}

lorhat <- function(n, y1,y2) {
  lor <- log(y1)-log(y2)+log(n-y2)-log(n-y1)
  lor[is.na(lor)] <- 0
  lor
}

v2inds <- function(v) {
  ans <- numeric()
  while (sum(v) > 0) {
    ans <- c(ans, (1:length(v))[v > 0])
    v[v > 0] <- v[v > 0]-1
  }
  ans
}



# test of above functions
n.its <- 0
n <- 1000
p1 <- .5
p2 <- .53
lor <- log(p1) - log(p2) - log(1-p1)+log(1-p2)
res <- numeric(n.its)
tau <- 0.1
for (i in 1:n.its) {
  res[i] <- counts2p(tau,n,rbinom(1,n,p1),rbinom(1,n,p2))
}
#plot(sort(res))
#title(paste(tau,lor))
#lines(c(1,n.its),0:1)


# data for distribution of haplotypes under logistic model
hap.distr <- function(haplos,n.blocks,len.block,v,v0) {
  haplos.ca <- haplos
  haplos.co <- haplos
  nn <- dim(haplos)[1]
  m <- n.blocks * len.block
  stnorm <- qnorm((1:1000 -.5)/1000)
  weights <- matrix(0,dim(haplos)[1],n.blocks)
  marginals.ca <- matrix(0,dim(haplos)[1],n.blocks)
  marginals.co <- matrix(0,dim(haplos)[1],n.blocks)
  for (j in 1:n.blocks) {
    weights[,j] <- haplos[,(j-1)*len.block+(1:len.block)] %*%
                     v[(j-1)*len.block+(1:len.block)]
  }
  mus <- apply(weights,2,mean)
  vars <- apply(weights,2,var)
  cmus <- sum(mus) - mus + v0
  cvars <- sum(vars) - vars
  for (j in 1:n.blocks) {
    temp <- (stnorm + cmus[j]) * sqrt(cvars[j])
    temp2 <- exp(t(temp) %x% t(t(rep(1,nn))))
    temp3 <- temp2 + exp(weights[,j])
    temp4 <- temp3/(1+temp3)
    temp5 <- apply(temp4,1,mean)
    marginals.ca[,j] <- temp5/sum(temp5)
    haplos.ca[,(j-1)*len.block+(1:len.block)] <-
      haplos.ca[,(j-1)*len.block+(1:len.block)]*marginals.ca[,j]
  }
  for (j in 1:n.blocks) {
    temp <- (stnorm + cmus[j]) * sqrt(cvars[j])
    temp2 <- exp(t(temp) %x% t(t(rep(1,nn))))
    temp3 <- temp2 + exp(weights[,j])
    temp4 <- 1/(1+temp3)
    temp5 <- apply(temp4,1,mean)
    marginals.co[,j] <- temp5/sum(temp5)
    haplos.co[,(j-1)*len.block+(1:len.block)] <-
      haplos.co[,(j-1)*len.block+(1:len.block)]*marginals.co[,j]
  }
  p.ca <- apply(haplos.ca,2,sum)
  p.co <- apply(haplos.co,2,sum)
  
  lor <- log(p.ca)-log(1-p.ca)+log(1-p.co)-log(p.co)
  lor[is.na(lor)] <- 0
  list(w=weights,m.ca=marginals.ca,m.co=marginals.co,lor=lor)
}

# sample cases and controls using logistic model
# v is vector controlling weighting
sam.haps <- function(haplos, n.blocks, len.block, v, v0, l, diagn=F) {
  nn <- dim(haplos)[1]
  m <- n.blocks * len.block
  cases <- matrix(0,l, n.blocks*len.block)
  controls <- matrix(0,l,n.blocks*len.block)
  i.ca <- 0
  i.co <- 0
  while (min(c(i.ca,i.co)) < l) {
    h <- haplos[cbind(floor(nn*runif(n.blocks)+1) %x% rep(1,len.block),1:m)]
    temp <- exp(h %*% v + v0)
    p <- temp/(1+temp)
    if (diagn) {print(p)}
    if (runif(1) < p) {
      i.ca <- i.ca + 1
      if (i.ca < l) { cases[i.ca,] <- h }
    }
    else {
      i.co <- i.co+1
      if (i.co < l) { controls[i.co,] <- h }
    }
  }
  list(ca=cases,co=controls)
}

# helper function for finding tag SNPS and imputing
proj.block<- function(hblock,hblock2,locs,val=F,nos=F) {
  record <- numeric(dim(hblock2)[1])
  cmean <- apply(hblock,2,mean)
  projected <- t(hblock2)
  projected <- 0*projected+cmean
  projected[locs,] <- t(hblock2[,locs])
  filt <- rep(1, dim(hblock)[1])
  subb <- t(t(hblock[,locs]))*2-1
  subb2 <- t(t(hblock2[,locs]))*2-1
  n.class <- 0
  while (sum(filt) > 0) {
    n.class <- n.class+1
    i <- which(filt==1)[1]
    temp <- t(t(subb)-subb[i,])
    p.temp <-  t(t(subb2)-subb[i,])
    filtt <- which(apply(abs(temp),1,sum)==0)
    p.filtt <- which(apply(abs(p.temp),1,sum)==0)
    filt[filtt] <- 0
    temp2 <- hblock[filtt,]
    if (length(filtt)==1) {temp2 <- t(temp2)}
    projected[,p.filtt] <- apply(temp2,2,mean)
    record[p.filtt] <- n.class
  }
  if (val) {
    return(sum((t(projected)-hblock2)^2))
  }
  if (nos) {
    return(record)
  }
  t(projected)
}

# returns a list of tag SNPs based on greedy method
find.tags <- function(haplos,n.blocks,len.block,n.tags=2) {
  taglist <- numeric(n.tags * n.blocks)
  for (j in 1:n.blocks) {
    hblock <- haplos[,(j-1)*len.block+(1:len.block)]
    locs <- c()
    for (k in 1:n.tags) {
      scores <- numeric(len.block)
      for (i in 1:len.block) {
        scores[i] <- proj.block(hblock,hblock,c(locs,i),T)
      }
      locs <- c(locs,order(scores)[1])
    }
    if (length(unique(locs)) < n.tags) {
      kk <- n.tags - length(unique(locs))
      rem <- (1:len.block)[-locs]
      locs <- c(unique(locs),rem[1:kk])
    }
    taglist[(j-1)*n.tags + (1:n.tags)] <- locs + (j-1)*len.block
  }
  taglist
}

# imputes the rest of the dat given the refernce (haplos)
impute.haps <- function(haplos,n.blocks,len.block, dat,taglist) {
  dat2 <- dat
  for (j in 1:n.blocks) {
    hblock <- haplos[,(j-1)*len.block+(1:len.block)]
    locs <- taglist[taglist <= j*len.block & taglist > (j-1)*len.block]-(j-1)*len.block
    hblock2 <- dat[,(j-1)*len.block+(1:len.block)]
    dat2[,(j-1)*len.block+(1:len.block)] <- proj.block(hblock,hblock2,locs)
  }
  dat2
}

#case and controls to p values
cc2ps <- function(tau,cases,controls) {
  l <- dim(cases)[1]
  y1 <- apply(cases,2,sum)
  y0 <- apply(controls,2,sum)
  counts2p(tau,l,y1,y0)
}

# sets v0 so that mean probability is bounded away from 0
setv0 <- function(haplos,n.blocks,len.block,v,eps,n.its=5) {
  nn <- dim(haplos)[1]
  m <- n.blocks * len.block
  v0 <- 0
  for (ii in 1:n.its) {
    ps <- numeric()  
    for (i in 1:20) {
      h <- haplos[cbind(floor(nn*runif(n.blocks)+1) %x% rep(1,len.block),1:m)]
      temp <- exp(h %*% v + v0)
      p <- temp/(1+temp)
      ps <- cbind(ps,p)
    }
    mu <- mean(ps)
    if (mu < eps || mu > (1-eps)) {
      targ <- eps
      if (mu > .5) {targ <- 1-eps}
      v0 <- v0 + log(targ/mu)
    }
  }
  v0
}

gen.haplos <- function(n.blocks,len.block,n.generations,pop.size,mut.rate) {
  haplos <- matrix(0,pop.size,m)
  m <- n.blocks * len.block
  for (ii in 1:n.generations) {
    mut.filt <- matrix(rbinom(pop.size*m,1,mut.rate),pop.size,m)
    haplos[mut.filt==1] <- 1-haplos[mut.filt==1]
    for (j in 1:n.blocks) {
      filt <- sample(pop.size,pop.size,T)
      haplos[,(j-1)*len.block+(1:len.block)] <- haplos[filt,(j-1)*len.block+(1:len.block)]
    }
  }
  haplos
}

sparse.v <- function(m,n.cause) {
  filt <- sample(m, n.cause)
  v <- numeric(m)
  v[filt] <- rnorm(n.cause)
  v
}


# MODEL

# number of haplotype blocks
n.blocks <- 200
len.block <- 5
m <- n.blocks * len.block

# population haplotypes
pop.size <- 20
haplos <- gen.haplos(n.blocks,len.block,30,pop.size,0.1)
v <- sparse.v(m,200)
v0 <- setv0(haplos,n.blocks,len.block,v,.2)
# approximate number of alternatives
info <- hap.distr(haplos,n.blocks,len.block,v,v0)
sum(abs(info$lor) > tau)

alts <- rep(F, m)
n.calib <- 1000
calib <- sam.haps(haplos,n.blocks,len.block,v,v0,n.calib)
# note: originally, calib was a 'calibration test sample'
# however, now it actually generates the case and control populations
lor.calib <- lorhat(n.calib, apply(calib$ca,2,sum),apply(calib$co,2,sum))

#plot(lor.calib,info$lor)
alts[abs(lor.calib) > tau] <- T
m1 <- sum(alts)

# SIMULATED EXPERIMENTS
tau <- 0.05
# get a smaller reference panel
ref.size <- floor(pop.size * .1)
refp <- haplos[sample(pop.size,ref.size),]
# tag SNPS
taglist <- find.tags(refp,n.blocks,len.block,2)
# large-sample properties
caca <- impute.haps(haplos,n.blocks,len.block,calib$ca,taglist)
caco <- impute.haps(haplos,n.blocks,len.block,calib$co,taglist)
lor.i <- lorhat(n.calib, apply(caca,2,sum),apply(caco,2,sum))
tfilt <- rep(F,m); tfilt[taglist] <- T
impute.error <- sum((abs(lor.i) < tau) & (abs(lor.calib) > tau) & !tfilt)
impute.power <- sum((abs(lor.i) > tau) & (abs(lor.calib) > tau) & !tfilt)

# get cases and controls
l <- 20000
l.val <- 100

ind.cases <- sample(n.calib,l,T)
ind.controls <- sample(n.calib,l,T)
cases <- calib$ca[ind.cases,]
controls <- calib$co[ind.controls,]
cases.val <- calib$ca[sample(n.calib,l.val,T),]
controls.val <- calib$co[sample(n.calib,l.val,T),]
cases.im <- caca[ind.cases,]
controls.im <- caco[ind.cases,]
ps <- cc2ps(tau, cases,controls)
tests <- bhq(ps,0.1)
sum(alts[tests$hs.rejected])
sum(1-alts[tests$hs.rejected])
ps.i <- cc2ps(tau, cases.im, controls.im)
tests.i <- bhq(ps.i,0.1)
sum(alts[setdiff(tests.i$hs.rejected,taglist)])
sum(1-alts[tests.i$hs.rejected])

# implement our method
n.v <- 500
inds.v <- sample((1:m)[-taglist],n.v)
ps.screen <- cc2ps(tau,calib$ca[sample(n.calib,l.val,T),],calib$co[sample(n.calib,l.val,T),])
screening <- bhq(ps.screen[inds.v],.9)
null.samp <- inds.v[screening$hs.accepted]
ps.rest <- cc2ps(tau,cases.im,controls.im)
p.est <- ps.rest[null.samp]
#hist(p.est)
#plot(sort(p.est))
phats <- finvert(p.est,ps.rest,.5)
#plot(sort(phats))
#lines(c(1,m),c(0,1))
ourtest <- bhq(phats,0.5)
length(ourtest$hs.rejected)
sum(alts[setdiff(ourtest$hs.rejected,taglist)])

