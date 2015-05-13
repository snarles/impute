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
## now with hypergem.
finvert <- function(p1, p2, alpha=.1) {
  p1 <- sort(p1)
  mm <- length(p1)
  nn <- length(p2)
  gg <- c((1:mm)/mm,1)
  ff <- mm/nn
  tvals <- rep(1,mm+1)
  for (i in 1:mm) {
    temp <- phyper(i,1:nn,nn-(1:nn),mm,T)
    tvals[i] <- max(which(temp >= (1-alpha)))/nn
  }
  # find threshold values
  #tvals <- c(qbeta(1-alpha, 1:mm, mm:1),1)
  #tvals <- c(((1:mm)^2 + ff^2)/(2*(1:mm)*ff+ff*(1-ff)*qnorm(1-alpha)^2),nn)/nn
  #tvals[tvals < gg] <- gg[tvals < gg]
  breaks <- c(-Inf, p1, Inf)
  p3 <- p2
  for (i in 1:(mm+1)) {
    filt <- (p2 <= breaks[i+1]) & (p2 >= breaks[i])
    p3[filt] <- tvals[i]
  }
  p3
}

# converts counts to p values or z scorss
counts2p <- function(tau, n, y1,y2, zs=F) {
  l <- log(y1) - log(n-y1) + log(n-y2) - log(y2)
  l[is.na(l)] <- 0
  #v <- y1 * (n - y1)/n *(1/y1 + 1/(n-y1))^2 + y2*(n-y2)/n*(1/y2+1/(n-y2))^2
  v <- 1/y1 + 1/(n-y1) + 1/y2 + 1/(n-y2)
  z <- (abs(l) - tau)/sqrt(v)
  z[z <0] <- 0
  if (zs) {return(z)}
  p <- 2-2*pnorm(z)
  p
}

eq <- function(a,b) {
  abs(a-b) < 1e-10
}

lorhat <- function(n, y1,y2) {
  if (length(n)==1) {n <- rep(n,length(y1))}
  filt <- eq(y1,n) | eq(y1,0) | eq(y2,n) | eq(y2,0)
  lor <- numeric(length(y1))
  lor[eq(y1,n) || eq(y2,0)] <- Inf
  lor[eq(y2,n) || eq(y1,0)] <- -Inf
  lor[!filt] <- log(y1[!filt])-log(y2[!filt])+log(n[!filt]-y2[!filt])-log(n[!filt]-y1[!filt])
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

# get a logistic sum to equal a target
findv0 <- function(vec,targ) {
  f <- function(x) { (targ - sum(exp(vec+x)/(1+exp(vec+x))))^2 }
  grid <- t(-100:100)^3
  scores <- apply(grid,2,f)
  inds <- order(scores)[1:3]
  lb <- min(grid[inds])
  ub <- max(grid[inds])
  res <- optimize(f,lower=lb,upper=ub)
  res$minimum
}

# form a population with phenotypes, with a certain case probability
phenopop <- function(haplos,n.blocks,len.block,v,nlarge,case.prob) {
  nn <- dim(haplos)[1]
  m <- n.blocks * len.block
  mat <- matrix(0,nlarge,m)
  for (j in 1:n.blocks) {
    mat[,(j-1)*len.block+(1:len.block)] <-
      haplos[sample(nn,nlarge,T),(j-1)*len.block+(1:len.block)]
  }
  temp <- mat %*% v
  v0 <- findv0(temp,case.prob*length(temp))
  ps <- exp(temp+v0)/(1+exp(temp+v0))
  p.case <- ps/sum(ps)
  snps.ca <- t(p.case) %*% mat
  p.control <- (1-ps)/sum(1-ps)
  snps.co <- t(p.control) %*% mat
  lor <- lorhat(1,snps.ca,snps.co)
  sam.ca <- function(l,inds=F) {
    ind.ca <- rmultinom(1,l,p.case)
    if (inds) {return(ind.ca)}
    mat[v2inds(ind.ca),]
  }
  sam.co <- function(l,inds=F) {
    ind.co <- rmultinom(1,l,p.control)
    if (inds) {return(ind.co)}
    mat[v2inds(ind.co),]
  }
  list(mat=mat,ps=ps,lor=lor,ca=sam.ca,co=sam.co,nb=n.blocks,
       lb=len.block,p.case=p.case,p.control=p.control)
}

# takes a phenopop object and a taglist, and returns imputation distributions
imputepop <- function(pp,taglist) {
  mat.i <- impute.haps(pp$mat,pp$nb,pp$lb,pp$mat,taglist)
  snps.ca.i <- t(pp$p.case) %*% mat.i
  snps.co.i <- t(pp$p.control) %*% mat.i
  lor.i <- lorhat(1,snps.ca.i,snps.co.i)
  list(mat.i=mat.i, lor.i=lor.i)
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

#cc2ps given matrix and counts
cc2ps2 <- function(tau,mat,ca.count,co.count,zs=F) {
  l <- sum(ca.count)
  y1 <- t(ca.count) %*% mat
  y0 <- t(co.count) %*% mat
  counts2p(tau,l,y1,y0,zs)
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

fdrplot <- function(ps,alts) {
  a <- alts[order(ps)]
  aa <- cumsum(a)
  plot(1:length(ps),1-aa/(1:length(ps)),type="l",xlab="No. rejections",ylab="fdp")
}

# MODEL

set.seed(1)

# number of haplotype blocks
n.blocks <- 100
len.block <- 10
m <- n.blocks * len.block

# population haplotypes
pop.size <- 50
haplos <- gen.haplos(n.blocks,len.block,30,pop.size,0.2)
v <- 10*sparse.v(m,2)
# generate case and control populations
pp <- phenopop(haplos,n.blocks,len.block,v,3000,0.2)
lor <- pp$lor
pdf("impute_g1.pdf")
plot(sort(abs(lor)),type="l")
points((1:m)[sort(abs(lor))>tau], sort(abs(lor))[sort(abs(lor))>tau])
lines(c(1,m),c(tau,tau),col="red")
title("True absolute log OR")
dev.off()

taglist <- find.tags(haplos,n.blocks,len.block,10)
pdf("impute_g2.pdf")
tagg <- numeric(10)
for (i in 1:10) {
  filt <- rep(10*(1:n.blocks -1),each=i)+rep(1:i,n.blocks)
  haplos.i <- impute.haps(haplos,n.blocks,len.block,haplos,taglist[filt])
  tagg[i] <- sum((haplos-haplos.i)^2)/length(haplos)
}
plot(tagg,xlab="tag SNPS per block",ylab="MSE",type="b")
title("Imputation Error vs. # of tag SNPs/block")
dev.off()

maf <- .5-abs(apply(haplos,2,sum)/pop.size-.5)

pdf("impute_g3.pdf")
hist(maf, main="Minor allele frequency in reference", xlab="Minor allele frequency")
dev.off()




# SIMULATED EXPERIMENTS
tau <- 0.1
alts <- (abs(lor) > tau)
m1 <- sum(alts)
# get a smaller reference panel
ref.size <- floor(pop.size * .5)
refp <- haplos[sample(pop.size,ref.size,T),]
refp <- haplos
# errors in reference panel
e.rate <- 0.0
mut.filt <- matrix(rbinom(ref.size,1,e.rate),ref.size,m)
refp[mut.filt==1] <- 1-refp[mut.filt==1]

# tag SNPS
taglist <- find.tags(refp,n.blocks,len.block,4)
alts2 <- alts[-taglist]
ppi <- imputepop(pp,taglist)
lor.i <- ppi$lor.i

pdf("impute_g7.pdf")
plot(lor,lor.i,ylab="Imputed",xlab="True",pch=".",xlim=3*tau*c(-1,1),ylim=3*tau*c(-1,1))
filt1 <- (abs(lor) > tau)
filt2 <- (abs(lor.i) > tau)
points(lor[!filt1 & filt2],lor.i[!filt1 & filt2],pch="x",col="violet")
points(lor[filt1 & !filt2],lor.i[filt1 & !filt2],pch="X",col="red")
points(lor[filt1 & filt2],lor.i[filt1 & filt2],pch="o",col="green")
title("True log-odds ratio vs. large-sample imputed LOR")
dev.off()

tfilt <- rep(F,m); tfilt[taglist] <- T
impute.error <- sum((abs(lor.i) < tau) & (abs(lor) > tau) & !tfilt)
impute.power <- sum((abs(lor.i) > tau) & (abs(lor) > tau) & !tfilt)

# get cases and controls
l <- 10000
l.val <- 1000
source("sub1.R")

ca.count <- pp$ca(l,T)
co.count <- pp$co(l,T)
ps <- cc2ps2(tau,pp$mat, ca.count,co.count)
#plot(ps,abs(lor))
tests <- bhq(ps,0.1)
sum(alts[tests$hs.rejected])
sum(1-alts[tests$hs.rejected])
ps.i <- cc2ps2(tau, ppi$mat.i, ca.count,co.count)
zs.i <- cc2ps2(tau, ppi$mat.i, ca.count,co.count,T)
tests.i <- bhq(ps.i[-taglist],0.2)
length(tests.i$hs.rejected)
sum(!alts2[tests.i$hs.rejected])

pdf("impute_g8.pdf")
fdrplot(ps,alts)
title("False discovery proportion: complete genotype data", sub="10000 cases, 10000 controls")
dev.off()

pdf("impute_g9.pdf")
fdrplot(ps.i[-taglist],alts2)
title("False discovery proportion: imputed genotype data", sub="10000 cases, 10000 controls")
dev.off()

# implement our method
n.v <- floor((m-length(taglist))*1)
inds.v <- sample((1:m)[-taglist],n.v)
ps.val <- cc2ps2(tau,pp$mat,pp$ca(l.val,T),pp$co(l.val,T))
#plot(abs(lor)[inds.v],ps.val[inds.v])
val.res <- bhq(ps.val[inds.v],.5)
#null.samp <- inds.v[order(-ps.val)[1:floor(n.v * .8)]]
null.samp <- inds.v[val.res$hs.accepted]
p.est <- zs.i[null.samp]
#hist(p.est)
#plot(sort(p.est))
phats <- finvert(-p.est,-zs.i,0.1)
#plot(sort(phats[!alts]))
plot(sort(phats[-taglist]))
lines(c(1,m-m1),c(0,1))
ourtest <- bhq(phats[-taglist],0.2)
length(ourtest$hs.rejected)
sum(!alts2[ourtest$hs.rejected])





