library(MASS)
library(corpcor)
library(tseries)
library(rJava)
library(xlsx)

#corr matrix from book 
cv=c(1,0.92,0.33,0.26,0.28,0.16,0.29,0.42,0.92,1,0.26,0.22,0.27,0.14,0.25,0.36,0.33,0.26,1,0.41,0.3,0.25,0.58,0.71,0.26,0.22,0.41,1,0.62,0.42,0.54,0.44,0.28,0.27,0.3,0.62,1,0.35,0.48,0.34,0.16,0.14,0.25,0.42,0.35,1,0.4,0.22,0.29,0.25,0.58,0.54,0.48,0.4,1,0.56,0.42,0.36,0.71,0.44,0.34,0.22,0.56,1)
cv=matrix(cv,ncol=8,nrow=8)
cv
mean=c(0.83,0.85,0.97,1.46,1.11,1.46,1.37,1.29)
SD=c(1.99,1.52,5.47,7,6.19,7.01,5.99,4.28)
#This is wrong because it uses corr matrix to generate data and not the covariance matrix
dgp=mvrnorm(216,mean,cv)
dgp=data.frame(dgp)
colnames(dgp)=c("Euro Bonds","US Bonds","Canada","France","Germany","Japan","UK","US")
dgp=as.matrix(dgp)
vcv=cov(dgp)
vcv
#Making cov mat out of corr mat
b=SD %*% t(SD)
acov=b*cv
acov

#perturbed truth
rdgp1=mvrnorm(216,mean,acov)
rdgp1=data.frame(rdgp1)
colnames(rdgp1)=c("Euro Bonds","US Bonds","Canada","France","Germany","Japan","UK","US")
rdgp1=as.matrix(rdgp1)



#drawing from truth
means=c(mean(rdgp1[,1]),mean(rdgp1[,2]),mean(rdgp1[,3]),mean(rdgp1[,4]),mean(rdgp1[,5]),mean(rdgp1[,6]),mean(rdgp1[,7]),mean(rdgp1[,8]))
means
covariance=cov(rdgp1)
covariance
#History 1
hist1=mvrnorm(216,mn,draw)
hist1

#RiskPort
library(RiskPortfolios)
a=covEstimation(hist1,control=list("bs"))
a
b=meanEstimation(hist1,control=list("bs"))
b

optimalPortfolio(a,b,control=list("bs"))

hist.returns=hist1
hist.means=apply(hist.returns,2,mean,na.rm=TRUE)
hist.cov=cov(hist.returns,use="complete.obs")
hist.cov
#final-opt.r
# Portfolio Optimizer
make.port <- function(means, 
                      covariance, 
                      max.allocation = .5,
                      min.allocation = 0,
                      risk.premium.up = .5, 
                      risk.increment=.005) {
  # A function to optimize a portfolio, using quadratic programming
  # Inputs:
  #       means           a vector of stock means
  #       covariance      the stocks' covariance matrix
  #       max.allocation  a parameter that controls the maximum weight in the final portfolio
  #       min.allocation  a parameter that controls the minimum weight in the final portfolio 
  #       risk.premium.up the maximum risk taken
  #       risk.increment  the step size
  # Outputs:
  #   The stock portfolio's efficent frontier 
  
  # Get the number of stocks
  n <- ncol(covariance)
  
  # Matrices to solve with quadprog, with max and min constraints
  Amat <- cbind(1, diag(n), -diag(n))
  bvec <- c(1, rep(min.allocation, n), rep(-max.allocation, n))
  meq <- 1
  
  # Set up the iterator
  iters <- seq(from=0, to=risk.premium.up, by=risk.increment) %>% setNames(nm = .)
  dvec <- lapply(iters, function(i) means * i)
  
  # Get the solutions
  sols <- lapply(dvec, solve.QP, Dmat = covariance, Amat = Amat, bvec = bvec, meq = meq)
  
  # Summarize and return the results
  ######################################
  # Generate a list of summary functions
  funs <- list("Std.Dev" = function(x) sqrt(x %*%  covariance %*% x),
               "Exp.Return" = function(x) x %*% means,
               "sharpe" = function(x) x %*% means / 
                 sqrt(x %*% covariance %*% x))
  
  # Apply to the solutions list
  summ <- lapply(funs, function(f) lapply(sols, 
                                          function(s) f(s$solution))) %>% 
    lapply(as.numeric) %>%
    do.call(what = "cbind")
  
  # Combine into a results data frame
  results <- lapply(sols, function(s) s$solution) %>% 
    do.call(what = "rbind") %>% 
    data.frame(summ)
  
  return(results)
}
#plot port and efficient frontier
# Plot port and efficient Frontier
eff.plot <- function(eff, eff.optimal.point) {
  p <- ggplot(eff, aes(x=Std.Dev, y=Exp.Return)) + 
    geom_line(col="Red") +
    geom_point(data=eff.optimal.point, aes(x=Std.Dev, y=Exp.Return, label=sharpe),
               size=5) +
    annotate(geom="text", x=eff.optimal.point$Std.Dev,
             y=eff.optimal.point$Exp.Return,
             label=paste("Risk: ",
                         round(eff.optimal.point$Std.Dev*100, digits=3),"\nReturn: ",
                         round(eff.optimal.point$Exp.Return*100, digits=4),"%\nSharpe: ",
                         round(eff.optimal.point$sharpe*100, digits=2), "%", sep=""),
             hjust=0, vjust=1.2) +
    ggtitle("Efficient Frontier and Optimal Portfolio") +
    labs(x = "Risk (standard deviation of portfolio)", y = "Return")
  return(p)
}

baseline_cs <- function(x) {
  # Takes a cumulative sum, with the beginning indexed at 1
  x[1] <- x[1] + 1
  cumsum(x)
}

# Portfolio performance
port.performance <- function(eff.optimal.point, test_port = port.test) {
  n <- dim(eff.optimal.point)[2] - 3
  weights <- eff.optimal.point[1,1:n] %>% as.numeric()
  tret <- test_port %>% extract2("R") %>% `%*%`(weights) %>% rev() %>% baseline_cs
  tlab <- sort(as.Date(row.names(port.test$R)),decreasing=FALSE)
  return(list(tret = tret, tlab = tlab, weights = weights))
}

# Quick summary
port.summary <- function(performance, test_port = port.test) {
  c(52 * mean(test_port$R %*% performance$weights), 
    sqrt(52) * sd(test_port$R %*%  performance$weights),
    52 * mean(test_port$R %*%  performance$weights)/ 
      (sqrt(52) * sd(test_port$R%*% performance$weights)))
}

# Descriptive Stats last six months and test portfolio
end <- 1
start <- end + 23

ss <- port.hist %>%
  extract2("R") %>% 
  extract((start:end), ) 

funs <- c("Asset", "Means", "StDev", "Total", "N")
c_summ_hist <- adply(ss, 2, each(mean, sd, sum, length)) %>% setNames(funs)
attach(c_summ_hist)
c_summ_new <- port.test$R %>% adply(2, each(mean, sd, sum, length)) %>% setNames(funs)
Cov <- cov(ss)
Cov.New <- cov(port.test$R)

# Plot of stock returns
##############################
# Get the market returns
market <- baseline.hist %>% extract2("R") %>% extract(start:end,) %>% c(1,.) %>% cumsum()

# Get the quote dates
hist.dates <- port.hist %>% extract2("R") %>% rownames() %>% extract((start + 1) : end) %>%
  as.Date()

# Get the portfolio returs
return.out <- port.hist %>% extract2("R") %>% 
  extract(start:end,) %>% rbind(1, .) %>% apply(2, cumsum)

# Put together for plotting
hist.plot <- data.frame(SP500 = market, return.out, Date=hist.dates) %>% melt(id = "Date")

# Make Plot
p.returns <- ggplot(hist.plot,aes(x=Date,y=value)) +
  geom_line(aes(group=variable,colour=variable)) + 
  ggtitle("Cumulative Returns, Previous 3 Months") +
  labs(y="Cumulative Return") + 
  scale_colour_discrete(name="Ticker Name") +
  theme(legend.position="bottom") + 
  guides(col = guide_legend(nrow = 2)) 

print(p.returns)

# Basic Stats, historical
stocks.table <- data.frame("Annualized Return" =  52 * c_summ_hist$Means * 100, 
                           "Annualized St. Dev." = sqrt(52) * c_summ_hist$StDev * 100, 
                           "Sharpe Ratio" = 52 * c_summ_hist$Means /(sqrt(52)* c_summ_hist$StDev),
                           check.names = FALSE, row.names = c_summ_hist$Asset) %>% 
  t

# Basic Stats, Test port
# Looking forward to get the means
test.table = data.frame("Annualized Return" =  52 * c_summ_new$Means * 100, 
                        "Annualized St. Dev." = sqrt(52) * c_summ_new$StDev * 100,
                        "Sharpe Ratio" = 52 * c_summ_new$Means /(sqrt(52)* c_summ_new$StDev),
                        check.names = FALSE, row.names = c_summ_new$Asset) %>% 
  t

# Make Naive Port
eff <- make.port(Means, Cov,.5,-.5)
eff.optimal.point <- subset(eff, sharpe == max(sharpe))

# Efficient Frontier and performance
p.eff.naive <- eff.plot(eff, eff.optimal.point)
print(p.eff.naive)
perf <- port.performance(eff.optimal.point)

# Assessing performance
Basic <- port.summary(perf)
# Looking forward to get the means
New.Means <- apply(port.test$R,2,mean)
New.Cov <- cov(port.test$R)
new.eff <- make.port(New.Means, New.Cov,.5,-.5)
new.eff.optimal.point <- subset(new.eff, sharpe == max(sharpe))

nperf <- port.performance(new.eff.optimal.point)
p.foresight <- eff.plot(eff, new.eff.optimal.point)
print(p.foresight)

# Summarizing for final tables
Market <- c(52*mean(baseline.test$R), sqrt(52)*sd(baseline.test$R), 
            52*mean(baseline.test$R)/(sqrt(52)*sd(baseline.test$R)))

Foresight <- port.summary(nperf)
# Bayesian Mean estimate, Gibbs Sampler, Known Variance
m = start-end
mh = as.numeric(dim(hist.returns)[1]) # Adjustment to prior cov
l0 = 1/mh*hist.cov                      # Historical variance
m0 = apply(hist.returns,2,mean)  # historical means
m0 = m0 + abs(m0*.2)             # pos adj. -- Analyst predicts growth
s = 1/m*Cov                       # sample Cov for mu
xbar = Means                    # sample means

mun = solve(solve(l0) + m*solve(s)) %*% (solve(l0) %*% m0 + m*solve(s) %*% xbar)
ln = solve(solve(l0) + m*solve(s))

gibbs.known.norm = function(burn.in, samp.size,mun, ln)
{
  k = length(mun)
  mu.sample = matrix(0,nrow=samp.size,ncol=k)
  mus = rep(0,k)
  betas = NULL
  for(i in 1:(burn.in + samp.size))
  {
    for(j in 1:k)
    {
      betas = ln[j,-j] %*% solve(ln[-j,-j])
      mj = mun[j] + betas %*% (mus[-j] - mun[-j])
      lj = ln[j,j] - ln[j,-j] %*% solve(ln[-j,-j]) %*% ln[-j,j]
      
      mus[j] = rnorm(1,mj,sqrt(lj))
    }
    if(i>burn.in)
    {
      mu.sample[i-burn.in,] = mus
    }
  }
  return(list(sample=mu.sample, mus=apply(mu.sample,2,mean)))
}

Bayes.Means = gibbs.known.norm(burn.in=5000,samp.size=10000,mun,ln)
b1.eff = make.port(Bayes.Means$mus, Cov,.5,-.5)
b1.optimal.point = b1.eff[b1.eff$sharpe==max(b1.eff$sharpe),]

# Efficient Frontier and performance
p = eff.plot(b1.eff, b1.optimal.point); plot(p)
b1.perf = port.performance(b1.optimal.point)
Bayes.Mean =port.summary(b1.perf) 

# Sum of squares matrix
S = matrix(0,nrow=dim(port.hist$R)[2],ncol=dim(port.hist$R)[2])

for(i in end:start)
{
  S = S + (port.hist$R[i,]-Means) %*% t(port.hist$R[i,]-Means)
}

# A function for the gibbs sampler
gibbs.unknown.ni = function(burn.in, samp.size,S,m,xbar)
{
  mus = rep(0,length(xbar))
  sigmas = matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  mu.sample = matrix(0,nrow=samp.size,ncol=length(xbar))
  sigma.sample = matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  for(i in 1:(burn.in + samp.size))
  {
    sigmas = riwish(m-1,S)
    mus = mvrnorm(n=1,mu=xbar,Sigma=sigmas/m)
    
    if(i>burn.in)
    {
      mu.sample[i-burn.in,] = mus
      sigma.sample = sigma.sample + sigmas
    }
  }
  return(list(mu.sample=mu.sample, mus=apply(mu.sample,2,mean),
              sigmas=1/samp.size * sigma.sample))
}

# Run the sampler
Bayes.unknown.ni = gibbs.unknown.ni(burn.in=5000,samp.size=10000,S,m,Means)

# Prepping for the portfolio optimization
cov.out = Bayes.unknown.ni$sigmas
row.names(cov.out) = port
colnames(cov.out) = port
b2.eff = make.port(Bayes.unknown.ni$mus, cov.out,.5,-.5)
b2.optimal.point = b2.eff[b2.eff$sharpe==max(b2.eff$sharpe),]

# Efficient Frontier and performance
p = eff.plot(b2.eff, b2.optimal.point); plot(p)
b2.perf = port.performance(b2.optimal.point)
Bayes.ni = port.summary(b2.perf)

# Both Unknown, Informative Prior
k0 = mh
m = start-end
v0 = k0-1
# m0 = apply(hist.returns,2,mean)  # historical means
# m0 = m0 + abs(m0*.2)             # pos adj. -- Analyst predicts growth
# l0 = 1/mh*hist.cov  # Historical variance
# s = 1/m*Cov                       # sample Cov for mu
# xbar = Means                    # sample means 

mun = k0/(k0 + m) * m0 + m/(k0 + m) *xbar
kn = k0 + m
vn = v0 + m
ln = l0 + S + k0 * m/(k0 + m)* (xbar - m0) %*% t(xbar - m0)

gibbs.unknown.ip = function(burn.in, samp.size,vn,ln,mun,kn)
{
  mus = rep(0,length(mun))
  sigmas = matrix(0, nrow=dim(ln)[1], ncol=dim(ln)[2])
  mu.sample = matrix(0,nrow=samp.size,ncol=length(mun))
  sigma.sample = matrix(0, nrow=dim(ln)[1], ncol=dim(ln)[2])
  
  for(i in 1:(burn.in + samp.size))
  {
    sigmas = riwish(vn,solve(ln))
    mus = mvrnorm(n=1,mu=mun,Sigma=sigmas/kn)
    
    if(i>burn.in)
    {
      mu.sample[i-burn.in,] = mus
      sigma.sample = sigma.sample + sigmas
    }
  }
  return(list(mu.sample=mu.sample, mus=apply(mu.sample,2,mean),
              sigmas=1/samp.size * sigma.sample))
}

Bayes.unknown.ip = gibbs.unknown.ip(burn.in=5000, samp.size=10000,vn,ln,mun,kn)

cov.out = Bayes.unknown.ip$sigmas
row.names(cov.out) = port
colnames(cov.out) = port
b3.eff = make.port(Bayes.unknown.ip$mus, cov.out,.5,-.5)
b3.optimal.point = b3.eff[b3.eff$sharpe==max(b3.eff$sharpe),]

# Efficient Frontier and performance
p = eff.plot(b3.eff, b3.optimal.point); print(p)
b3.perf = port.performance(b3.optimal.point)
Bayes.ip = port.summary(b3.perf)

# Final Comparison
compare = rbind(Market=Market, 
                Basic= Basic, 
                Foresight = Foresight, 
                "Unkown Mean" = Bayes.Mean,
                "Noninformative Prior" = Bayes.ni,
                "Informative Prior" = Bayes.ip)

compare = round(100*compare,2)
colnames(compare) = c("Expected Return", "Risk", "Annualized Sharpe")

# Comparing Weights
final.weights = rbind(Baseline= eff.optimal.point, 
                      Foresight = new.eff.optimal.point, 
                      "Unkown Mean" = b1.optimal.point,
                      "Noninformative Prior" = b2.optimal.point,
                      "Informative Prior" = b3.optimal.point)

# Cleaning up
final.weights = round(100*final.weights[,1:8],2)

# Simulated probabilities
market = rev(baseline.test$R)
market[1] = baseline.test$R[1]+1
market = cumsum(market)
final.returns = data.frame(Market = market,
                           b1 = b1.perf$tret,
                           b2 = b2.perf$tret,
                           b3 = b3.perf$tret,
                           naive = perf$tret,
                           Foresight = nperf$tret, 
                           Date = rev(as.Date(row.names((baseline.test$R)))))

ret.out = melt(final.returns, id.vars="Date")

p.cumret = ggplot(ret.out, aes(x=Date,y=value))
p.cumret = p.cumret + geom_line(aes(group=variable, colour =variable)) + 
  scale_colour_discrete(name="Test Portfolio",
                        labels=c("Market","Unknown Mean", "Noninformative Prior",
                                 "Informative Prior", "Baseline","Foresight")) +
  ggtitle("Cumulative Returns for Each Tested Portfolio") +
  labs(y="Cumulative Return") +
  theme(legend.position="bottom") + guides(col = guide_legend(nrow = 2)) +
  labs(y="Cumulative Return")
print(p.cumret)

# Simulating distributions of returns calculating probabilities
# Weights-matrix
w.matrix = t(rbind(eff.optimal.point[1:8],
                   new.eff.optimal.point[1:8],
                   b1.optimal.point[1:8],
                   b2.optimal.point[1:8],
                   b3.optimal.point[1:8]))

colnames(w.matrix) = row.names(final.weights)

sim.means = c(New.Means, Market=mean(baseline.test$R))
sim.Cov = cov(cbind(port.test$R, Market=baseline.test$R))

# 5000 Samples
n = 5000
m = dim(port.test$R)[1]
out= matrix(0, nrow = n, ncol = dim(w.matrix)[2]+1)
contribution = matrix(0, nrow = n * dim(w.matrix)[1], ncol = dim(w.matrix)[2])

for(i in 1:n)
{
  sample = rmvnorm(m,sim.means,sim.Cov)
  market = sample[,9]
  market[1] = market[1] + 1
  
  # Contribution
  cmat = matrix(0, nrow= dim(w.matrix)[1], ncol = dim(w.matrix)[2])
  
  for(j in 1: dim(w.matrix)[2])
  {
    cmat[,j] =  diag(w.matrix[,j]) %*% t(sample[,-9]) %*% 
      +     matrix(1, nrow = m, ncol=1)
  }
  
  # ready for output
  ports = sample[,-9] %*% w.matrix
  tot = apply(ports,2,sum)
  
  contribution[((i-1)*dim(w.matrix)[1] + 1):(i*dim(w.matrix)[1]),] = 
    cmat %*% solve(diag(tot))
  out[i,] = c(tot + 1,Market=sum(market))
}

# Add an index to contribution

colnames(contribution) = colnames(w.matrix)
idx = rep(row.names(w.matrix),n)

cont.means = aggregate(contribution,list(var=idx),mean)
row.names(cont.means) = cont.means[,1]
cont.means = cont.means[,-1]

cont.sds = aggregate(contribution,list(var=idx),sd)
row.names(cont.sds) = cont.sds[,1]
cont.sds = cont.sds[,-1]

# Sample Stats
sim.r.means = apply(out,2,mean)
names(sim.r.means) = c(colnames(w.matrix),"Market")
sim.r.cov = cov(out)
colnames(sim.r.cov) = c(colnames(w.matrix),"Market")
row.names(sim.r.cov) = c(colnames(w.matrix),"Market")
sd.r = sqrt(diag(sim.r.cov))

# Calculate Probs
# P(X > Y) = P(Y - X < 0)
probs = rep(0, dim(out)[2]-1)
names(probs) = c(colnames(w.matrix))
for(i in 1:(dim(out)[2]-1))
{
  diff = c(-1,1)
  d = c(i,6)
  m = sim.r.means[d]
  v = sim.r.cov[d,d]
  m1 = m %*% diff
  v1 = diff %*% v %*% diff
  probs[i] = pnorm(0,mean=m1,sd=sqrt(v1))
}

# Validation results
validate.table = cbind(Return = sim.r.means, Risk = sd.r, "P > Market " = c(probs,NA))
# Distribution plots
colnames(out) = c(colnames(w.matrix),"Market")
dist.plots = melt(out)
df.out = as.data.frame(out)
name = colnames(out)[1]

p = ggplot(df.out, aes(y=Baseline, x=Market)) + geom_point(alpha=.5) +
  geom_density2d() + annotate("segment", x=0.75,y=0.75,xend=2,yend=2, 
                              colour="blue")
print(p)

plots = dist.plots[dist.plots$Var2 %in% c(name,"Market"),]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_discrete(name="Portfolio") + 
  ggtitle(paste0("Densities for Simulated Portfolios, 104 Weeks in Test Period\n",
                 "P(>Market = ", round(probs[1],4),")"))
print(p)

plots = dist.plots[dist.plots$Var2 %in% name,]
plots[,3] = dist.plots[dist.plots$Var2 %in% "Market",3] -
  dist.plots[dist.plots$Var2 %in% name,3]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_manual(values="blue",guide=FALSE) + 
  ggtitle("Difference in Densities")
print(p)

name = colnames(out)[2]
plots = dist.plots[dist.plots$Var2 %in% c(name,"Market"),]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_discrete(name="Portfolio") + 
  ggtitle("Densities for Simulated Portfolios, 104 Weeks in Test Period")
print(p)

plots = dist.plots[dist.plots$Var2 %in% name,]
plots[,3] = dist.plots[dist.plots$Var2 %in% "Market",3] -
  dist.plots[dist.plots$Var2 %in% name,3]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_manual(values="blue",name="Portfolio") + 
  ggtitle("Difference in Densities")
print(p)

name = colnames(out)[3]
plots = dist.plots[dist.plots$Var2 %in% c(name,"Market"),]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_discrete(name="Portfolio") + 
  ggtitle("Densities for Simulated Portfolios, 104 Weeks in Test Period")
print(p)

plots = dist.plots[dist.plots$Var2 %in% name,]
plots[,3] = dist.plots[dist.plots$Var2 %in% "Market",3] -
  dist.plots[dist.plots$Var2 %in% name,3]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_manual(values="blue",name="Portfolio") + 
  ggtitle("Difference in Densities")
print(p)

name = colnames(out)[4]
plots = dist.plots[dist.plots$Var2 %in% c(name,"Market"),]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_discrete(name="Portfolio") + 
  ggtitle("Densities for Simulated Portfolios, 104 Weeks in Test Period")
print(p)

plots = dist.plots[dist.plots$Var2 %in% name,]
plots[,3] = dist.plots[dist.plots$Var2 %in% "Market",3] -
  dist.plots[dist.plots$Var2 %in% name,3]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_manual(values="blue",name="Portfolio") + 
  ggtitle("Difference in Densities")
print(p)

colnames(df.out)[5] = "IP"
p.ipdist = ggplot(df.out, aes(y=IP, x=Market)) + geom_point(alpha=.5) + geom_density2d() + 
  annotate("segment", x=0.75,y=0.75,xend=2,yend=2, colour="blue") +
  ylab("Bayesian, Informative Prior") + 
  ggtitle("Joint Distribution") + 
  theme(legend.position="bottom")
print(p.ipdist)

name = colnames(out)[5]
plots = dist.plots[dist.plots$Var2 %in% c(name,"Market"),]
p.ipdens = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_discrete(name="Portfolio") + 
  ggtitle("Densities") + 
  theme(legend.position="bottom")
print(p.ipdens)

plots = dist.plots[dist.plots$Var2 %in% name,]
plots[,3] = dist.plots[dist.plots$Var2 %in% "Market",3] -
  dist.plots[dist.plots$Var2 %in% name,3]
p = ggplot(plots, aes(x=value, fill=Var2)) + geom_density(alpha=.3) +
  scale_fill_manual(values="blue",name="Portfolio") + 
  ggtitle("Difference in Densities")
print(p)