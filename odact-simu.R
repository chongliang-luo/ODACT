
## ODACT simulation study (JAMIA revision 20241005)  
## add ODACT with small leading site, and ODACT-avg, and tabulate bias and MSE

require(RcppArmadillo)
require(Rcpp)
library(survival)
library(mvtnorm)
library(beeswarm)
library(vioplot)
library(data.table)
library(ggplot2)

source('odact-engine.R')

################################## SIMULATION #################################
## setup parameters
# rm(list=ls())
n <- 1000
K <- 10    
px <- 2
h <- 0.2
cens <- 15  # cens smaller, more events 
evalt <- seq(0, 0.8, 0.2)
nt = length(evalt)
c0 = -2
nl <- n/K 
n.site <- c(rep(nl*1.5, K/2), rep(nl*0.5, K/2)) # half/half big/small sites
# id.site = sample(rep(1:K, n.site))

## simulate data
# cens=beta shape2, larger => censor smaller => event rate lower
dd = data.table(coxlcgen(n=n, cens=cens, c0=-2) ) 
plot(survfit(Surv(yo, do)~1, data=data.frame(dd)))
plot(density(dd$yo[dd$do==1])) 
# event rate at each evalt 
sapply(evalt,  function(a) sum(dd$do[dd$yo< a+h & dd$yo>= a-h])/n)


## run replicates
nrep = 50
res <- replicate(nrep, try(do.one.lc(n=n, K=K, n.site=n.site, cens=cens, c0=-2,
                                    evalt = evalt, h=h, init='local')), simplify=FALSE)
# res[[1]]$gold
# res[[1]]$locals
# res[[1]]$meta
# res[[1]]$sl2
# apply(res[[1]]$sl2, 1:2, mean)
 

# remove a few NA replicates
id.na = which(unlist(lapply(res, length)) ==1) 
length(id.na)
if(length(id.na) > 0)  res = res[-id.na]
 

## extract all estimates as a dataframe
ball = data.frame()
for(ix in 1:px){
  for(i in 1:length(evalt)){  # 0  to 0.8
    dp = c(
      do.call(rbind, lapply(res, function(x) x$gold[ix,]))[,i], 
      do.call(rbind, lapply(res, function(x) x$meta[ix,]))[,i],
      do.call(rbind, lapply(res, function(x) x$sl2[ix,,1]))[,i], 
      do.call(rbind, lapply(res, function(x) x$sl2[ix,,K]))[,i],
      do.call(rbind, lapply(res, function(x) apply(x$sl2,1:2,mean)[ix,]))[,i]  
    )
    ball = rbind(ball, 
                 data.frame(method=rep(c('pooled', 'meta', 
                                         'ODACT (large)', 'ODACT (small)', 
                                         'ODACT (avg)'), each=nrep-length(id.na)), 
                            ix=ix, it=i,X=paste0('X',ix), beta.t = c(-2,0)[ix]+c(1,0)[ix]*evalt[i],
                            time=paste0('time = ',evalt[i]), est=dp) )
  }
}
ball$method = factor(ball$method, levels = c('pooled', 'meta', 
                                             'ODACT (large)', 'ODACT (small)', 
                                             'ODACT (avg)'))
ball = data.table(ball)

## Figure 1, boxplots
# pdf('Fig1.pdf', height=9, width=12)
ggplot(ball, aes(x=method, y=est, color=method)) +
  geom_boxplot(stat = "boxplot", notch = TRUE) +
  facet_grid(X~time)+   
  geom_hline(data=ball, aes(yintercept = beta.t), linetype='dashed', col='red') +
  coord_cartesian(ylim = c(-4, 4))+ # ylim(-40, 40) +
  labs(title='', x='', y="log hazard ratio" ) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)
        , axis.ticks.x = element_blank() 
        , panel.background = element_blank()
        , panel.grid.minor = element_blank()
        , panel.grid.major = element_blank()
        # , panel.border = element_rect()
        , axis.title = element_text(size=14)
        , axis.text.y = element_text(size=14)
        , axis.text.x = element_blank()
        , strip.text.x = element_text(size=14)
        , strip.text.y = element_text(size=14, angle=0)
        , legend.position = c(0.08, 0.9) # 'none' # 
        # , legend.key = element_rect(size=0.3)
        # , legend.key.height = unit(0.5, "cm")
  )   + 
  scale_color_manual(values=c("orange",'black', "blue","lightblue","green"))
# dev.off()



## RMSE, and t.test of bias to the pooled est
# par(mfrow=c(2, nt), mar=c(3,2,2,0), mgp=c(1.25,0.5,0))
bias = data.frame()
for(ix in 1:px){
  for(i in 1:(nt)){
    dp = c(
      do.call(rbind, lapply(res, function(x) x$meta[ix,]-x$gold[ix,]))[,i],
      do.call(rbind, lapply(res, function(x) x$sl2[ix,, 1]-x$gold[ix,]))[,i],
      do.call(rbind, lapply(res, function(x) x$sl2[ix,,10]-x$gold[ix,]))[,i],
      do.call(rbind, lapply(res, function(x) apply(x$sl2,1:2,mean)[ix,]-x$gold[ix,]))[,i]  
    )
    bias = rbind(bias, 
                 data.frame(method=rep(c('meta', 'ODACT (large)', 'ODACT (small)', 'ODACT (avg)'), 
                                       each=nrep-length(id.na)), 
                           ix=ix, it=i,X=paste0('X',ix),  
                           time=paste0('time = ',evalt[i]), est=dp) ) 
  }
}

## Table S1
bias = data.table(bias)
bias[,MSE:=est^2]
sd.trim = function(a,n=c(0,20)) sd(a[order(a)[(n[1]+1):(length(a)-n[2])]])
mean.trim = function(a,n=c(0,25)) mean(a[order(a)[(n[1]+1):(length(a)-n[2])]]) 
tabs1 = bias[,.(RMSE=round(sqrt(mean.trim(MSE)),2),
                MAD=round(median(sqrt(MSE)),2),
                pv=round(t.test(est)$p.val,3)),by=.(X,method,time)]  
tabs1[,MAD.pv := paste0(MAD, ' (', pv, ')')]
tabs1 = cbind(tabs1$method[1:4], rbind(tabs1$time[tabs1$method=='meta'],matrix(tabs1$MAD.pv, nrow=4)))
# write.csv(tabs1, file='TableS1_rev.csv')
############################## END: SIMULATION #################################
