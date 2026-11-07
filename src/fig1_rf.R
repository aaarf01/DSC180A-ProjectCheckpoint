# Adapted from random_forest_analysis_v2.R by L. Ustick.  args: <data dir> <out dir>
# Runs the sat_data list both as written in the repo (slot 1 = 1d, LOG D2) and corrected (slot 1 = 1p).
library(randomForest)
a <- commandArgs(TRUE); d <- a[1]; o <- a[2]
x <- read.csv(file.path(d,"data_insitu_v2_25-04.csv"))
rd <- function(s) read.csv(file.path(d,sprintf("genomic_OCmatchups_v3_medians_%s.csv",s)))
fixed <- lapply(c("1p","3p","5p","1d","2d","5d"),rd); asrepo <- fixed; asrepo[[1]] <- fixed[[4]]
pred <- function(s) cbind(nitri=x$insitu_and_woa_nutricline[s],p_lim_h=x$Omega_P_High[s],n_lim_h=x$Omega_N_High[s],
  n_lim_m=x$Omega_N_Med[s],fe_lim_h=x$Omega_Fe_High[s],fe_lim_m=x$Omega_Fe_Med[s])
r2 <- function(m){n<-length(m$y);100*(1-m$mse[length(m$mse)]/(var(m$y)*(n-1)/n))}
grid <- function(L){G<-N<-matrix(NA,15,6);for(t in 1:15)for(s in 1:6){y<-L[[s]][,t];k<-!is.na(y);set.seed(123)
  G[t,s]<-r2(randomForest(pred(k),y[k]));N[t,s]<-sum(k)};list(G=G,N=N)}
lab <- list(paste0("pm",seq(2,30,2),"D"),c("1P","3x3P","5x5P","1x1D","2x2D","5x5D"))
for(nm in c("asrepo","fixed")){g<-grid(get(nm));dimnames(g$G)<-dimnames(g$N)<-lab
 write.csv(round(g$G,2),file.path(o,paste0("figS4_r2_",nm,".csv")));write.csv(g$N,file.path(o,paste0("figS4_n_",nm,".csv")))}
y<-fixed[[5]][,10];k<-!is.na(y);set.seed(123);m1<-randomForest(as.data.frame(pred(k)[,"nitri",drop=FALSE]),y[k],ntree=500)
set.seed(123);m2<-randomForest(pred(k),y[k])
write.csv(data.frame(model=c("nutricline","nutricline+genes"),r2=c(r2(m1),r2(m2))/100,n=sum(k)),file.path(o,"fig1B_r2.csv"),row.names=FALSE)
write.csv(data.frame(theta=y[k],nut_pred=predict(m1),full_pred=predict(m2),pred(k)),file.path(o,"rf_oob.csv"),row.names=FALSE)
