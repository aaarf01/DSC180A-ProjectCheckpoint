# Fig 1C partial dependence at the RF used in Fig 1B. args: data dir, out dir
# Paper does not state how bars were built; this reports PD at each biomarker's own 95th percentile,
# minus PD at its median (LOG A11).
library(randomForest); library(pdp)
a<-commandArgs(TRUE); x<-read.csv(file.path(a[1],"data_insitu_v2_25-04.csv"))
y<-read.csv(file.path(a[1],"genomic_OCmatchups_v3_medians_2d.csv"))[,10];k<-!is.na(y)
P<-data.frame(nitri=x$insitu_and_woa_nutricline[k],Omega_P_High=x$Omega_P_High[k],Omega_N_High=x$Omega_N_High[k],Omega_N_Med=x$Omega_N_Med[k],Omega_Fe_High=x$Omega_Fe_High[k],Omega_Fe_Med=x$Omega_Fe_Med[k])
set.seed(123); m<-randomForest(P,y[k],importance=TRUE); out<-NULL
for(v in names(P)[-1]){q<-quantile(P[[v]],c(.5,.95)); g<-data.frame(q); names(g)<-v
 pd<-partial(m,pred.var=v,pred.grid=g,train=P); out<-rbind(out,data.frame(var=v,x_med=q[1],x_95=q[2],pd_med=pd$yhat[1],pd_95=pd$yhat[2]))}
write.csv(out,file.path(a[2],"fig1C_pdp.csv"),row.names=FALSE)
write.csv(data.frame(importance(m)),file.path(a[2],"rf_importance.csv"))
for(v in names(P)[-1]){pd<-partial(m,pred.var=v,grid.resolution=25,train=P);write.csv(pd,file.path(a[2],paste0("pd_",v,".csv")),row.names=FALSE)}
