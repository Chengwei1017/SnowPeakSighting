#!/usr/bin/env Rscript
print("MCSNOW V20220717")
print("Monte Carlo Simulation for Snow peak sighting forecast model.")
library(reshape2)
library(dplyr)
library(doParallel)
if(!require(Rmpi)) rMPI=F else rMPI=T
HPC=T
USERANID=T
if(HPC){
SNOWWDR="/public/home/cdaes_luchengwei/LCW/Models/MCSNOW/"
OUTDIR="/public/home/cdaes_luchengwei/LCW/Models/MCSNOW/"
}else{
SNOWWDR="/home/lcw/LCW/Models/MCSNOW/"
OUTDIR="/home/lcw/LCW/Models/MCSNOW/"
}
SNOWMONT=c("01","02","03","04","05","06","07","08","09","10","11","12")
SNOWFILE=paste0(SNOWWDR,"原始数据/2020",SNOWMONT,"2/SNOW/GRSNOWPROFILE_YMF.csv")
SNOWSHAD=read.csv(file=paste0(SNOWWDR,"SNOWSHADE.csv"),header=T)[,c(3,4,7)]
SIMTIME=100000
RUNCORE=30
TESTONLY=F

if(rMPI){
  rID=mpi.comm.rank(0)
  rLN=mpi.comm.size(0)
  if(USERANID) rID=round(runif(1,1000,9999),0)
  print("Using MPI...")
  print(paste0("Job ID:",rID))
  print(paste0("Job LN:",rLN))

  if(rID==0){
    SIMTIME=SIMTIME/rLN
  }
  SIMTIME=mpi.bcast(SIMTIME,1,0,0)
  FTAG=paste0(Sys.Date(),"_",substr(Sys.time(),12,13),"_R",SIMTIME,"_J",rID)
  print(paste0("Run times for this thread:",SIMTIME,"."))
}else{
  FTAG=paste0(Sys.Date(),"_",substr(Sys.time(),12,13),"_R",SIMTIME)
}

print(paste0("Processing ",FTAG,"..."))

count_true2=function(PK,RS){
  return(length(PK[PK==T])/length(RS[RS==T]))
}

acc_calc=function(daily_rep,limits=0.3){
  daily_for=daily_rep
  SPDAY=read.csv(file = paste0(SNOWWDR,"YMF.csv"),header=T)
  daily_for$ifPK=0
  daily_for$ifPK[daily_for$Rate>limits]=1
  dailycmp=data.frame(Date=daily_for$DATE,FOR=daily_for$ifPK)
  dailycmp=base::merge(dailycmp,SPDAY,by.x="Date",by.y="日期",all.x=T)
  dailycmp$可见幺妹峰[is.na(dailycmp$可见幺妹峰)]=0
  TP=nrow(subset(dailycmp,可见幺妹峰==1 & FOR==1))
  TN=nrow(subset(dailycmp,可见幺妹峰==0 & FOR==0))
  FP=nrow(subset(dailycmp,可见幺妹峰==0 & FOR==1))
  FN=nrow(subset(dailycmp,可见幺妹峰==1 & FOR==0))
  TD=nrow(dailycmp)
  Pr=TP/(TP+FP)*100
  Ac=(TP+TN)/(TP+FP+FN+TN)*100
  Re=TP/(TP+FN)*100
  F1=2*(Pr*Re)/(Pr+Re)
  return(data.frame(TP,TN,FP,FN,Pr,Ac,Re,F1))
}

seq_comp=function(daily_rep,limits=0.3){
  spst=subset(daily_rep,Rate>limits)
  daily_for=daily_rep
  daily_for$ifPK=0
  daily_for$ifPK[daily_for$Rate>limits]=1
  SPDAY=read.csv(file = paste0(SNOWWDR,"YMF.csv"),header=T)
  dailycmp=data.frame(Date=daily_for$DATE,FOR=daily_for$ifPK)
  dailycmp=base::merge(dailycmp,SPDAY,by.x="Date",by.y="日期",all.x=T)
  dailycmp$可见幺妹峰[is.na(dailycmp$可见幺妹峰)]=0
  dailycmp$Month=substr(dailycmp$Date,6,7)
  names(dailycmp)=c("Date","Modeled","Posted","Month")
  dcp=reshape2::melt(dailycmp,id=c("Date","Month"))
  
  forln=nrow(dailycmp)
  hitdf=subset(dailycmp,Posted==1)
  repdf=subset(dailycmp,Modeled==1 & Posted==0)
  
  rept=data.frame(
    hitrate=length(hitdf$Modeled[hitdf$Modeled==1])/nrow(hitdf),
    miscnt=length(repdf$Modeled[hitdf$Modeled==1])
  )
  return(rept)
}

daily_snowindex=function(GRIN,fixTopo=T,SPthreshold=0.3){
  cRate=function(rate,srate,sightrate_t=SPthreshold,SRT=0.5){
    iRep=data.frame(rate,srate)
    iRep=subset(iRep,srate<SRT)
    iCnt=length(iRep$rate[iRep$rate>=sightrate_t])/length(iRep$rate)
    return(iCnt)
  }
  SNOWRATE=as.data.frame(dplyr::summarise(dplyr::group_by(GRIN,across(all_of(c("LAT","LON","X","Y","Date")))),.groups="drop",rate=count_true2(ifPK,ifPRS0),count=n(),HoS=sum(PKH)))
  SNOWRATE=base::merge(SNOWRATE,SNOWSHAD,by=c("X","Y"),all.x=T)
  SNOWRATE$TOPOFIXED=F
  if(fixTopo){
    SNOWRATE$rate=SNOWRATE$rate*(1-SNOWRATE$ShadeRate)
    SNOWRATE$HoS=SNOWRATE$HoS*(1-SNOWRATE$ShadeRate)
  }
  
  drep=as.data.frame(dplyr::summarise(dplyr::group_by(SNOWRATE,Date),.groups="drop",Rate=cRate(rate,ShadeRate)))
  drep$SPthreshold=SPthreshold
  rm(SNOWRATE)
  names(drep)=c("DATE","Rate","SPthreshold")
  return(drep)
}

mcm_snow=function(parid){
  SR=determ_snowidx(GRPROF,MCMThres[parid,])
  daily_rep=daily_snowindex(SR)
  rm(SR)
  rep=seq_comp(daily_rep)
  sco=acc_calc(daily_rep)
  return(data.frame(PARID=parid,rep,sco))
}

gen_thres=function(ranlen=1000000){
  qcthread=runif(ranlen, 0, 1.5)
  pcldthread=runif(ranlen, 0, 100)
  qvthread=runif(ranlen, 0, 22)
  qvathread=runif(ranlen, 0, 15)
  qvsthread=runif(ranlen, 0, 526)
  pmthread=runif(ranlen, 0, 150)
  pmathread=runif(ranlen, 0, 30)
  pmmthread=runif(ranlen, 0, 170)
  rsthread=runif(ranlen, 0, 370)
  thres_df=data.frame(ID=c(1:ranlen),qcthread,pcldthread,qvthread,qvathread,qvsthread,pmthread,pmathread,pmmthread,rsthread)
  return(thres_df)
}

determ_snowidx=function(snowprofile,thres,LRS=F){
  qcthread=thres$qcthread  
  pcldthread=thres$pcldthread
  qvthread=thres$qvthread  
  qvathread=thres$qvathread 
  qvsthread=thres$qvsthread 
  pmthread=thres$pmthread  
  pmathread=thres$pmathread 
  pmmthread=thres$pmmthread 
  rsthread=thres$rsthread 
  tsnowidx=data.frame(Time=snowprofile$Time,Date=substr(snowprofile$Time,1,10),MONTH=substr(snowprofile$Time,6,7),LAT=snowprofile$LAT,LON=snowprofile$LON,X=snowprofile$X,Y=snowprofile$Y,ifQC=snowprofile$QCS<=qcthread,ifCLD=snowprofile$peakcl<=pcldthread,ifQV=snowprofile$QVS<=qvthread,ifQVM=snowprofile$QVM<=qvathread,ifQVS=snowprofile$QVSum<=qvsthread,ifPM=snowprofile$localpm<pmthread,ifPMA=snowprofile$PMA<pmathread,ifPMM=snowprofile$PMM<pmmthread,ifRS=snowprofile$peakrs>=rsthread,ifPRS0=snowprofile$peakrs>0,ifLRS=snowprofile$localrs<=snowprofile$peakrs)
  if(LRS){
  tsnowidx$ifPK=(tsnowidx$ifQC & tsnowidx$ifQV & tsnowidx$ifPM & tsnowidx$ifRS & tsnowidx$ifLRS & tsnowidx$ifCLD & tsnowidx$ifQVM & tsnowidx$ifPMA & tsnowidx$ifPMM)
  }else{
  tsnowidx$ifPK=(tsnowidx$ifQC & tsnowidx$ifQV & tsnowidx$ifPM & tsnowidx$ifRS & tsnowidx$ifCLD & tsnowidx$ifQVM & tsnowidx$ifPMA & tsnowidx$ifPMM)
  }
  tsnowidx$PKH=0
  tsnowidx$PKH[tsnowidx$ifPK]=1
  return(tsnowidx)
}

crt_thres=data.frame(pmthread=50,
                     pmathread=30,
                     pmmthread=50,
                     rsthread=50,
                     shadthread=100,
                     pcldthread=40,
                     qvthread=30,
                     qvathread=25,
                     qvsthread=200,
                     qcthread=0)

if(file.exists("RAW.RData")){
print("Reading raw data...")
load("RAW.RData")
}else{
GRPROF=data.frame()
pb=txtProgressBar(max=length(SNOWFILE),style=3)
for(xi in 1:length(SNOWFILE)){
  GRDF=data.frame(MONTH=SNOWMONT[xi],read.csv(file=SNOWFILE[xi],header=T))
  GRPROF=rbind(GRPROF,GRDF)
  setTxtProgressBar(pb,xi)
}
close(pb)
save(GRPROF,file="RAW.RData")
}

MCMThres=gen_thres(SIMTIME)

if(!TESTONLY){
registerDoParallel(cores=RUNCORE)
pres=foreach(
  i=1:nrow(MCMThres), 
  .combine=rbind, 
  .packages = c("dplyr"),
  .export=c("GRPROF","mcm_snow"),
  .verbose=T
) %dopar% mcm_snow(i)

write.csv(file=paste0(OUTDIR,"RESULTS/MCM_RESULTS_",FTAG,".csv"),pres)
write.csv(file=paste0(OUTDIR,"RESULTS/MCM_THRESVL_",FTAG,".csv"),MCMThres)
}else{
print("Single job test...")
system.time({
SR=determ_snowidx(GRPROF,crt_thres)
daily_rep=daily_snowindex(SR)
rm(SR)
rep=seq_comp(daily_rep)
sco=acc_calc(daily_rep)
print(data.frame(rep,sco))
})
}

if(rMPI) mpi.quit()
