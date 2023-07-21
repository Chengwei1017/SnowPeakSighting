#!/usr/bin/env Rscript
options(warn=-1)
options(info=-1)

print("SNOWPEAK V20210907")
#V20210118 Parallel process.
#V20210222 CN grids process and UL forecasts.
#V20210508 Pre-estimation

today="2017-07-01"
CASENAME="00"
extGrid=F
CMAQDOMS="2"
ProHGrids=F
USEIDW=F
CFSF=F
GFSF=F
CMAF=F
CHNF=F
ULF=F

USECAPLOC=T
USEPRECALCSHADE=T
Args=commandArgs()
if(length(Args)>5){
  today=Args[6]
  CASENAME=Args[7]
  FORERUN=F
  SENARUN=F
  SENAFOR=F
  if(Args[8]=='FORE') FORERUN=T
  if(Args[8]=='SCEN') SENARUN=T
  if(Args[8]=="SCFO") SENAFOR=T
  
  if(Args[8]=="ULF")  ULF=T
  if(Args[8]=="CFSF") CFSF=T
  if(Args[8]=="GFSF") GFSF=T
  if(Args[8]=="CMAF") CMAF=T
  if(Args[8]=="CHNF") CHNF=T
  
  CMAQDOMS=Args[9]
  extGrid=F
  ProHGrids=F
}else{
  print("USAGE: EXTCMAQ STARTMON(yyyy-mm) OUTPUTTAG FORE/SCEN/NORM/CFSF/GFSF/CMAF/CHNF/ULF CMAQDOM USEIDW(I)/ProHGrids(H)*")
  quit()
}

if(length(Args)>9) {
  if(Args[10]=='H') ProHGrids=T
  if(Args[10]=='I') USEIDW=T
}


library(RNetCDF)
library(lcwcmaq)
library(ggplot2)
library(reshape2)
library(dplyr)
library(foreign)
library(lcwmars)
if(!require(doParallel)) PARAPRO=F else PARAPRO=T
if(USECAPLOC){
GetLoc=function(latdf,londf,latlon){
  x=length(londf[,1])
  y=length(latdf[1,])
  xp=0
  yp=0
  lat=latlon[1]
  lon=latlon[2]
  for(i in 1:(y-1)){
    for(j in 1:(x-1)){
      if(latdf[j,i]<lat && latdf[j,(i+1)]>lat && londf[j,i]<lon && londf[(j+1),i]>lon){
        xp=i
        yp=j
        break
      }
    }
  }
  return(c(yp,xp))
}
}

GetPK2D=function(locx,locy,vhgt,rhgt,pcshade){
    VH=ZFF[locx,locy,,1]
    chgt=vhgt-rhgt
    if(debug) print(paste(vhgt,rhgt,chgt,locx,locy))
    vid=0
    for(vi in 2:length(VH)){
      if(VH[vi]>chgt & VH[vi-1]<chgt) vid=vi
    }
    if(chgt<VH[1]) vid=1
    shade=0
    if(USEPRECALCSHADE){
      if(pcshade=="TRUE") shade=1
    }else{
    if(chgt<0 & abs(chgt)>shadthread ){
      shade=1
      print(paste0("Warning: Peak shading found for this point.(",chgt,")"))
    }
    }

    lengqc=length(QC[1,1,1,])
    if(debug) print(paste0("Vid for this point:",vid," over ",chgt,"..."))
    odf=data.frame(
      T2=TEMP2[locx,locy,],
      WD=WDIR10[locx,locy,],
      WS=WSPD10[locx,locy,],
      PBL=PBL[locx,locy,],
      PSFC=PRSFC[locx,locy,],
      RH2=RH[locx,locy,],
      PCPN=PCPN[locx,locy,],
      SR=RGRND[locx,locy,],
      CLD=CFRAC[locx,locy,],
      VIS=3.9/(10*exp(VISM[locx,locy,]/10)),
      QC=QC[locx,locy,vid,2:lengqc]*1000, #g/kg,
      QV=QV[locx,locy,vid,2:lengqc]*1000,
      PM25=PM25_UGM3[locx,locy,],
      ifShade=shade
     )
  
  todf=cbind(Time=Timeget(),odf)
  return(todf)
}
#######################################################################################
qxgrids=read.csv(file="/home/lcw/LCW/Models/RMODELS/SNOWINX/snowobs.csv",header=T)
pkgrids=read.csv(file="/home/lcw/LCW/Models/RMODELS/SNOWINX/snowpeak.csv",header=T)
#######################################################################################
APPL=paste0(substr(today,1,4),substr(today,6,7),"")
if(FORERUN) APPL="SCASTOBSD"
if(CFSF) APPL="CFSCSTOBSD"
if(GFSF) APPL="GFSCSTOBSD"
if(CMAF) APPL="CMACSTOBSD"
if(CHNF) APPL="CHNCSTOBSD"

if(FORERUN) {print("WARNING:FORERUN TRUE!")}
if(SENARUN) {print("WARNING:SENARUN TRUE!")}
if(SENAFOR) {print("WARNING:SENAFOR TRUE!")}

MFNAME=paste("//NAS//MODELOUT//CMAQ//mcip//",APPL,"OBSD",CMAQDOMS,"//METCRO2D",sep="")
GFNAME=paste("//NAS//MODELOUT//CMAQ//mcip//",APPL,"OBSD",CMAQDOMS,"//GRIDCRO2D",sep="")
GDNAME=paste("//NAS//MODELOUT//CMAQ//mcip//",APPL,"OBSD",CMAQDOMS,"//GRIDDOT2D",sep="")
GF3D2=paste("//NAS//MODELOUT//CMAQ//mcip//",APPL,"OBSD",CMAQDOMS,"//METCRO3D",sep="")
SMKNAME=paste("//NAS//MODELOUT//SMOKE//",APPL,"OBSD",CMAQDOMS,"//emis3d.ncf",sep="")

FNAME=paste("//NAS//MODELOUT//CMAQ//cctm//",APPL,"OBSD",CMAQDOMS,"//CCTM_P.",APPL,"OBSD",CMAQDOMS,".AIRVIEW",".ncf",sep="")
PANAME=paste("//NAS//MODELOUT//CMAQ//cctm//",APPL,"OBSD",CMAQDOMS,"//CCTM_P.",APPL,"OBSD",CMAQDOMS,".PA_1",".ncf",sep="")
SANAME=paste("//NAS//MODELOUT//CMAQ//cctm//",APPL,"OBSD",CMAQDOMS,"//CCTM_P.",APPL,"OBSD",CMAQDOMS,".SA_ACONC",".ncf",sep="")

if(SENARUN){
  OUTDIR=paste0("//NAS//MODELOUT//AIRCSV//CMAQSIMU//",APPL,"//") 
  FNAME= paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",APPL,"OBSD",CMAQDOMS,"_AIRVIEW")
  if(!file.exists(FNAME)) FNAME= paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",APPL,"OBSD","_AIRVIEW") #OLD RULE.
  PANAME=paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",APPL,"OBSD",CMAQDOMS,"_IPRCONC")
  SANAME=paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",APPL,"OBSD",CMAQDOMS,"_SAACONC")
  GFNAME=paste("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//GRIDCRO2D",CMAQDOMS,sep="")
  GDNAME=paste("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//GRIDDOT2D",CMAQDOMS,sep="")
}

if(SENAFOR){
  CASTTP=Sys.getenv("CASTTPYE")
  if(CASTTP==""){
    print("No CASTTYPE EV found, use GFS as default.")
    CASTTP="GFS"
  }
  print(paste0("Forecast type for this scenfore: ",CASTTP,"."))
  OUTDIR=paste0("//NAS//MODELOUT//AIRCSV//CMAQSIMU//",APPL,"//") 
  FNAME= paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",CASTTP,"CSTOBSD",CMAQDOMS,"_AIRVIEW")
  print(FNAME)
  if(!file.exists(FNAME)) FNAME= paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",CASTTP,"CSTOBSD","_AIRVIEW") #OLD RULE.
  PANAME=paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",APPL,"OBSD",CMAQDOMS,"_IPRCONC")
  SANAME=paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",APPL,"OBSD",CMAQDOMS,"_SAACONC")
  MFNAME=paste("//NAS//MODELOUT//CMAQ//mcip//",CASTTP,"CSTOBSD",CMAQDOMS,"//METCRO2D",sep="")
  GFNAME=paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//GRIDCRO2D",CMAQDOMS)
  GDNAME=paste("//NAS//MODELOUT//CMAQ//mcip//",CASTTP,"CSTOBSD",CMAQDOMS,"//GRIDDOT2D",sep="")
  GF3D2=paste("//NAS//MODELOUT//CMAQ//mcip//",CASTTP,"CSTOBSD",CMAQDOMS,"//METCRO3D",sep="")
  SMKNAME=paste("//NAS//MODELOUT//SMOKE//",CASTTP,"CSTOBSD",CMAQDOMS,"//emis3d.ncf",sep="")
}

if(FORERUN | CFSF | GFSF | CMAF | CHNF) {
  FNAME=paste("//NAS//MODELOUT//CMAQ//cctm//",APPL,CMAQDOMS,"//CCTM_P.",APPL,CMAQDOMS,".AIRVIEW",".ncf",sep="")
  PANAME=paste("//NAS//MODELOUT//CMAQ//cctm//",APPL,CMAQDOMS,"//CCTM_P.",APPL,CMAQDOMS,".PA_1",".ncf",sep="")
  SANAME=paste("//NAS//MODELOUT//CMAQ//cctm//",APPL,CMAQDOMS,"//CCTM_P.",APPL,CMAQDOMS,".SA_ACONC",".ncf",sep="")
  MFNAME=paste("//NAS//MODELOUT//CMAQ//mcip//",APPL,CMAQDOMS,"//METCRO2D",sep="")
  GFNAME=paste("//NAS//MODELOUT//CMAQ//mcip//",APPL,CMAQDOMS,"//GRIDCRO2D",sep="")
  GDNAME=paste("//NAS//MODELOUT//CMAQ//mcip//",APPL,CMAQDOMS,"//GRIDDOT2D",sep="")
  GF3D2=paste("//NAS//MODELOUT//CMAQ//mcip//",APPL,CMAQDOMS,"//METCRO3D",sep="")
  SMKNAME=paste("//NAS//MODELOUT//SMOKE//",APPL,CMAQDOMS,"//emis3d.ncf",sep="")
}

if(ULF){
  OUTDIR=paste0("//NAS//MODELOUT//AIRCSV//CMAQSIMU//",APPL,"//") 
  FNAME= paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//CFSCSTOBSD",CMAQDOMS,"_AIRVIEW")
  PANAME=paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",APPL,"OBSD",CMAQDOMS,"_IPRCONC")
  SANAME=paste0("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//",APPL,"OBSD",CMAQDOMS,"_SAACONC")
  MFNAME=paste("//NAS//MODELOUT//CMAQ//mcip//CFSCSTOBSD",CMAQDOMS,"//METCRO2D",sep="")
  GFNAME=paste("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//GRIDCRO2D",CMAQDOMS,sep="")
  GDNAME=paste("//NAS//MODELOUT//CMAQ//cctm//SENA//",CASENAME,"//GRIDDOT2D",CMAQDOMS,sep="")
  GF3D2=paste("//NAS//MODELOUT//CMAQ//mcip//CFSCSTOBSD",CMAQDOMS,"//METCRO3D",sep="")
  SMKNAME=paste("//NAS//MODELOUT//SMOKE//CFSCSTOBSD",CMAQDOMS,"//emis3d.ncf",sep="")
}


if(file.exists(PANAME)) PA=T
if(file.exists(SANAME)) ProSA=T


print("File checking..")
print(paste("FNAME:",FNAME,file.exists(FNAME)))
print(paste("MFNAME:",MFNAME,file.exists(MFNAME)))
print(paste("GFNAME:",GFNAME,file.exists(GFNAME)))
print(paste("GDNAME:",GDNAME,file.exists(GDNAME)))
print(paste("GF3D2:",GF3D2,file.exists(GF3D2)))
print(paste("SMKNAME:",SMKNAME,file.exists(SMKNAME)))
print("Done.")

MFNAME3D2=open.nc(GF3D2)
ZF=var.get.nc(MFNAME3D2,"ZF")[,,1,]
Sdata=length(ZF[1,1,])-1
TA=var.get.nc(MFNAME3D2,"TA")-273.15
ZFF=var.get.nc(MFNAME3D2,"ZF")
QC=var.get.nc(MFNAME3D2,"QC")
QV=var.get.nc(MFNAME3D2,"QV")


close.nc(MFNAME3D2)

aq=open.nc(FNAME)
lays=dim.inq.nc(aq,"LAY")$length
rows=dim.inq.nc(aq,"ROW")$length
cols=dim.inq.nc(aq,"COL")$length
cmaqtime=var.get.nc(aq,"TFLAG")[,,1][1,1]

today=as.Date(paste0(substr(cmaqtime,1,4),"-01-01"))+as.numeric(substr(cmaqtime,5,7))-1
print(paste0("Start Date from File:",today))
close.nc(aq)

jultoday=JulDate(today)
STARTTIME=FDate(paste(strptime(FDate(today), "%m/%d/%y %H:%M:%S")+9*3600,"",sep=""))
RDATE=paste("20",substr(STARTTIME,7,8),"-",substr(STARTTIME,1,2),"-",substr(STARTTIME,4,5),sep="")

if(SENARUN) APPL=paste0(substr(today,1,4),substr(today,6,7),"")

OUTDIR=paste0("//NAS//MODELOUT//AIRCSV//CMAQSIMU//",APPL,"//")
if(FORERUN | SENAFOR){
  OUTDIR=paste0("//NAS//MODELOUT//AIRCSV//CMAQSIMU//FORECAST//",today,"//")  
}

aq=open.nc(FNAME)
mt=open.nc(MFNAME)
gf=open.nc(GFNAME)
gd=open.nc(GDNAME)

LATC=var.get.nc(gf,"LAT")
LONC=var.get.nc(gf,"LON")

LATD=var.get.nc(gd,"LATD")
LOND=var.get.nc(gd,"LOND")
HGT=var.get.nc(gf,"HT")

latrg=range(LATD)+c(-1,1)
lonrg=range(LOND)+c(-1,1)

numsites1=length(qxgrids$POINTS)
qxgrids=subset(qxgrids,LONG>lonrg[1] & LONG<lonrg[2] & LAT>latrg[1] & LAT<latrg[2])
numsites2=length(qxgrids$POINTS)

print(paste0("Keep ", numsites2," sites of ",numsites1,"."))

#LATDU=var.get.nc(gd,"LATU")
#LONDU=var.get.nc(gd,"LONU")

close.nc(gf)
SDataO=length(var.get.nc(aq,"NO2_UGM3")[1,1,]) # a safty method
Sdata=SDataO
print("AQ Data Loading...")
NO2_UGM3=var.get.nc(aq,"NO2_UGM3")[,,1:SDataO]
PM10_UGM3=var.get.nc(aq,"PM10")[,,1:SDataO]
PM25_UGM3=var.get.nc(aq,"PM25_TOT")[,,1:SDataO]
SO2_UGM3=var.get.nc(aq,"SO2_UGM3")[,,1:SDataO]
O3_UGM3=var.get.nc(aq,"O3_UGM3")[,,1:SDataO]
CO_MGM3=var.get.nc(aq,"CO_MGM3")[,,1:SDataO]
NH3_UGM3=var.get.nc(aq,"NH3_UGM3")[,,1:SDataO]

TEMP2=var.get.nc(aq,"TEMP2")[,,1:SDataO]
WDIR10=var.get.nc(aq,"WDIR10")[,,1:SDataO]
WSPD10=var.get.nc(aq,"WSPD10")[,,1:SDataO]

PBL=var.get.nc(aq,"PBLH")[,,1:SDataO]
PRSFC=var.get.nc(aq,"PSFC")[,,1:SDataO]
RGRND=var.get.nc(aq,"SR")[,,1:SDataO]
CFRAC=var.get.nc(aq,"CFRAC")[,,1:SDataO]

PCPN=var.get.nc(aq,"RT")[,,1:SDataO]
RH=var.get.nc(aq,"RH")[,,1:SDataO]
VISM=var.get.nc(aq,"VIS")[,,1:SDataO]

TLEN=length(PM25_UGM3[1,1,])-24
close.nc(aq)
close.nc(mt)
close.nc(gd)

print("All data loaded.")
#
#######################################################################################
PQX=T
PGR=T
ProWind=T
CW=0.2
dir.create(OUTDIR,recursive=T)
dir.create(paste0(OUTDIR,CASENAME,"//"))
dir.create(paste0(OUTDIR,CASENAME,"//SNOWINFO"))
##################################WINDLANE#####################
wdir2wdir=function(wdirdf,wspddf){
  sector <- cut(wdirdf, breaks = seq(22.5, 382.5, 45),
                labels = c("NE", "E", "SE", "S", "SW", "W","NW", "N"))
  sector[is.na(sector)]="N"
  sector=as.vector(sector)
  sector[wspddf<CW]="C"
  sechis=as.data.frame(table(sector))
  mainwd=sechis$sector[which(sechis$Freq==max(sechis$Freq))][1]
  cr=length(sector[sector=="C"])/length(sector)
  mr=length(sector[sector==mainwd])/length(sector)
  ms=mean(wspddf[sector==mainwd],na.rm=T)
  as=mean(wspddf,na.rm=T)
  wlr=data.frame(MD=mainwd,MR=mr,MS=ms,CR=cr,AS=as)
  return(wlr)
}

#ignore project
extract_snowinfo=function(pkll,obll,obname,obcode,obshade,breaks=200){
  print(paste0("Processing ",obname,"..."))
#  print(obll)
#  print(pkll)
  pkloc=GetLoc(LATD,LOND,pkll)
  obloc=GetLoc(LATD,LOND,obll)
  if(obloc[1]==0 | obloc[2]==0){
    print("Observation point out of range.")
    return(data.frame())
  }
  deltall=(obll-pkll)/breaks
  obseq=data.frame()
  llseq=data.frame()
  for(xi in 1:breaks){
    locid=GetLoc(LATD,LOND,c(obll[1]-(xi-1)*deltall[1],obll[2]-(xi-1)*deltall[2]))
    obseq=rbind(obseq,data.frame(XI=locid[1],YI=locid[2]))
    llseq=rbind(llseq,data.frame(LAT=LATD[locid[1],locid[2]],LON=LOND[locid[1],locid[2]],HGT=HGT[locid[1],locid[2]]))
  }
  if(debug){
  print(obseq)
  print(data.frame(XI=pkloc[1],YI=pkloc[2]))
  }
  obseq=subset(obseq,XI>0)
  obseq=subset(obseq,YI>0)
  obseq=rbind(obseq,data.frame(XI=pkloc[1],YI=pkloc[2]))
  obseq=unique(obseq)
  llseq=unique(llseq)
  if(debug) print(llseq)
  lldis=data.frame()
  for(xi in 1:length(llseq$HGT)){
    lldis=rbind(lldis,data.frame(DIST=geosphere::distm(c(llseq$LON[xi],llseq$LAT[xi]),c(llseq$LON[1],llseq$LAT[1]))))
  }
  llseq=cbind(llseq,lldis)
  deltahgt=llseq$HGT[length(llseq$HGT)]-llseq$HGT[1]
  deltadis=llseq$DIST[length(llseq$HGT)]-llseq$DIST[1]
  llhgt=data.frame()
  for(xi in 1:length(llseq$HGT)){
    llhgt=rbind(llhgt,data.frame(VHGT=llseq$HGT[1]+deltahgt/deltadis*llseq$DIST[xi]))
  }
  if(debug) print(paste(nrow(obseq),nrow(llseq),nrow(llhgt),nrow(lldis)))
  llseq=cbind(obseq,llseq,llhgt)
  
  snowinfo=data.frame()
  for(xi in 1:length(llseq$HGT)){
    #print(llseq$VHGT[xi])
    pk2d=GetPK2D(llseq$XI[xi],llseq$YI[xi],llseq$VHGT[xi],llseq$HGT[xi],obshade)
    if(debug) print(pk2d)
    pk2d$LON=llseq$LON[xi]
    pk2d$LAT=llseq$LAT[xi]
    pk2d$AGL=llseq$HGT[xi]
    pk2d$EXTID=xi
    pk2d$DIST=llseq$DIST[xi]
    snowinfo=rbind(snowinfo,pk2d)
  }
  #write.csv(file=paste0(OUTDIR,CASENAME,"//SNOWINFO//",obname,".csv"),snowinfo)
  for(pn in names(snowinfo)[2:ncol(snowinfo)]){
  snowpsdf=data.frame(Time=Timeget())
  for(xi in 1:length(llseq$HGT)){
    snowps=subset(snowinfo,EXTID==xi)
    snowpsdf=cbind(snowpsdf,VL=snowps[[pn]])
  }
  names(snowpsdf)=c("Time",paste0("ID",1:length(llseq$HGT)))
  #write.csv(file=paste0(OUTDIR,CASENAME,"//SNOWINFO//",obname,"_",pn,".csv"),snowpsdf)
  }
  
  snowpsdf=data.frame(Time=Timeget())
  snowqvdf=data.frame(Time=Timeget())
  ifShad=F
  for(xi in 1:length(llseq$HGT)){
    snowqc=subset(snowinfo,EXTID==xi)[["QC"]]
    snowqv=subset(snowinfo,EXTID==xi)[["QV"]]   
    if(max(subset(snowinfo,EXTID==xi)[["ifShade"]])==1) ifShad=T
    snowpsdf=cbind(snowpsdf,data.frame(VL=snowqc))
    snowqvdf=cbind(snowqvdf,data.frame(VL=snowqv))
  }
  snowpsdf=snowpsdf[,-1]
  snowqvdf=snowqvdf[,-1]
  snowqcline=rowSums(snowpsdf)
  snowqvline=apply(snowqvdf,1,mean)
  localpm=subset(snowinfo,EXTID==1)[["PM25"]]
  localrs=subset(snowinfo,EXTID==1)[["SR"]]
  peakrs=subset(snowinfo,EXTID==length(llseq$HGT))[["SR"]]
  peakcl=subset(snowinfo,EXTID==length(llseq$HGT))[["CLD"]]
  snowprofile=data.frame(Time=Timeget(),QCS=snowqcline,QVS=snowqvline,localpm,peakrs,peakcl,localrs)
  write.csv(file=paste0(OUTDIR,CASENAME,"//SNOWINFO//",obname,"_PEAKINDEX.csv"),snowprofile)
  
  snowidx=data.frame(Time=Timeget(),ifQC=snowprofile$QCS<=qcthread,ifCLD=peakcl<=pcldthread,ifQV=snowprofile$QVS<=qvthread,ifPM=snowprofile$localpm<pmthread,ifRS=peakrs>rsthread,ifLRS=localrs<=peakrs)
  snowidx$Shaded=ifShad
  snowidx$QXCODE=obcode
  snowidx$ifPK=(snowidx$ifQC & snowidx$ifQV & snowidx$ifPM & snowidx$ifRS & snowidx$ifLRS & snowidx$ifCLD)
  if(debug) print(snowidx)
  if(debug) print(subset(snowidx,ifPK==T))
  return(snowidx)
}

toForeStr=function(spki){
  spki=subset(spki,ifPK==T)
  #print(head(spki))
  Days=as.vector(unique(substr(spki$Time,1,10)))
  if(length(Days)>0){
  spkistrout=""
  for(iday in Days){
    ispki=subset(spki,substr(spki$Time,1,10)==iday)
    hours=substr(ispki$Time,12,16)
    hourstr=paste(hours,collapse=",")
    resstr=paste0(iday,"best view time:",hourstr,",")
    resstr=substr(resstr,1,nchar(resstr)-1)
  }
  spkistrout=paste(spkistrout,resstr,sep=";")

  }else{
    spkistrout="No valid days"
  }
  print(Days)
  return(spkistrout)
}

drawspindex=function(cdsnowpeak,startid=16,endid=184){
  qxs=as.vector(unique(cdsnowpeak$QXCODE))
  spdf=data.frame()
  for(iqx in qxs){
    sptmp=subset(cdsnowpeak,QXCODE==iqx)[startid:endid,]
    #print(sptmp)
    spdf=rbind(spdf,sptmp)
  }
  p=ggplot()+geom_tile(aes(x=Time,y=QXCODE,fill=ifPK),color="white",spdf,show.legend=F)+coord_fixed()+theme_bw()+scale_fill_manual(values=c("TRUE"="blue","FALSE"="gray"))+scale_x_discrete(breaks=sptmp$Time[seq(1,length(sptmp$Time),3)])
  p=p+geom_point(aes(x=Time,y=QXCODE,color=Shaded),data=spdf,show.legend = F,size=1,shape=17)+scale_color_manual(values=c("TRUE"="darkgreen","FALSE"="white"))
  p=p+theme(axis.text.x = element_text(angle=90),axis.title.y=element_blank(),axis.title.x=element_blank())+ylab("")+xlab("")
  ggsave(filename=paste0(OUTDIR,CASENAME,"//SNOWINFO//CDPEAKINDEX.png"),p,width=22,height=4.5)
  print(paste0("Saveing:",OUTDIR,CASENAME,"//SNOWINFO//CDPEAKINDEX.png"))
}

drawsprate=function(cdsnowpeak,startid=16,endid=184){
  qxs=as.vector(unique(cdsnowpeak$QXCODE))
  spdf=data.frame()
  for(iqx in qxs){
    sptmp=subset(cdsnowpeak,QXCODE==iqx)[startid:endid,]
    #print(sptmp)
    spdf=rbind(spdf,sptmp)
  }
  tms=as.vector(unique(spdf$Time))
  sprate=data.frame()
  for(itime in tms){
    sptmp=subset(spdf,Time==itime)
    tpeak=length(sptmp$ifPK)
    vpeak=length(subset(sptmp,ifPK==T)$ifPK)
    rpeak=round(vpeak*100/tpeak,2)
    dpeak=data.frame(Time=itime,PRate=rpeak)
    sprate=rbind(sprate,dpeak)
  }

  p=ggplot()+geom_bar(aes(x=Time,y=PRate),color="white",fill="lightblue",sprate,show.legend=F,stat="identity",position="stack")+theme_bw()+scale_x_discrete(breaks=sprate$Time[seq(1,length(sprate$Time),3)])
  p=p+theme(axis.text.x = element_text(angle=90),axis.title.y=element_blank(),axis.title.x=element_blank())+ylab("Snowmountain Index")+xlab("")+ylim(c(0,100))
  ggsave(filename=paste0(OUTDIR,CASENAME,"//SNOWINFO//CDPEAKRATE.png"),p,width=11,height=4.0)
  print(paste0("Saveing:",OUTDIR,CASENAME,"//SNOWINFO//CDPEAKRATE.png"))
}



pmthread=100
rsthread=350
shadthread=100
pcldthread=30
qvthread=10
qcthread=0
skipprocess=F
cdsnowpeak=data.frame()
debug=F
if(!skipprocess){
for(pt in 1:length(qxgrids$NAMES)){
  spki=extract_snowinfo(c(as.vector(pkgrids[1,]$LAT),as.vector(pkgrids[1,]$LONG)),c(as.vector(qxgrids[pt,]$LAT),as.vector(qxgrids[pt,]$LONG)),obname=as.vector(qxgrids[pt,]$NAMES),obcode=as.vector(qxgrids[pt,]$POINTS),obshade=as.vector(qxgrids[pt,]$SHADED))
  spkiStr=toForeStr(spki)
  cdsnowpeak=rbind(cdsnowpeak,spki)
}
write.csv(file=paste0(OUTDIR,CASENAME,"//SNOWINFO//CDPEAKINDEX.csv"),cdsnowpeak)
}
cdsnowpeak=read.csv(file=paste0(OUTDIR,CASENAME,"//SNOWINFO//CDPEAKINDEX.csv"),header=T)
drawspindex(cdsnowpeak)
drawsprate(cdsnowpeak)
