#SNOWSHADE CALC
#20220606
library(inlmisc)
library(raster)
library(geosphere)
library(sp)
library(doParallel)
library(RNetCDF)
library(ggplot2)

SNOWWDR="G:\\OneDrive\\工作文档\\大气工作\\20210920-雪山指数预报\\"
mapcd=fortify(rgdal::readOGR(dsn="G:\\Projects\\成都边界区县",layer="cd"))
mapbd=fortify(rgdal::readOGR(dsn="G:\\OneDrive\\工作文档\\大气工作\\20190227-CAMx网格制作\\CDBRD",layer="Chengdu_frame"))

gd=open.nc(paste0(SNOWWDR,"GRIDDOT2D"))
gr=open.nc(paste0(SNOWWDR,"GRIDCRO2D"))
LATD=var.get.nc(gd,"LATD")
LOND=var.get.nc(gd,"LOND")
HGT=var.get.nc(gr,"HT")
LAT=var.get.nc(gr,"LAT")
LON=var.get.nc(gr,"LON")
RECALC=F

PFMT="wmf"
PDPI=96

#coordinates for snow peak
latp=31.107
lonp=102.905

windowsFonts(HT=windowsFont("黑体"),
             RMN=windowsFont("Times New Roman"),
             ST=windowsFont("宋体"))
plotthemev=theme(
  axis.text.x = element_text(family="RMN",size=12,angle=90),
  axis.text.y = element_text(family="RMN",size=12,angle=0),
  axis.title=element_text(family="RMN",size=13, colour="black",face="bold"),
  strip.text =element_text(family="RMN",size=12, colour="black"),
  strip.background = element_rect(fill = "gray90"),
  legend.title=element_blank(), #element_text(family="HT", size=size_t,angle=0,face="plain"),
  legend.text=element_text(family="RMN",size=10, colour="black",face="plain"),
  legend.position="right", 
  legend.direction="vertical",
  legend.key.width=unit(0.3,"cm"),
  legend.key.height=unit(1.0,"cm")
)

plotthemeh=theme(
  axis.text.x = element_text(family="RMN",size=12,angle=90),
  axis.text.y = element_text(family="RMN",size=12,angle=0),
  axis.title=element_text(family="RMN",size=13, colour="black",face="bold"),
  strip.text =element_text(family="RMN",size=12, colour="black"),
  strip.background = element_rect(fill = "gray90"),
  legend.title=element_blank(), #element_text(family="HT", size=size_t,angle=0,face="plain"),
  legend.text=element_text(family="RMN",size=10, colour="black",face="plain"),
  legend.position="top", 
  legend.direction="horizontal",
  legend.key.width=unit(0.3,"cm"),
  legend.key.height=unit(1.0,"cm")
)

dem=raster(paste0(SNOWWDR,"sc_srtm_ne250.tif"))
names(dem)="value"
demcd=crop(dem,rgdal::readOGR(dsn="G:\\OneDrive\\工作文档\\大气工作\\20190227-CAMx网格制作\\CDBRD",layer="Chengdu_frame"))
nras=raster(ext=extent(demcd),crs=crs(demcd),res=res(demcd)*6)
demr=resample(demcd,nras,method="bilinear")

dempt=as.data.frame(as(demr, "SpatialPixelsDataFrame"))

ext_dem=function(latv,lonv,minhgt=5000,showplot=F){
  coords=rbind(c(latv,lonv),c(latp,lonp))
  lines=SpatialLines(list(Lines(list(Line(list(x=c(lonv,lonp),y=c(latv,latp)))),ID=1)))
  crs(lines)=crs(dem)
  heights=as.data.frame(extract(dem,lines,cellnumbers=T,along=T)[[1]])
  coord_xy=coordinates(dem)[heights[,1],]
  dists_xy=distGeo(p1=c(lonv,latv),p2=coord_xy)
  delh=heights$value[length(heights$value)]-heights$value[1]
  deld=dists_xy[length(dists_xy)]-dists_xy[1]
  calh=delh/deld*dists_xy+heights$value[1]
  height_xy=data.frame(dist=dists_xy,rhgt=heights$value,chgt=calh,dhgt=calh-heights$value)
  plot(height_xy$dist,height_xy$rhgt,type="l",ylab="Height",xlab="Distance")
  lines(height_xy$dist,height_xy$chgt,col="green")
  Shaded=F
  shaddf=subset(height_xy,dhgt<0)
  if(nrow(shaddf)>0){
    if(shaddf[1,]$rhgt<minhgt) Shaded=T
    points(shaddf$dist,shaddf$rhgt,col="red")
  }
  if(showplot){
    p=ggplot()
    p=p+geom_path(aes(x=long,y=lat,group=group),mapcd,color="gray50")
    p=p+geom_point(aes(x=lonp,y=latp,shape="OBS"),color="black",size=2)
    p=p+geom_text(aes(x=lonp*1.0015,y=latp*1.002,label="Yaomei Peak"),color="black",size=3.0,family="RMN")
    p=p+geom_point(aes(x=lonv,y=latv,shape="PEK"),color="black",size=2)
    p=p+geom_segment(aes(x=lonv,y=latv,xend=lonp,yend=latp),lineend = "butt",color="red")
    p=p#+coord_map()
    p=p+scale_shape_manual(name=NA,
                           breaks=c('OBS', 'PEK'),
                           values=c('OBS'=17, 'PEK'=16),
                           labels=c("Observation point","Peak point"),
                           guide = guide_legend(override.aes = list(
                             linetype = c("blank", "blank"),
                             shape = c(16, 17),
                             size=c(2,2))))
    p=p+theme_bw()+plotthemeh
    p=p+labs(x="Longitude(°)",y="Latitude(°)")
    p
    
    q=ggplot()
    q=q+geom_point(aes(x=height_xy$dist[1],y=height_xy$rhgt[1],shape="PEK"),color="black",size=2)
    q=q+geom_point(aes(x=height_xy$dist[nrow(height_xy)],y=height_xy$rhgt[nrow(height_xy)],shape="OBS"),color="black",size=2)
    q=q+geom_line(aes(x=dist,y=rhgt,group=1),height_xy)
    q=q+geom_line(aes(x=dist,y=chgt,group=1),height_xy,color="green",linetype="dashed")
    #q=q+geom_segment(aes(x=height_xy$dist[1],y=height_xy$rhgt[1],xend=height_xy$dist[nrow(height_xy)],yend=height_xy$rhgt[1]),lineend = "butt",color="red")
    q=q+scale_shape_manual(name=NA,
                           breaks=c('OBS', 'PEK'),
                           values=c('OBS'=17, 'PEK'=16),
                           labels=c("Observation point","Peak point"),
                           guide = guide_legend(override.aes = list(
                             linetype = c("blank", "blank"),
                             shape = c(16, 17),
                             size=c(2,2))))
    if(Shaded){
      q=q+geom_point(aes(x=dist,y=rhgt),shaddf,color="red",shape=1)
      q=q+geom_hline(yintercept = shaddf[1,]$rhgt,color="red",linetype="dashed")
    }
    q=q+theme_bw()+plotthemeh
    q=q+labs(x="Distance(m)",y="Elevation(m)")
    q
    
    pl=ggpubr::ggarrange(p,
                         q,
                         nrow=1,
                         common.legend = T,
                         align = "hv",
                         legend = "top")
    return(pl)
  }else{
    return(data.frame(LAT=latv,LON=lonv,Shaded=Shaded))
  }
}

gen_shad_map=function(x,y,n=5,showplot=F){
  latrang=range(LATD[x,y],LATD[x+1,y],LATD[x,y+1],LATD[x+1,y+1])
  lonrang=range(LOND[x,y],LOND[x+1,y],LOND[x,y+1],LOND[x+1,y+1])
  latcent=LAT[x,y]
  loncent=LON[x,y]
  dlat=(latrang[2]-latrang[1])/n/2
  dlon=(lonrang[2]-lonrang[1])/n/2
  
  latseq=seq(latrang[1]+dlat,latrang[2]-dlat,length.out=n-1)
  lonseq=seq(lonrang[1]+dlon,lonrang[2]-dlon,length.out=n-1)
  newll=expand.grid(LAT=latseq,LON=lonseq)
  
  plot(y=c(LATD[x,y],LATD[x+1,y],LATD[x,y+1],LATD[x+1,y+1]),
       x=c(LOND[x,y],LOND[x+1,y],LOND[x,y+1],LOND[x+1,y+1]),
       xlab="Longitude",
       ylab="Latitude")
  points(x=newll$LON,
         y=newll$LAT,
         col="red")
  points(x=loncent,
         y=latcent,
         col="blue")
  
  if(showplot){
    gridpoly=data.frame(y=c(LATD[x,y],LATD[x+1,y],LATD[x+1,y+1],LATD[x,y+1],LATD[x,y]),
                     x=c(LOND[x,y],LOND[x+1,y],LOND[x+1,y+1],LOND[x,y+1],LOND[x,y]))
    
    p=ggplot()
    p=p+geom_path(aes(x=x,y=y),gridpoly,size=1)
    p=p+geom_point(aes(x=LON,y=LAT),newll,color="red",shape=1,size=2)
    p=p+coord_map()
    p=p+theme_bw()+plotthemev
    p=p+labs(x="Longitude(°)",y="Latitude(°)")
    p
    return(p)
  }else{
    return(newll)
  }

}

check_grid=function(x,y){
  grid_ll=gen_shad_map(x,y)
  check_by_id=function(xi){
    grdf=grid_ll[xi,]
    return(ext_dem(grdf$LAT,grdf$LON))
  }
  
  cl<- makeCluster(10)      
  registerDoParallel(cl)
  
  chekdf=foreach(
    i=1:nrow(grid_ll), 
    .combine=rbind, 
    .packages = c("raster", "geosphere"),
    .export=c("grid_ll","check_by_id","ext_dem","latp","lonp","dem")
  ) %dopar% check_by_id(i)
  stopCluster(cl)
  
  shadrate=length(chekdf$Shaded[chekdf$Shaded==T])/length(chekdf$Shaded)
  return(list(dist=chekdf,rate=shadrate))
}

if(RECALC){
  distdf=data.frame()
  shadre=data.frame()
  pb=txtProgressBar(max=length(HGT[,1]),style=3)
  for(xi in 1:length(HGT[,1])){
    for(yi in 1:length(HGT[1,])){
      checkres=check_grid(xi,yi)
      distdf=rbind(distdf,checkres$dist)
      shadre=rbind(shadre,
                   data.frame(X=xi,Y=yi,LAT=LAT[xi,yi],LON=LON[xi,yi],ShadeRate=checkres$rate))
    }
    setTxtProgressBar(pb,xi)
  }
  
  close(pb)
  write.csv(file=paste0(SNOWWDR,"SNOWSHADEF.csv"),shadre)
  write.csv(file=paste0(SNOWWDR,"SNOWDISTRF.csv"),distdf)
  pip=point.in.polygon(point.x=shadre$LON,point.y=shadre$LAT,pol.x=mapbd$long,pol.y=mapbd$lat)
  shadre2=shadre[pip==1,]
  pip=point.in.polygon(point.x=distdf$LON,point.y=distdf$LAT,pol.x=mapbd$long,pol.y=mapbd$lat)
  distdf2=distdf[pip==1,]
  write.csv(file=paste0(SNOWWDR,"SNOWSHADE.csv"),shadre2)
  write.csv(file=paste0(SNOWWDR,"SNOWDISTR.csv"),distdf2)
}else{
  shadre=read.csv(file=paste0(SNOWWDR,"SNOWSHADEF_YMF.csv"))
  pip=point.in.polygon(point.x=shadre$LON,point.y=shadre$LAT,pol.x=mapbd$long,pol.y=mapbd$lat)
  shadre2=shadre[pip==1,]
  write.csv(file=paste0(SNOWWDR,"SNOWSHADE.csv"),shadre2)
}

plot_shade=function(){
  p=ggplot()
  p=p+geom_point(aes(x=LON,y=LAT,size=ShadeRate*100,color=ShadeRate*100),shadre2)+scale_size_area(max_size=2.5)
  p=p+scale_color_gradientn("地形遮挡概率",colours = hcl.colors(30, palette = "Blue-Red"))
  p=p+geom_path(aes(x=long,y=lat,group=group),mapcd,color="gray50")
  p=p+geom_point(aes(x=lonp,y=latp),color="black",shape=17,size=2)
  p=p+theme_bw()+plotthemev+coord_map()+guides(size = "none") 
  p=p+labs(x="Longitude(°)",y="Latitude(°)")
  return(p)
}

plot_dem=function(icrop=T){
  if(icrop){
    pip=point.in.polygon(point.x=dempt$x,point.y=dempt$y,pol.x=mapbd$long,pol.y=mapbd$lat)
    dempt=dempt[pip==1,]
  }
  p=ggplot()
  p=p+geom_tile(aes(x=x,y=y,fill=value),dempt)
  p=p+scale_fill_gradientn("地形遮挡概率",colours = terrain.colors(50))
  p=p+geom_path(aes(x=long,y=lat,group=group),mapcd,color="gray50")
  p=p+geom_point(aes(x=lonp,y=latp),color="black",shape=17,size=2)
  p=p+theme_bw()+plotthemev+coord_map()+guides(size = "none") 
  p=p+labs(x="Longitude(°)",y="Latitude(°)")
  return(p)
}

terr_occ=ext_dem(30.7,103.6,showplot=T)
ggsave(filename = paste0(SNOWWDR,"论文绘图\\TERR_OCC.",PFMT),dpi=PDPI,width=8,height=4.0,terr_occ)
samp_grd=gen_shad_map(44,44,showplot=T)
ggsave(filename = paste0(SNOWWDR,"论文绘图\\SAMP_GRD.",PFMT),dpi=PDPI,width=8,height=7,samp_grd)
terr_ocd=plot_shade()
ggsave(filename = paste0(SNOWWDR,"论文绘图\\TERR_GRD.",PFMT),dpi=PDPI,width=8,height=6,terr_ocd)
cdter=plot_dem(F)
ggsave(filename = paste0(SNOWWDR,"论文绘图\\CDDEM.",PFMT),dpi=PDPI,width=8,height=6,cdter)

ter_mgr=ggpubr::ggarrange(cdter,
                     terr_ocd,
                     nrow=1,
                     common.legend = F,
                     align = "hv",
                     legend = "right")

ggsave(filename = paste0(SNOWWDR,"论文绘图\\CDDEMMG.",PFMT),dpi=PDPI,width=8,height=3.6,ter_mgr)
