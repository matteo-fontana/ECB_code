#rolling functional boxplot analysis on Santa Barbara Vulcanelli

#Load packages

library(fda)
library(roahd)
library(ggplot2)
library(raster)
library(svglite)
library(cluster)
library(fpc)
library(data.table)
library(lubridate)
library(dendextend)
library(progress)



####Load data and preprocessing######

load('caltanissetta_envisat.rdata')
head(dt)

#compute geodesic distance

distance_grid=c(750)
source('custom_fbplot_x.r', echo=TRUE)
library(ggmap)
environment(fbplot.custom_x) <- asNamespace('roahd')


for(distance in distance_grid){

  distance=750
#distances are in meters

#define vulcanelli area when radius<distance
dt$is_santabarbara=dt$d_vulcanelli<distance


#Plot data
dt_hgt = dt[dt$is_santabarbara,7:38]*10
dt=dt[dt$is_santabarbara,]
time=as.Date(colnames(dt_hgt),format = "%d/%m/%Y")

filename=paste('paper_works/recentered/Functional_Boxplot_Exploration/Final/small_area/Data_unsmoothed_',distance,'.pdf',sep='')
pdf(file=filename,onefile = T)
par(mar=c(6,4,4,4)+.1)
matplot(time,t(dt_hgt),type='l',xaxt="n",xlab='',ylab='Vertical Displacement (mm)',main='Displacement Data - Vulcanelli - D<1000m')
axis(1,at=time,labels=time,las=2)
dev.off()

time=as.Date(colnames(dt_hgt),format = "%d/%m/%Y")


#select smoothing

#Bsplines? do we want to check for periodicity

basis=create.bspline.basis(breaks=time,norder = 6)
plot(basis,main='B-Spline Basis')


#Code for selection of smoothing parameter

#cairo_pdf(filename = 'new_EDA/bsplinebasis.pdf',height = 8.27 ,width = 11.69)
#plot(basis,main='B-Spline Basis')
#dev.off()

#svglite('new_EDA/bsplinebasis.svg',height = 8.27,width = 11.69)
#plot(basis,main='B-Spline Basis')
#dev.off()


#set up lambda vectors
# lambda=seq(0,15)
# lambda=10^lambda
# 
# gcv=numeric(length(lambda))
# 
# for (i in 1:length(lambda)) {
#   functionalPar = fdPar(fdobj=basis, Lfdobj=3, lambda=lambda[i])  
#   gcv[i] = sum(smooth.basis(time, t(as.matrix(dt_hgt)), functionalPar)$gcv)
# }
# 
# plot(log10(lambda),gcv, type='l')
# lambda[which.min(gcv)]
# 
# lambda_dt=data.table(loglambda=log10(lambda),gcv=gcv)
# 
# fd_obj = fdPar(fdobj=basis, Lfdobj=3, lambda=10^10)  
# 
# disp_fd=smooth.basis(time, t(as.matrix(dt_hgt)), fd_obj)
# plot(disp_fd)

###new map
#define bounding box


#compute zoom
lonrange=extendrange(dt$Lon,f=.05)
latrange=extendrange(dt$Lat,f=.05)

center_lon=((range(dt$Lon)[2]-range(dt$Lon)[1])/2) +  range(dt$Lon)[1]
center_lat=((range(dt$Lat)[2]-range(dt$Lat)[1])/2) +  range(dt$Lat)[1]

#download background map

map=get_stamenmap(bbox = c(left=lonrange[1],right=lonrange[2],top=latrange[2],bottom=latrange[1]),maptype = 'terrain',zoom=15)


filename=paste('paper_works/recentered/Functional_Boxplot_Exploration/Final/PAPER_GRADE_WORK/Data_map_',distance,'.pdf',sep='')
pdf(file=filename,onefile = T)
map_area=ggmap(map) + 
  geom_point(data=dt,aes(x=Lon,y=Lat),alpha=.7,size=1)+
  labs(title='Area of Interest - ENVISAT',y='Latitude (WGS84)',x='Longitude (WGS84)')
print(map_area)
dev.off()

#now implement the rolling system.
#start from the first six observations: from 2004 onwards
starting_obs=6

windowsize=12
recenter=T


#windowsize is in months
#let's compute 2 cuts: 1 for incremental smoothing, the second one for actual processing:

#plot data 
time_str=format(time,'%d-%m-%Y')
list_rollingdata_cut_s=vector(mode='list', length=length(time)-starting_obs+1)
list_time_s=vector(mode='list', length=length(time)-starting_obs+1)
list_time=list_time_s
#recentered displacement data.

for (index in 0:(length(time)-starting_obs)) {
  
  mask=time>(time[starting_obs+index]-months(windowsize))
  start=match(T,mask)
  vel_eval_temp=as.matrix(dt_hgt[,1:(starting_obs+index)])
  list_rollingdata_cut_s[[index+1]]=vel_eval_temp
  list_time_s[[index+1]]=time[1:(starting_obs+index)]
  list_time[[index+1]]=time[start:(starting_obs+index)]
}

#smoothing, eval and recenter:

list_rollingdata_disp=vector(mode='list', length=length(time)-starting_obs+1)
list_rollingdata_smooth=vector(mode='list', length=length(time)-starting_obs+1)

for(i in 1:length(list_rollingdata_disp)){
  
  basis=create.bspline.basis(breaks=list_time_s[[i]],norder = 6)
  functionalPar = fdPar(fdobj=basis, Lfdobj=3, lambda=10^10) 
  disp_fd=smooth.basis(list_time_s[[i]], t(as.matrix(list_rollingdata_cut_s[[i]])), functionalPar)
  list_rollingdata_smooth[[i]]=disp_fd
  temp=t(eval.fd(list_time[[i]],disp_fd$fd))
  list_rollingdata_disp[[i]]=temp-temp[,1]
  
}



#Code for publication plot
# dt_bplot=dt_hgt
# dt_bplot[,1:5]=0
# 
# 
# 
# filename=paste('paper_works/recentered/Functional_Boxplot_Exploration/Final/small_area/scalar/sequential_bplot',distance,'.pdf',sep='')
# pdf(file=filename,onefile = T)
# boxplot(dt_bplot,xrange=range(time), main='Time-Indexed Boxplots',xlab='')
# dev.off()


filename=paste('paper_works/recentered/Functional_Boxplot_Exploration/Final/PAPER_GRADE_WORK/unsmoothed_rolling_data_',distance,'.pdf',sep='')
pdf(file=filename,onefile = T)

for(i in 1:length(list_rollingdata_cut)){
  
  title=paste('Data @',time_str[starting_obs+i-1],sep = ' ')
  matplot(list_time[[i]],t(list_rollingdata_cut[[i]]),main=title,type='l',xlim = range(time),xlab,ylim = c(-20,25))
}
dev.off()


# for (index in 0:(length(time)-starting_obs)) {
#   
#   mask=time>(time[(starting_obs+index)])
#   start=match(T,mask)
#   vel_eval_temp=t(disp_eval)[,start:(starting_obs+index)]
#   if(recenter) {list_rollingdata_disp[[index+1]]=vel_eval_temp - vel_eval_temp[,1]
#   } else {list_rollingdata_disp[[index+1]]=vel_eval_temp}
#   list_time[[index+1]]=time[start:(starting_obs+index)]
# }
# 
# pdf(file='paper_works/recentered/Functional_Boxplot_Exploration/Final/recentered/disp_data_point.pdf',onefile = T)
# matplot(time,disp_eval,main='Displacement Data',type='l',ylim = c(-75,110),xlim = range(time))
# for(i in 1:length(list_rollingdata_disp)){
#   
#   title=paste('Data @',time_str[starting_obs+i-1],sep = ' ')
#   scatter(list_time[[i]],t(list_rollingdata_disp[[i]]),main=title,ylim = c(-75,110),xlim = range(time),xlab)
# }
# dev.off()


#now, smoothing:


filename=paste('paper_works/recentered/Functional_Boxplot_Exploration/Final/PAPER_GRADE_WORK/smoothed_rolling_data_',distance,'.pdf',sep='')
pdf(file=filename,onefile = T)

for(i in 1:length(list_rollingdata_cut_s)){
  
  title=paste('Data @',time_str[starting_obs+i-1],sep = ' ')
  matplot(list_time[[i]],t(list_rollingdata_disp[[i]]),main=title,type='l',xlim = range(time),xlab,ylim = c(-20,25))
}
dev.off()


#let's cycle on a grid of FValues, then:
fvalue_grid=c(3.25)

for (fvalue in fvalue_grid) {
  

filepath_time=paste('paper_works/recentered/Functional_Boxplot_Exploration/Final/PAPER_GRADE_WORK/disp_fboxplot_12m_',fvalue,'_Fval_',distance,'.pdf',sep='')

cairo_pdf(filename=filepath_time,onefile = T)

list_boxplot=vector(mode='list',length = length(list_rollingdata_disp))

for(i in 1:length(list_rollingdata_disp)){
  
  title=paste('Functional Boxplot @',time_str[starting_obs+i-1],sep = ' ')
  test=list_rollingdata_disp[[i]]
  
  test_fdata=fData(as.numeric(list_time[[i]]),test)
  fbox=fbplot.custom_x(test_fdata,Fvalue=fvalue ,main = title,ygrid = list_time[[i]],x_range = range(time),ylim=c(-20,25))
  list_boxplot[[i]]=fbox
}
dev.off()



filepath_map=paste('paper_works/recentered/Functional_Boxplot_Exploration/Final/PAPER_GRADE_WORK/disp_out_map_',fvalue,'_Fval_',distance,'.pdf',sep='')

cairo_pdf(filename=filepath_map,onefile = T)
pb=progress_bar$new(total = length(list_rollingdata_disp))
pb$tick(0)
for(i in 1:length(list_rollingdata_disp)){
  
  title=paste('Outlier Geographical Location @', time_str[starting_obs+i-1],sep = ' ')
  test=list_rollingdata_disp[[i]]
  median_candidates=test[list_boxplot[[i]]$Depth==max(list_boxplot[[i]]$Depth)]
  med=colMeans(matrix(median_candidates,ncol=dim(test)[2]))
  check_down=t(t(test)<med)
  check_down_vec=rowSums(check_down)>=ncol(test)/2
  check_down_vec=check_down_vec+1
  
  is_outlier=rep(0,nrow(dt_hgt))
  col_plotting=rep("black",nrow(dt_hgt))
  col_plotting[list_boxplot[[i]]$ID_outliers]=list_boxplot[[i]]$Col_outlying
  is_outlier[list_boxplot[[i]]$ID_outliers]=check_down_vec[list_boxplot[[i]]$ID_outliers]
  is_outlier=factor(is_outlier,levels = c(0,1,2))
  
  plot=ggmap(map) +
    geom_point(data=dt,aes(x=Lon,y=Lat,col=is_outlier,shape=is_outlier,alpha=is_outlier,size=is_outlier),col=col_plotting)+
    scale_alpha_manual(values=c(0.45,0.8,0.8),drop=F)+
    scale_size_manual(values=c(4,8,8),drop=F)+
    scale_shape_manual(values=c("\u25CF","\u25B2","\u25BC"),drop=F)+
    labs(title=title ,y='Latitude (WGS84)',x='Longitude (WGS84)')
  print(plot)
  
  pb$tick()
}
dev.off()

}
}


#outlier detection for CC?


#let's confront everything with univariate boxplot

filepath='paper_works/recentered/Functional_Boxplot_Exploration/Final/scalar/scalar_boxplotR.pdf'

cairo_pdf(filename=filepath,onefile = T)


for(i in 1:length(list_rollingdata_disp)){

temp_data=t(disp_eval)[,starting_obs+i-1]
bp=boxplot(temp_data, range=3 )

is_outlier=rep(0,33)
bp_out_pos=bp$out[bp$out>=0]
bp_out_neg=bp$out[bp$out<0]
is_outlier[which(temp_data %in% bp_out_pos)]=1
is_outlier[which(temp_data %in% bp_out_neg)]=-1

is_outlier=factor(is_outlier,levels = c(-1,0,1))


title=paste('Zoom on Vulcanelli - Scalar @', time_str[starting_obs+i-1],sep = ' ')
plot=ggmap(map) +
  geom_point(data=dt_proc,aes(x=Lon,y=Lat,col=is_outlier,shape=is_outlier,size=is_outlier,alpha=is_outlier))+
  labs(title=title ,y='Latitude (WGS84)',x='Longitude (WGS84)')+
  scale_color_manual(values=c('navyblue',"black",'red'),drop=F) + 
  scale_size_manual(values=c(8,4,8),drop=F) + 
  scale_shape_manual(values=c("\u25BC","\u25CF","\u25B2"),drop=F)+
  scale_alpha_manual(values=c(0.95,.45,0.95),drop=F)
print(plot)
print(i)

}

dev.off()

