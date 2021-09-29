path='/Users/taoxianming/Documents/research/gongan/page_methylation'
###
load(paste0(path,'/3Dface/matGPA_luhj2500xxyyzz.RData'))
phe=read.csv(paste0(path,'/phe/TangKun3Dand2500.csv'))[,c('luhj3D','Me.id','age')]
##
#testPhe=phe[1:1000,c(1,3)]##
trainPhe=na.omit(phe[1001:nrow(phe),c(1,3)])
###train set 1463 samples
matTrainUse=merge(trainPhe,data.frame(rownames(matGPA),matGPA),by.x=1,by.y=1)
#matTestUse=na.omit(merge(testPhe,data.frame(rownames(matGPA),matGPA),by.x=1,by.y=1))
#matUse=rbind(matTrainUse,matTestUse)
####
matUse=matTrainUse
age=matUse[,2]
rownames(matUse)=matUse[,1];matUse=matUse[,-c(1,2)];
colnames(matUse)=1:ncol(matUse)
##969 3D test samples, 988 age, 959
correlations=apply(matUse,2,function(x) cor(age,x))## -1 -- 1

###get result
resPath='/Users/taoxianming/Documents/research/gongan/page_methylation/3Dface'
samples=c('L01171','mean')
sample=samples[2]
##get x y z coordinates
v=read.csv(paste0(resPath,'/',sample,'_xxyyzz.csv'))[,1]
num=length(v)/3
xCol=1:num;yCol=(num+1):(num*2);zCol=(num*2+1):(num*3);
v=data.frame(x=v[xCol],y=v[yCol],z=v[zCol])
##-80 80;-100 100; -50 50
##get feature importance
imp=correlations
imp=data.frame(ximp=imp[xCol],yimp=imp[yCol],zimp=imp[zCol])
#x y z
#min 0.00000000 0.07108802 0.02635692 
#max  65.61339  70.79111 100.00000
#mean 15.95408 14.86488 37.02508
#median 7.936389  8.811817 30.941206
## get color by color panel
colNum=200
library("RColorBrewer")
Rcolor<-function(col,imp){
  index=ceiling(imp*100+100);##0-200
  if(index==0)index=1;
  color=col[index]##i:1-100;imp:0-100
  return(color)
}
##
reds=colorRampPalette(c("blue", "white","red"))(colNum)
greens=colorRampPalette(c("blue", "white","red"))(colNum)
blues=colorRampPalette(c("blue", "white","red"))(colNum)
xColor=unlist(lapply(imp$ximp,function(x) Rcolor(reds,x)))
yColor=unlist(lapply(imp$yimp,function(x) Rcolor(greens,x)))
zColor=unlist(lapply(imp$zimp,function(x) Rcolor(blues,x)))
##get all color or by coordinates
colSum=col2rgb(xColor)+col2rgb(yColor)+col2rgb(zColor)
mmax=max(colSum);mmin=min(colSum);##
#mmin=255;mmax=765;
colSum=(colSum-mmin)/(mmax-mmin);
allColor=rgb(t(colSum))

##
colors=list(xColor,yColor,zColor,allColor)
names(colors)=c('xColor','yColor','zColor','allColor')
##plot 3D
library(rgl)
degrees=c(45,-0.01,-45)
degrees2=c(1,-0.01,-1)
sides=c('left','center','right');
zoomValue=0.5;
###
#for(ci in 1:length(colors)){
#ci=1  ##run seperately successful
par3d(windowRect = c(0,0,1920, 1080))
#layout3d(matrix(1:3, 1, 3))
layout3d(matrix(1:4,1,4))
 #rgl.close();#rgl.open()
for (i in 2:2){
  for(ci in 1:4){
  next3d();rgl.viewpoint(zoom=zoomValue);
  par3d(userMatrix = rotationMatrix(degrees[i]*pi/180,
    degrees2[i]*pi/180,1, 0),zoom=0.6)
  plot3d(v,size=4,col=colors[[ci]],xlim=c(-80,80), ylim=c(-100,100), 
    zlim=c(-50,50),box=F,aspect=F,type='p',axes=F,
    xlab='',ylab='',zlab='')#,lit=F);#ps,axes
  #axes3d(c('x', 'y'))
  #text3d(0,100,0,sides[i],font=1)#,cex=2  000
 }
}


jpgFile=paste0(resPath,'/xyz.cor.center.',zoomValue,'.jpg')
rgl.snapshot(jpgFile)
#}



