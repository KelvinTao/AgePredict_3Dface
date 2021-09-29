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
sample='mean'
##get x y z coordinates
v=read.csv(paste0(resPath,'/',sample,'_xxyyzz.csv'))[,1]
num=length(v)/3
xCol=1:num;yCol=(num+1):(num*2);zCol=(num*2+1):(num*3);
v=data.frame(x=v[xCol],y=v[yCol],z=v[zCol])
##-80 80;-100 100; -50 50
##get feature importance
imp=correlations
imp=data.frame(ximp=imp[xCol],yimp=imp[yCol],zimp=imp[zCol])
##
#col=1
#m1=0
#m2=m1+0.1
#sum(abs(imp[,col])>=m1 & abs(imp[,col])<m2)

##sum(abs(imp[,1])<0.1)
#apply(imp,2,min)
#apply(imp,2,max)
#apply(imp,2,mean)
#apply(imp,2,median)
#x y z
#min -0.5197644 -0.4495358 -0.6108473
#max 0.5021587 0.5397384 0.6323678 
#mean -0.002008602  0.018142269  0.034390545
#median -0.0004056428 -0.0115948623 -0.0167715484
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
blues=greens=reds=colorRampPalette(c("blue","white","red"))(colNum)
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
par3d(windowRect = c(0,0,1920, 1080))
layout3d(t(matrix(1:12, 3, 4)))
for(ci in 1:4){
for (i in seq_along(degrees)){
  next3d();rgl.viewpoint(zoom=zoomValue);
  par3d(userMatrix = rotationMatrix(degrees[i]*pi/180,
    degrees2[i]*pi/180,1, 0),zoom=0.6)
  plot3d(v,size=4,col=colors[[ci]],xlim=c(-80,80), ylim=c(-100,100), 
    zlim=c(-50,50),box=F,aspect=F,type='p',axes=F,
    xlab='',ylab='',zlab='')#,lit=F);#ps,axes
 }
}

jpgFile=paste0(resPath,'/xyz.cor.all.',zoomValue,'.jpg')
rgl.snapshot(jpgFile)



