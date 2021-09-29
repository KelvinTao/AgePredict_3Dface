
###get result
path='/Users/taoxianming/Documents/research/gongan/page_methylation'
##get x y z coordinates
v=read.csv(paste0(path,'/3Dface/mean_xxyyzz.csv'))[,1]
num=length(v)/3
xCol=1:num;yCol=(num+1):(num*2);zCol=(num*2+1):(num*3);
v=data.frame(x=v[xCol],y=v[yCol],z=v[zCol])
##coordinates: -80 80;-100 100; -50 50
##get feature importance
load(paste0(path,'/result/3Dpred2500/model/gongan1462luhj.model.imp.RData'))
imp=impRes$importance[,1]#/100*255
imp=data.frame(ximp=imp[xCol],yimp=imp[yCol],zimp=imp[zCol])
#x y z
#apply(imp,2,mean)
#min 0.00000000 0.06188304 0.03932612
#max  66.21611  72.39724 100.00000
#mean 16.07237 15.10039 37.05791
#median 8.043487  8.884362 30.915796
## get color by color panel
colNum=100
library("RColorBrewer")
Rcolor<-function(col,imp){
  index=ceiling(imp);
  if(index==0)index=1;
  color=col[index]##i:1-100;imp:0-100
  return(color)
}
##
reds=colorRampPalette(c("white", "red"))(colNum)
greens=colorRampPalette(c("white", "green"))(colNum)
blues=colorRampPalette(c("white", "blue"))(colNum)
xColor=unlist(lapply(imp$ximp,function(x) Rcolor(reds,x)))
yColor=unlist(lapply(imp$yimp,function(x) Rcolor(greens,x)))
zColor=unlist(lapply(imp$zimp,function(x) Rcolor(blues,x)))
##get all color or by coordinates
colSum=col2rgb(xColor)+col2rgb(yColor)+col2rgb(zColor)
#mmax=max(colSum);mmin=min(colSum);##from circos
mmin=255;mmax=765;
colSum=(colSum-mmin)/(mmax-mmin);
allColor=rgb(t(colSum))
##
#xyzRGB=rgb(imp/100)
##
colors=list(xColor,yColor,zColor,allColor)#,xyzRGB)
names(colors)=c('xColor','yColor','zColor','allColor')#,'xyzRGB')
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
layout3d(matrix(1:4, 1, 4))
#layout3d(t(matrix(1:12, 4, 3)))
 #rgl.close();#rgl.open()
#for (i in seq_along(degrees)){
for (i in 2:2){
  for(ci in 1:4){
  next3d();rgl.viewpoint(zoom=zoomValue);
  par3d(userMatrix = rotationMatrix(degrees[i]*pi/180,
    degrees2[i]*pi/180,1, 0),zoom=0.5)##zoom little, image bigger
  plot3d(v,size=4,col=colors[[ci]],xlim=c(-80,80), ylim=c(-100,100), 
    zlim=c(-50,50),box=F,aspect=F,type='p',axes=F,
    xlab='',ylab='',zlab='')#,lit=F);#ps,axes
  #axes3d(c('x', 'y'))
  #text3d(0,100,0,sides[i],font=1)#,cex=2  000
 }
}


jpgFile=paste0(path,'/3Dface/xyz.importance.center.',zoomValue,'.jpg')
rgl.snapshot(jpgFile)




