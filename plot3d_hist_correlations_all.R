path='/Users/taoxianming/Documents/research/gongan/page_methylation'
###
load(paste0(path,'/3Dface/matGPA_luhj2500xxyyzz.RData'))
phe=read.csv(paste0(path,'/phe/TangKun3Dand2500.csv'))[,c('luhj3D','Me.id','age')]
##
trainPhe=na.omit(phe[1001:nrow(phe),c(1,3)])
###train set 1462 samples
matTrainUse=merge(trainPhe,data.frame(rownames(matGPA),matGPA),by.x=1,by.y=1)
####
matUse=matTrainUse
age=matUse[,2]
rownames(matUse)=matUse[,1];matUse=matUse[,-c(1,2)];
colnames(matUse)=1:ncol(matUse)
##969 3D test samples, 988 age, 959
correlations=apply(matUse,2,function(x) cor(age,x))## -1 -- 1
###get result
resPath=paste0(path,'/3Dface')
sample='mean'
##get x y z coordinates
v=read.csv(paste0(resPath,'/',sample,'_xxyyzz.csv'))[,1]
num=length(v)/3
xCol=1:num;yCol=(num+1):(num*2);zCol=(num*2+1):(num*3);
v=data.frame(x=v[xCol],y=v[yCol],z=v[zCol])
##xyz: -80 80;-100 100; -50 50
imp=data.frame(ximp=correlations[xCol],yimp=correlations[yCol],zimp=correlations[zCol])
##
col=1
m1=0.05
m2=m1+0.1
sum(abs(imp[,col])>=m1 & abs(imp[,col])<m2)
# get color by color panel
colNum=20
color=colorRampPalette(c("blue","white","red"))(colNum)
##
name=c('x','y','z')
xyzCor=NULL
for(i in 1:3){
  pdf(paste0(resPath,'/',name[i],'.cor.hist.pdf'))
  a=hist(x=imp[,i],breaks=seq(-1,1,0.1),plot=F)
      heights=a$counts/length(imp[,i])
      barplot(height = heights,col=color,xaxt ="n",space=0,
        ylab = "Proportion",#ylim=c(0,0.2),
        xlab = "Correlation",
        main = "",font.lab=2,font.axis=1,cex.axis=2,cex.lab=1.5)
  axis(1,at=seq(0,20,5),labels=c(-1,-0.5,0,0.5,1),las=1,font=1,cex.axis=2)

  if(F){
  if(i<3){
  a=hist(x=imp[,i],breaks=seq(-1,1,0.1),col=color,ylab='Counts',
  xlab='correlations',main='',
  font.lab=2,font.axis=1,cex.axis=1.4,cex.lab=1.5)
  }else{
    a=hist(x=imp[,i],breaks=seq(-1,1,0.1),col=color,ylab='Counts',
  xlab='correlations',main='',ylim=c(0,1000),
  font.lab=2,font.axis=1,cex.axis=1.4,cex.lab=1.5)
  }
  }


  if(i==1)xyzCor=a$mids
  xyzCor=rbind(xyzCor,a$counts,round(a$counts/length(imp[,i]),2))
  nr=nrow(xyzCor)
  rownames(xyzCor)[(nr-1):nr]=paste0(name[i],c('Count','Frequency'))
  dev.off()
}
colnames(xyzCor)=xyzCor[1,]
xyzCor=xyzCor[-1,]
write.csv(xyzCor,file=paste0(resPath,'/xyz_cor_hist.csv'),quote=F)

