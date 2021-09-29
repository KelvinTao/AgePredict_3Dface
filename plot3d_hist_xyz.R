
###get result
path='/Users/taoxianming/Documents/research/gongan/page_methylation'
resPath=paste0(path,'/3Dface')
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

## get color by color panel
colNum=10
##
colors=list(
reds=colorRampPalette(c("white", "red"))(colNum),
greens=colorRampPalette(c("white", "green"))(colNum),
blues=colorRampPalette(c("white", "blue"))(colNum)
)

##
name=c('x','y','z')
for(i in 1:3){
  pdf(paste0(resPath,'/',name[i],'.imp.hist.pdf'))
   a=hist(x=imp[,i],breaks=seq(0,100,10),plot=F)
      heights=a$counts/length(imp[,i])
      barplot(height = heights,col=colors[[i]],xaxt ="n",space=0,
        ylab = "Proportion",#ylim=c(0,0.2),
        xlab = "Feature importance",
        main = "",font.lab=2,font.axis=1,cex.axis=2,cex.lab=1.5)
      axis(1,at=seq(0,10,2),labels=seq(0,100,20),las=1,font=1,cex.axis=2)
  if(F){
  a=hist(x=imp[,i],breaks=breaks,col=colors[[i]],ylab='Counts',
  xlab='importance',main='',
  font.lab=2,font.axis=1,cex.axis=1.4,cex.lab=1.5)
  }
  #ylim=c(0,1500),
  if(i==1)xyzImp=a$mids
  xyzImp=rbind(xyzImp,a$counts,round(a$counts/length(imp[,i]),2))
  nr=nrow(xyzImp)
  rownames(xyzImp)[(nr-1):nr]=paste0(name[i],c('Count','Frequency'))
  dev.off()
}
colnames(xyzImp)=xyzImp[1,]
xyzImp=xyzImp[-1,]
write.csv(xyzImp,file=paste0(resPath,'/xyz_imp_hist.csv'),quote=F)



