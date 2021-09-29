
##color bars
resPath='/Users/taoxianming/Documents/research/gongan/page_methylation/3Dface'
##produce color bars
num=200
blues=greens=reds=colorRampPalette(c("blue","white","red"))
colors=list(reds,greens,blues)
varNames=c('X','Y','Z')
minXYZ=c(-0.52,-0.45,-0.61)
maxXYZ=c(0.50,0.54,0.63)
#plot and save
for(i in 1:3){
  pdf(file=paste0(resPath,'/',varNames[i],'.bar.pdf'))
  #par(pin=c(1,5),cex.axis=1,col.axis=cols[i])#
  par(pin=c(1,5),cex.axis=1)
  x=rep(1, num)
  barplot(x,col= colors[[i]](num),horiz=T,border=NA,space=0,xlim=c(0,1),xaxt ="n",yaxt = "n")
  axis(2,at=seq(0,200,50),labels=(seq(0,200,50)-100)/100,las=1,font=1,cex.axis=2)
  if(i==1){#axis(2,at=maxXYZ[i]*100+100,labels=paste0('max: ',maxXYZ[i],''),las=1,font=2,cex.axis=2)
    axis(2,at=maxXYZ[i]*100+100,labels=paste0('max: 0.50        '),las=1,font=2,cex.axis=2)
    axis(2,at=minXYZ[i]*100+100,labels=paste0('min: ',minXYZ[i],'        '),las=1,font=2,cex.axis=2)}
  if(i==2){axis(2,at=maxXYZ[i]*100+100,labels=paste0('max: ',maxXYZ[i],'       '),las=1,font=2,cex.axis=2)
    axis(2,at=minXYZ[i]*100+100,labels=paste0('min: ',minXYZ[i],'        '),las=1,font=2,cex.axis=2)}
  if(i==3){axis(2,at=maxXYZ[i]*100+100,labels=paste0('max: ',maxXYZ[i],''),las=1,font=2,cex.axis=2)
    #axis(2,at=minXYZ[i]*100+100,labels=paste0('min: ',minXYZ[i],''),las=1,font=2,cex.axis=2)
    axis(2,at=minXYZ[i]*100+100,labels=paste0('min: -0.61',''),las=1,font=2,cex.axis=2)}

  mtext(paste0('correlations of ',varNames[i]),side=3,font=2,line=-3.5,adj=0.48,outer=T,cex=2)
  dev.off()
}
