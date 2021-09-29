
##color bars
resPath='/Users/taoxianming/Documents/research/gongan/page_methylation/3Dface'
##produce color bars
numN=100
colors=list(reds=colorRampPalette(c("white", "red")),
greens=colorRampPalette(c("white", "green")),
blues=colorRampPalette(c("white", "blue")))
#coloNames=names(colors)
varNames=c('X','Y','Z')
#cols=c('red','green','blue')

##get feature importance
load(paste0(path,'/result/3Dpred2500/model/gongan1462luhj.model.imp.RData'))
imp=impRes$importance[,1]#/100*255
num=length(imp)/3
xCol=1:num;yCol=(num+1):(num*2);zCol=(num*2+1):(num*3);
imp=data.frame(ximp=imp[xCol],yimp=imp[yCol],zimp=imp[zCol])
min6=round(apply(imp,2,min),2)
max6=round(apply(imp,2,max),2)
#plot and save
for(i in 1:3){
  distMin1=min(abs(min6[i]-seq(0,100,20)))
  #distMin2=min(abs(abs(min6[i])-seq(0,100,20)))
  distMax1=min(abs(max6[i]-seq(0,100,20)))
  #distMax2=max(abs(abs(max6[i])-seq(0,100,20)))
  dist=9
  pdf(file=paste0(resPath,'/',varNames[i],'.bar.pdf'))
  par(pin=c(1,5),cex.axis=1)
  x=rep(1, numN)
  barplot(x,col= colors[[i]](numN),horiz=T,border=NA,space=0,xlim=c(0,1),xaxt ="n",yaxt = "n")
  min=as.character(min6[i]);
  max=as.character(max6[i]);
  if(i==2){max=paste0(max,'0')}
  ##
  axis(2,seq(0,100,20),las=1,font=1,cex.axis=2)
  if ((distMin1 >=dist) | distMin1==0 ){
      axis(2,at=min6[i],labels=paste0('min: ',min,''),las=1,font=2,cex.axis=2)
    }else{
      axis(2,at=min6[i],labels=paste0('min: ',min,'        '),las=1,font=2,cex.axis=2)
  }
  if ((distMax1 >=dist) |distMax1==0){
      axis(2,at=max6[i],labels=paste0('max: ',max,''),las=1,font=2,cex.axis=2)
    }else{
      axis(2,at=max6[i],labels=paste0('max: ',max,'        '),las=1,font=2,cex.axis=2)
  }

  #axis(2,seq(0,100,20),las=1,font=1,cex.axis=2)
  #axis(2,at=maxXYZ[i],labels=paste0('max:',maxXYZ[i]),las=1,font=2,cex.axis=2)
  #text(2,66,labels=maxXYZ[1],pos=2,cex=1,col=cols[i])
  #if(i==1)axis(2,at=c(0,20,40,60,66,80,100),las=1,font=2)
  #if(i==2)axis(2,at=c(0,20,40,60,71,80,100),las=1,font=2)
  #if(i==3)axis(2,at=c(0,20,40,60,80,100),las=1,font=2)
  mtext(paste0('importance of ',varNames[i]),side=3,font=2,line=-3.5,adj=0.48,outer=T,cex=2)
  dev.off()
}
