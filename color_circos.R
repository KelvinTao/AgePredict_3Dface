resPath='/Users/taoxianming/Documents/research/gongan/page_methylation/3Dface'
##color circos
library(circlize)
## produce data
num=10#
NO=num-1
x=NULL;for(i in 0:NO){x=c(x,rep(i,num*num))}
y=NULL;
for(i0 in 0:NO){
	y0=NULL
	for(i in 0:NO){
		y0=c(y0,rep(i,num))
	}
	y=c(y,y0)
}
z=rep(rep(0:NO,num),num)
xyz=data.frame(x,y,z)
##calculate color
reds=colorRampPalette(c("white", "red"));
greens=colorRampPalette(c("white", "green"));
blues=colorRampPalette(c("white", "blue"));
##
xColor=reds(num);yColor=greens(num);zColor=blues(num)
colSum=col2rgb(xColor)+col2rgb(yColor)+col2rgb(zColor)
mmax=max(colSum);mmin=min(colSum);len=mmax-mmin;
##calculate each point color
getColor<-function(x){
	colSum=col2rgb(xColor[x[1]])+col2rgb(yColor[x[2]])+col2rgb(zColor[x[3]])
    return(rgb(t((colSum-mmin)/len)))
}
xyzColor=apply(xyz+1,1,getColor)

###circos heatmap
factors=as.factor(t(matrix(rep(0:NO,num),num,num)))
row=col=num;
col_mat_list=vector("list", num)##each red
num2D=num*num;
for(i in 1:num){
  seqUse=(i*num2D):((i-1)*num2D+1)
  mat0=matrix(xyzColor[seqUse],row,col)
  mat=mat0
    for(k in 1:col){
      mat[,k]=mat0[,col+1-k]
    }
  col_mat_list[[i]]=mat
}

##circos heatmap plot
pdf(file=paste0(resPath,'/XYZ_RGB_mixed.color.circos.pdf'))
circos.par(track.height = 0.35,cell.padding = c(0, 0, 0, 0), gap.after=5.9,
	start.degree=105,track.margin=c(0.05,0.05))#
circos.initialize(factors,xlim = cbind(rep(0,num), table(factors)))
circos.track(ylim=c(0,row), bg.border = NA, panel.fun = function(x, y){
	sector.index = as.numeric(CELL_META$sector.index)+1
  #if(sector.index==7)break
	col_mat=col_mat_list[[sector.index]]
    #add color
    for(i in 1:row){
        circos.rect(1:col-1, rep(row-i, col), 
            1:col, rep(row-i+1, col), 
            border = col_mat[i, ], col = col_mat[i, ])
    }
      xi=sector.index
      if(sector.index!=7)circos.text(5, 11.2, paste0('Red ',xi*10),cex=1.5,font=2,col='red')
      if(sector.index==7)circos.text(5, 11.2, paste0('Red(Max) ',xi*10),cex=1.5,font=2,col='red')
      circos.axis('bottom',major.at=c(0.1,4.5,10),labels=c(10,50,100),
      	labels.cex=1,direction = "inside",major.tick=F,major.tick.length=0.1,
      	labels.col='green')#
      circos.yaxis(at=c(0.5,4.5,9.5),labels=c(10,50,100),
      	labels.cex=1,tick.length=0.1,,tick=F,labels.col='blue')
  }
)
title('Colors of merged importances of X (Red), Y (Green) and Z (Blue)',cex=3)
circos.clear()
dev.off()
