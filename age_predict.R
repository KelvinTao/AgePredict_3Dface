##
#path='/data/taoxm/gongan/page_methylation'
path='/Users/taoxianming/Documents/research/gongan/page_methylation'
##
resPath=paste0(path,'/result/3Dpred2500');
modelPath=paste0(resPath,'/model');predPath=paste0(resPath,'/pred')
dir.create(modelPath,recursive=T);dir.create(predPath,recursive=T);
###mat of xxyyzz
#load(paste0(path,'/data/luhj3D2500/qincheng2500mat.RData'))#mat
###GPA process
if(F){
rowName=rownames(mat);
###GPA, put L01171 as first target of GPA
library(vegan)
rowOne=mat[1,];mat[1,]=mat[1226,];mat[1226,]=rowOne;
matAll=t(mat)
pNum=7160
xRow=1:pNum;yRow=(pNum+1):(pNum*2);zRow=(pNum*2+1):(pNum*3);
###use 3D L01171 in matTest as the target for procrustes
gpa12=procrustes(cbind(matAll[xRow,1],matAll[yRow,1],matAll[zRow,1]),cbind(matAll[xRow,2],matAll[yRow,2],matAll[zRow,2]))
matXYZ=cbind(gpa12$X,gpa12$Yrot)
target=gpa12$X
#GPA process
for(i in 3:ncol(matAll)){
  gpa=procrustes(target,cbind(matAll[xRow,i],matAll[yRow,i],matAll[zRow,i]))
  matXYZ=cbind(matXYZ,gpa$Yrot)
  print(i)
}
#reshape to XXXYYYZZZ
matGPA=NULL
for(i in seq(1,ncol(matXYZ),3)){
    matGPA=rbind(matGPA,as.vector(matXYZ[,i:(i+2)]))
    print(i)
}
###put L01171 back
rowOne=matGPA[1,];matGPA[1,]=matGPA[1226,];matGPA[1226,]=rowOne;
###
rownames(matGPA)=rowName;
colnames(matGPA)=1:ncol(matGPA)
###
save(matGPA,file=paste0(resPath,'/matGPA_luhj2500xxyyzz.RData'))
}

###train and test set
load(paste0(resPath,'/matGPA_luhj2500xxyyzz.RData'))
##mat
#phe=read.table(paste0(path,'/phe/geno2500_idAge.xls'),head=T)[,-2]
phe=read.csv(paste0(path,'/phe/TangKun3Dand2500.csv'))
phe=phe[,c('luhj3D','Me.id','age')]
###methylation samples
methyPhe=phe[1:1000,]
##
trainPhe=na.omit(phe[1001:nrow(phe),c(1,3)])
testPhe=methyPhe[,c(1,3)]##
###train set 1463 samples
matTrainUse=merge(trainPhe,data.frame(rownames(matGPA),matGPA),by.x=1,by.y=1)
rownames(matTrainUse)=matTrainUse[,1];matTrainUse=matTrainUse[,-1];
colnames(matTrainUse)=c('age',1:(ncol(matTrainUse)-1))
##969 3D test samples, 988 age, 959 both
matTestUse=merge(testPhe,data.frame(rownames(matGPA),matGPA),by.x=1,by.y=1)
idAge=matTestUse[,1:2]
rownames(matTestUse)=matTestUse[,1];matTestUse=matTestUse[,-c(1,2)];
colnames(matTestUse)=1:ncol(matTestUse)

###tuning
library(caret);library(parallel);library(doParallel);
#########functions
source(paste0('/data/taoxm/gongan/page_methylation/script/3Dpredict/function.R'))
###model parameters
modelMethodArray=c('svmRadial')
tuneLength=60;##for parameter tuning, no change
fitControl=getControl(tuneLength)
clusters=10;
###tuning and prediction
set.seed(100);
for (mi in seq_along(modelMethodArray)){
    mat0=na.omit(matTrainUse);##remove NA of age
    modelMethod=modelMethodArray[mi];
    ##training
    print('training')
    print(modelMethod)
    model=tune(mat0,modelMethod,tuneLength,fitControl,clusters)
    result=summary_model(mat0,model,modelMethod);#cor(model$validation$pred[,,as.numeric(result[[1]])],mat0[,1])
    ##save result
    modelFile=paste0(modelPath,'/gongan1463luhj.model.RData')
    save(model,file=modelFile)
    predFile=paste0(predPath,'/gongan1463luhj.pred.RData')
    save(result,file=predFile)
    ##get important score
    impRes=calImp(mat0,model,modelMethod);
    impFile=gsub('RData','imp.RData',modelFile)
    save(impRes,file=impFile)
    print(date());print('finish!')
    ###
    #load(modelFile)
    pred=data.frame(idAge,pred=predict(model,matTestUse))
    write.csv(pred,file=paste0(predPath,'/tune_predBy1463luhj-xyz.csv'),row.names=F,quote=F);
}
