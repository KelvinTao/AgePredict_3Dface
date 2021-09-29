if(F){
mergeMat<-function(pheMat,matUse){
  if(dim(pheMat)[2]>2){pheMat=pheMat[,-1]}
  pheMat=as.data.frame(pheMat)
  pheMat[,2]=as.numeric(pheMat[,2])
  pheMat=pheMat[!is.na(pheMat[,2]),]
  matUse=matUse[,-1]##remove imgid
  #merge phe and matUse: matUse:rownames: ERGOID//imageID,colomns:h s v (or x y z )
  matAll=merge(pheMat,data.frame(rownames(matUse),matUse),by.x=1,by.y=1)#merge by IID or ImageID
  #shuffle 3 time-------
  set.seed(100)
  matAll=matAll[sample(1:nrow(matAll),nrow(matAll)),]
  #IID to rownames and remove
  rownames(matAll)=matAll[,1];matAll=matAll[,-1];#remove IID or ImageID,phe x y z ...
  colnames(matAll)[2:ncol(matAll)]=1:(ncol(matAll)-1)#remember colomn NO of 3D
  return(matAll)
}
mergeMat_hybrid<-function(pheMat,matUse){
  if(dim(pheMat)[2]>2){pheMat=pheMat[,-1]}
  pheMat=as.data.frame(pheMat)
  pheMat[,2]=as.numeric(pheMat[,2])
  pheMat=pheMat[!is.na(pheMat[,2]),]
  #matUse[1:10,1:3]
  #merge phe and matUse: matUse:rownames: ERGOID//imageID,colomns:h s v (or x y z )
  matAll=merge(pheMat,data.frame(rownames(matUse),matUse),by.x=1,by.y=1)#merge by IID or ImageID
  #shuffle 3 time-------
  set.seed(100)
  matAll=matAll[sample(1:nrow(matAll),nrow(matAll)),]
  #IID to rownames and remove
  rownames(matAll)=matAll[,1];matAll=matAll[,-1];#remove IID or ImageID,phe x y z ...
  colnames(matAll)[2:ncol(matAll)]=1:(ncol(matAll)-1)#remember colomn NO of 3D
  return(matAll)
}


###remove highly correlated factors
remCor<-function(matUse,cutOff){
  library(caret)
  phe=matUse[,1];matData=matUse[,-1]
  # find attributes that are highly corrected (ideally >0.9)
  highCor=findCorrelation(cor(matData), cutOff)# need return
  if (length(highCor)==0){return(list(matUseR=cbind(phe,matData),highCor=0))
  }else{return(list(matUseR=cbind(phe,matData[,-highCor]),highCor=highCor))}
}

prepare<-function(pheFile,matFile,ifRemCor,cutCor){ 
  load(pheFile) #pheMat: ID phe 
  load(matFile) #mat: ID xxx...yyy...zzz...or xyzxyz
  mat=mergeMat(pheMat,mat);
  ##remove high correlation features
  if(ifRemCor=='yes')mat=remCor(mat,cutCor)$matUseR
  return(mat)
}
prepare_hybrid<-function(pheFile,matFile1,matFile2,ifRemCor,cutCor){ 
  load(matFile1);mat=mat[,-1];mat1=mat;rm(mat);#IID imgID leave one
  load(matFile2);mat=mat[,-1];mat2=mat;rm(mat);
  mat=merge(data.frame(rownames(mat1),mat1),data.frame(rownames(mat2),mat2),by.x=1,by.y=1)
  rownames(mat)=mat[,1];mat=mat[,-1];
  load(pheFile)
  mat=mergeMat_hybrid(pheMat,mat)
  ##remove high correlation features
  if(ifRemCor=='yes')mat=remCor(mat,max(cutCor))$matUseR
  ##merge hsv and xyz
  return(mat)
}
}
summary<-function(mat,model,modelMethod){
  ##get imgID and match
  matID=data.frame(imgID=rownames(mat),rowIndex=(1:dim(mat)[1]))
  jud=0.00001
  if (modelMethod=='glmnet'){
    alpha = model$bestTune$alpha;lambda= model$bestTune$lambda
    pred0=model$pred
    predAge=pred0[which((abs(pred0$alpha-alpha)<jud)&(abs(pred0$lambda-lambda)<jud)),][,1:3]
  }
  if (modelMethod=='svmRadial'){
    sigma=model$bestTune$sigma;C=model$bestTune$C
    pred0=model$pred
    predAge=pred0[which((abs(pred0$sigma-sigma)<jud)&(abs(pred0$C-C)<jud)),][,1:3]
  }
  pred=merge(matID,predAge,by.x=2,by.y=3)[,-1]
  return(list(bestTune=model$bestTune,pred))
}
summary_pls<-function(mat,model){
  set.seed(100)
  optNcomp=pls::selectNcomp(model,"randomization")
  bestTune=data.frame(optNcomp)
  predAge=model$validation$pred[,,optNcomp]
  matID=data.frame(rownames(mat),obs=mat[,1])
  pred=merge(data.frame(imgID=names(predAge),pred=predAge),matID,by.x=1,by.y=1)
  return(list(bestTune,pred))
}
if(F){
summary_pls<-function(mat,model){##summary by min RMSE
  all=30;rAll=rep(0,all)
  for (ncomp in 1:all){rAll[ncomp]=Metrics::rmse(model$validation$pred[,,ncomp],mat[,1])}
  optNcomp=which(rAll==min(rAll), arr.ind =T)
  bestTune=data.frame(optNcomp)
  predAge=model$validation$pred[,,optNcomp]
  matID=data.frame(rownames(mat),obs=mat[,1])
  pred=merge(data.frame(imgID=names(predAge),pred=predAge),matID,by.x=1,by.y=1)
  return(list(bestTune,pred))
}
}

##hyperparameter tune set up
####seeds
getControl<-function(tuneLength=30){
  kFolds=10
  set.seed(100)
  seeds <- vector(mode = "list", length = (kFolds+1))
  for(i in 1:kFolds) seeds[[i]] <- sample.int(1000, tuneLength)
  ## For the last model:
  seeds[[kFolds+1]] <- sample.int(1000, 1)
  ##running set
  fitControl <- trainControl(method = "cv",
                           number = kFolds,
                           returnResamp='all',
                           savePredictions='all',
                           search = "random",
                           seeds=seeds)
  return(fitControl)
}

training<-function(mat,modelMethod,tuneLength,fitControl){
  set.seed(1000)
  model=train(mat[,-1],mat[,1],
  #model=train(mat[1:200,2:300],mat[1:200,1],
        method=modelMethod,metric='RMSE',
        preProcess=c("center",'scale'),
        trControl=fitControl,
        tuneLength=tuneLength)

  return(model)
}

tune<-function(mat,modelMethod,tuneLength,fitControl,clusters){
  if(modelMethod!='oscorespls'){
    cl= makePSOCKcluster(clusters);registerDoParallel(cl);
    model=training(mat,modelMethod,tuneLength,fitControl)
    stopCluster(cl)# parallel computing
  }else{
    pls.options(parallel = makeCluster(clusters, type = "PSOCK"))# parallel computing
    colnames(mat)[1]='Age'
    model=plsr(Age~.,data=mat,scale=T,method='oscorespls',validation='CV')#scale by sd in each CV segment
    stopCluster(pls.options()$parallel)
  }
  return(model)
}
summary_model<-function(mat,model,modelMethod){
  if(modelMethod!='oscorespls'){
    return(summary(mat,model,modelMethod))
  }else{
    return(summary_pls(mat,model))
  }
}
getImp<-function(model){
  impRes=varImp(model);
  imp=impRes$importance;
  imp$colname=as.numeric(gsub('`','',rownames(imp)))
  impRes$impRes=imp;
  #impSort=sort(imp[,1],decreasing = T,index.return =T)
  #impCol=impSort$ix+1
  return(impRes)
}

calImp<-function(mat,model,modelMethod){
  if(modelMethod=='oscorespls'){
    ncomp=as.numeric(summary_model(mat,model,modelMethod)[[1]]);
    cl= makePSOCKcluster(clusters);registerDoParallel(cl);
    set.seed(100);modelMethod='pls';
    #fitControl=trainControl(method='cv',number=10,returnResamp='all',savePredictions='all');tuneGrid=data.frame(ncomp=ncomp);
    fitControl=trainControl(method='none');tuneGrid=data.frame(ncomp=ncomp);
    model=train(mat[,-1],mat[,1],method=modelMethod,metric='RMSE',preProcess=c("center",'scale'),trControl=fitControl,tuneGrid=tuneGrid);
    stopCluster(cl)# parallel computing
  }      
  impRes=getImp(model);
  return(impRes)
}

