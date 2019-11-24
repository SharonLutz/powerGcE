powerGcE <-
function(nCase=407,nControl=376,MAF=0.49,meanE=0,varE=0.99,beta0=-0.32,betaSNP=0.17,betaE=0.97,betaI=seq(-1,-0.75,by=0.05),nSim=1000,alpha=0.00000005,plot.output=TRUE,plot.name="powerGcE.pdf",seed=1){ 
	
  ####################################
  # input parameters
  ####################################
  # nCase is the number of cases
  # nControl is the number of controls
  # MAF= minor allele frequency for the SNP
  # meanE is the mean of the normally distributed environmental exposure
  # varE is the variance of the normally distributed environmental exposure
  # logit(P(Y=1))=Beta0+BetaSNP*SNP+BetaE*E+BetaI*SNP*E
  
  # nSim is the number of simulations
  # alpha is the alpha level, default=0.05
  
  ####################################
  # Error Checks
  ####################################
  # Check length(betaI)
  if(length(betaI)<2){stop("length(betaI) must be 2 or greater.")}
  if(length(beta0)!=1){stop("length(beta0) must equal 1.")}
  if(length(betaSNP)!=1){stop("length(betaSNP) must equal 1.")}
  if(length(betaE)!=1){stop("length(betaE) must equal 1.")}
  if(length(meanE)!=1){stop("length(meanE) must equal 1.")}
  if(length(varE)!=1){stop("length(varE) must equal 1.")}
  if(length(nCase)!=1){stop("length(nCase) must equal 1.")}
  if(length(nControl)!=1){stop("length(nControl) must equal 1.")}
  if(length(MAF)!=1){stop("length(MAF) must equal 1.")}
  if(length(nSim)!=1){stop("length(nSim) must equal 1.")}
  if(length(alpha)!=1){stop("length(alpha) must equal 1.")}
  
  # Check nCase, nControl, nSim are integers greater than 0
  if(!(nCase>0)){stop("nCase must be greater than 0.")}
  if(!(nControl>0)){stop("nControl must be greater than 0.")}
  if(!(nSim>0)){stop("nSim must be greater than 0.")}
  if(nCase%%1!=0){stop("nCase must be an integer.")}
  if(nControl%%1!=0){stop("nControl must be an integer.")}
  if(nSim%%1!=0){stop("nSim must be an integer.")}
  
  # Check varE > 0 
  if(!(varE>0)|varE==0){stop("varE must be greater than 0.")}
  
  # Check alpha>0 & alpha<1 and MAF
  if(alpha<0 | alpha>1){stop("alpha must be between 0 and 1.")}
  if(MAF<0 | MAF>1){stop("MAF must be between 0 and 1.")}
  
  ####################################
  # Store Results
  #################################### 
  rejectH0<-matrix(0,nrow=length(betaI),ncol=2)
  colnames(rejectH0)<-c("BetaI","power")
  rejectH0[,"BetaI"]<-betaI
    
     set.seed(seed)
  ####################################
  # Run Simulations
  #################################### 
  for(GLOBALVAR in 1:nSim){
    if(floor(GLOBALVAR/500)==ceiling(GLOBALVAR/500)){print(paste(GLOBALVAR,"of",nSim,"Simulations"))}
    ####################################
    # Generate Data
    ####################################
    # CYCLE through values of betaI
    for(bb in 1:length(betaI)){
      betaIv<-betaI[bb]
      
      ####################################
      # simulate data
      ####################################
      # overall sample size
      n<-nCase+nControl
      
      # generate the matrix of SNPs
      SNP<-rbinom(n,2,MAF)
      
      # generate the environment E
      E<- rnorm(n,mean=meanE,sd=sqrt(varE)) 
      
      # generate the binary outcome Y
      logitP<-beta0+betaSNP*SNP+betaE*E+betaIv*SNP*E
      Phat<-exp(logitP)/(1+exp(logitP))
      
      matY<-cbind(SNP,E,rbinom(n,1,Phat))
      while(length(matY[matY[,3]==1,3])<nCase|length(matY[matY[,3]==0,3])<nControl){
      matY<-rbind(matY,cbind(SNP,E,rbinom(n,1,Phat)))
      }
      matY1<-matY[matY[,3]==1,]
      matY0<-matY[matY[,3]==0,]
      matYr<-rbind(matY1[1:nCase,],matY0[1:nControl,])
      
      y<-matYr[,3]
      snp<-matYr[,1]
      e<-matYr[,2]
      
      ####################################
      # logistic regression for interaction
      ####################################
      modelI<-summary(glm(y~snp+e+snp*e,family=binomial()))$coef
      if(nrow(modelI)<(4)){stop("Not valid estimates in logistic regression. This could be due to lack of vairability.")}
      
      if(modelI[4,4]<alpha){rejectH0[bb,"power"]<-rejectH0[bb,"power"]+1}
        
      ####################################
      # Compile results 
      ####################################
      }#end betaI loop
      }#end globalvar
  
     rejectMat<-cbind(rejectH0[,"BetaI"],rejectH0[,"power"]/nSim)
  colnames(rejectMat)<-c("BetaI","power")
  
  ####################################
  # Create plot
  #################################### 
    rejectS<-rejectMat[order(rejectMat[,"BetaI"],decreasing=TRUE),]
    
  if(plot.output==TRUE){
  pdf(plot.name)
  plot(-1,-1,xlim=c(min(betaI),max(betaI)),ylim=c(0,1),xlab=expression(beta[I]),ylab="Emperical Power",main="")
    lines(rejectS[,1],rejectS[,2])
  dev.off()
  }
  ####################################
  # End function
  ####################################
  list(rejectMat)}
