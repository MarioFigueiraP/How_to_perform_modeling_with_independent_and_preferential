source("LibrariesAndAuxiliaryFunctions.R")
source("IndependentModel.R")
source("DependentModel.R")
source("DependentModelComplex.R")
# source("PointProcessInference.R")

library("gridExtra")
library("cowplot")

combinatorialFunction <- function(list){
  DF.old <- list
  DF.new <- list
  for(i in 2:length(A)){
    DF.new[[i-1]] <- rep(DF.old[[i-1]], length(DF.old[[i]]))
    DF.new[[i]] <- rep(DF.old[[i]], each=length(DF.old[[i-1]]))
    DF.old <- DF.new
  }
  DF <- as.data.frame(DF.new)
  return(DF)
}

range_Spatial <- c(0.2, 0.5, 0.8)
sigma_Spatial <- c(0.5)
# nTot <- c(80, 150, 200, 400)
# nProp <- c(0.1, 0.25, 0.5, 0.75, 0.9)
nTot <- c(60, 100, 160, 200)
nProp <- c(0.1, 0.25, 0.5, 0.75, 0.9)

# formula <- c("0.5 + (x) + 0.8*(y-0.5)**2")

formula <- c("0.5 + (x) + 0.8*(y-0.5)**2")

# formula <- c("0.5 + (x) + 0.8*(y-0.5)**2")#,
             # "-1+1.7*(x)**2 + 0.5*((y-0.5)**2+0.3)**-1",
             # "0.5 + (x) + 0.8*(y-0.5)**2")

# formula <- c("0.5 + (x) + 0.8*(y-0.5)**2")

A <- list(range=range_Spatial, sigma=sigma_Spatial, nTot=nTot, nProp=nProp, formula=formula)
DF <- combinatorialFunction(list=A)
DF <- lapply(X=DF, MARGIN=2, FUN=rep, 4)
# DF <- as.data.frame(apply(X=DF, MARGIN=2, FUN=rep, 5))

# limlattice <- c(0,1)
# lengthlattice <- 25
# xlattice0 <-  seq(limlattice[1], limlattice[2], length.out=lengthlattice)
# ylattice0 <- seq(limlattice[1], limlattice[2], length.out=lengthlattice)
# lattice0 <-  inla.mesh.lattice(xlattice0, ylattice0)
# meshInf <-  inla.mesh.2d(loc=lattice0$loc, max.edge=c(0.05,0.1))

rscale <- rep(1.6, length(DF[[1]]))
# rscale[DF$sigma==1&DF$range==0.2] <- 0.8
# rscale[DF$sigma==1&DF$range==0.5] <- 1.2
# rscale[DF$sigma==1&DF$range==0.8] <- 1.6

betaModel <- c(-2,2)

WAIC <- data.frame(IND=NA, DEP=NA, MIX=NA)
RMSEMean <- data.frame(IND=NA, DEP=NA, MIX=NA)
IntegratedValues <- data.frame(SIM=NA, IND=NA, DEP=NA, MIX=NA)

t1 <- Sys.time()
# for(i in 1:nrow(DF)){
for(i in 122:length(DF[[1]])){#seq_along(DF[[1]])){#length(DF[[1]])
  # i <- 1
  print(i)
  globalseed <- i
  Sim <- Simulation(globalseed=globalseed, range0=DF$range[i], sigma0=DF$sigma[i], nsamplesInd=DF$nTot[i]*DF$nProp[i], nsamplesDep=DF$nTot[i]*(1-DF$nProp[i]), rscale=rscale[i], cov.formula=DF$formula[i], beta.vec=betaModel)
  Sim$DFSampleInd$kind <- "ind"; Sim$DFSampleDep$kind <- "dep"
  DFSample <- rbind(Sim$DFSampleInd, Sim$DFSampleDep)
  # ggTestInd <- ggplot() + geom_tile(data=Sim$DFSim, mapping=aes(x=x,y=y,fill=ysim)) + scale_fill_viridis_c(option="turbo") +
  #   geom_point(data=Sim$DFSampleInd, mapping=aes(x=x,y=y))
  # print(ggTestInd)
  # ggTestDep <- ggplot() + geom_tile(data=Sim$DFSim, mapping=aes(x=x,y=y,fill=ysim)) + scale_fill_viridis_c(option="turbo") +
  #   geom_point(data=Sim$DFSampleDep, mapping=aes(x=x,y=y))
  # print(ggTestDep)
  ggTestDep <- ggplot() + geom_tile(data=Sim$DFSim, mapping=aes(x=x,y=y,fill=ysim)) + scale_fill_viridis_c(option="turbo") +
    geom_point(data=DFSample, mapping=aes(x=x,y=y))
  print(ggTestDep)

  meshInf <-  inla.mesh.2d(loc=rbind(Sim$DFSampleInd, Sim$DFSampleDep)[,1:2],
                           loc.domain = matrix(c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE),
                           max.edge=c(0.05,0.1), cutoff=0.01, offset=c(0.15,0.3))

  priormatern <- list(range=c(0.5,0.5), sigma=c(1,0.5)) #This First
  prior.fixed <- list(meanInterceptor=0, meanCov=0, precInterceptor=0.001, precCov=0.001)
  IndependentModelPCMatern <- IndependentModel(dataSample=DFSample, dataSim=Sim$DFSim, mesh=meshInf,
                                               priormaterntype="pcmatern", priormatern=priormatern,
                                               feedbacktype="moments", prior.fixed=prior.fixed)
  
  # DFIndPCMatern <- data.frame(x=Sim$DFSim$x, y=Sim$DFSim$y, ysim=IndependentModelPCMatern$IndepedentModel$summary.fitted.values[IndependentModelPCMatern$index.pred, "mean"])
  # ggplot() + geom_tile(data=DFIndPCMatern, mapping=aes(x=x,y=y,fill=ysim)) + scale_fill_viridis_c(option="turbo")
  
  priormatern <- list(range=c(0.5,0.5), sigma=c(1,0.5)) #This Second
  prior.fixed <- list(meanInterceptor=0, meanCov=0, precInterceptor=0.001, precCov=0.001)
  DependentModelPCMatern <- DependentModel(dataSample=DFSample, dataSim=Sim$DFSim, mesh=meshInf,
                                   priormaterntype="pcmatern", priormatern=priormatern,
                                   feedbacktype="moments", prior.fixed=prior.fixed,
                                   cov.formula=DF$formula[i], prior.scale.factor=c(0, 1))
  
  priormatern <- list(range=c(0.5,0.5), sigma=c(1,0.5)) #This third
  prior.fixed <- list(meanInterceptor=0, meanCov=0, precInterceptor=0.001, precCov=0.001)
  DependentModelPCMaternComplex <- DependentModelComplex(dataSample=DFSample, dataSim=Sim$DFSim, mesh=meshInf,
                                           priormaterntype="pcmatern", priormatern=priormatern,
                                           feedbacktype="moments", prior.fixed=prior.fixed,
                                           cov.formula=DF$formula[i], prior.scale.factor=c(0, 1))
  
  # DFDepPCMatern <- data.frame(x=Sim$DFSim$x, y=Sim$DFSim$y, ysim=DependentModelPCMatern$DepedentModel$summary.fitted.values[DependentModelPCMatern$index.pred, "mean"])
  # ggplot() + geom_tile(data=DFDepPCMatern, mapping=aes(x=x,y=y,fill=ysim)) + scale_fill_viridis_c(option="turbo")
  
  WAIC[i,1] <- sum(IndependentModelPCMatern$IndepedentModel$waic$local.waic[which(IndependentModelPCMatern$IndepedentModel$dic$family==1)])
  WAIC[i,2] <- sum(DependentModelPCMatern$DepedentModel$waic$local.waic[which(DependentModelPCMatern$DepedentModel$dic$family==1)])
  WAIC[i,3] <- sum(DependentModelPCMaternComplex$DepedentModel$waic$local.waic[which(DependentModelPCMaternComplex$DepedentModel$dic$family==1)])
  
  RMSEMean[i,1] <- IndependentModelPCMatern$RMSEmean
  RMSEMean[i,2] <- DependentModelPCMatern$RMSEmean
  RMSEMean[i,3] <- DependentModelPCMaternComplex$RMSEmean
  
  IntegratedValues[i,1] <- sum(Sim$DFSim$ysim)/nrow(Sim$DFSim)
  IntegratedValues[i,2] <- sum(IndependentModelPCMatern$IndepedentModel$summary.fitted.values[IndependentModelPCMatern$index.pred, "mean"])/nrow(Sim$DFSim)
  IntegratedValues[i,3] <- sum(DependentModelPCMatern$DepedentModel$summary.fitted.values[DependentModelPCMatern$index.pred, "mean"])/nrow(Sim$DFSim)
  IntegratedValues[i,4] <- sum(DependentModelPCMaternComplex$DepedentModel$summary.fitted.values[DependentModelPCMaternComplex$index.pred, "mean"])/nrow(Sim$DFSim)
}

 condition <- function(x,i){return(DF[[x]]==unique(DF[[x]])[i])}
mean(DIC[condition("nTot", 1)&condition("nProp", 3),1]-DIC[condition("nTot", 1)&condition("nProp", 3),2]<(-3))
mean(DIC[condition("nTot", 1)&condition("nProp", 1),1]-DIC[condition("nTot", 1)&condition("nProp", 1),2]>3)

boxplot(DIC[condition("nTot", 4),1]-DIC[condition("nTot", 4),2])

t2 <- Sys.time()

# saveRDS(DIC,"DIC_2.rds")
saveRDS(WAIC,"WAIC.rds")
saveRDS(RMSEMean, "RMSEmean.rds")
saveRDS(IntegratedValues, "IntegratedValues.rds")
# saveRDS(RMSEMedian, "RMSEmedian_2.rds")
# saveRDS(MAPEMean, "MAPEMean_2.rds")
# saveRDS(MAPEMedian, "MAPEMedian_2.rds")

DIC <- readRDS("DIC.rds")
WAIC <- readRDS("WAIC.rds")
RMSEMean <- readRDS("RMSEMean.rds")
RMSEMedian <- readRDS("RMSEMedian.rds")
MAPEMean <- readRDS("MAPEMean.rds")
MAPEMedian <- readRDS("MAPEMedian.rds")

ggplot(RMSEMean[order(RMSEMean$IND),]) + geom_line(aes(x=1:nrow(RMSEMean),y=IND)) +
  geom_line(aes(x=1:nrow(RMSEMean),y=DEP), col="red") +
  geom_line(aes(x=1:nrow(RMSEMean),y=MIX), col="blue")

# sum(DFRMSEMedian$FDEPM[73:108]<DFRMSEMedian$DEPM[73:108])/nrow(DFRMSEMean[73:108,])

k <- 3
range_Spatial[k]
sum(DFRMSEMedian$FDEPM[DF$range[1:nrow(DFRMSEMean)]!=range_Spatial[k]]<DFRMSEMedian$DEPM[DF$range[1:nrow(DFRMSEMean)]!=range_Spatial[k]])/nrow(DFRMSEMean[DF$range[1:nrow(DFRMSEMean)]!=range_Spatial[k],])

k <- 3
sigma_Spatial[k]
sum(DFRMSEMedian$FDEPM[DF$sigma[1:nrow(DFRMSEMean)]==sigma_Spatial[k]]<DFRMSEMedian$DEPM[DF$sigma[1:nrow(DFRMSEMean)]==sigma_Spatial[k]])/nrow(DFRMSEMean[DF$sigma[1:nrow(DFRMSEMean)]==sigma_Spatial[k],])

k <- 4
nTotal[k]
sum(DFRMSEMedian$FDEPM[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[k]]<DFRMSEMedian$DEPM[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[k]])/nrow(DFRMSEMean[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[k],])
sum(DFRMSEMedian$FDEPPCM[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[k]]<DFRMSEMedian$DEPPCM[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[k]])/nrow(DFRMSEMean[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[k],])

k <- 3
formula[k]
sum(DFRMSEMedian$FDEPM[DF$formula[1:nrow(DFRMSEMean)]==formula[k]]<DFRMSEMedian$DEPM[DF$formula[1:nrow(DFRMSEMean)]==formula[k]])/nrow(DFRMSEMean[DF$formula[1:nrow(DFRMSEMean)]==formula[k],])


sum(DFRMSEMedian$FDEPPCM<DFRMSEMedian$DEPM)/nrow(DFRMSEMean)


sum(DFRMSEMedian$FINDM[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[1]]<DFRMSEMedian$INDM[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[1]])/nrow(DFRMSEMean[DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[1],])


DF$nTotal[1:nrow(DFRMSEMean)]==nTotal[1]

# saveRDS(object=DF, file="./DFFinal3.rds")
# saveRDS(object=DFRMSEMean, file="./DFRMSEMeanFinal3.rds")
# saveRDS(object=DFRMSEMedian, file="./DFRMSEMedianFinal3.rds")
# saveRDS(object=DFBiasMean, file="./DFBiasMeanFinal3.rds")
# saveRDS(object=DFBiasMedian, file="./DFBiasMedianFinal3.rds")

DF <- readRDS(file="./DFFinal3.rds")
DFRMSEMean <- readRDS(file="./DFRMSEMeanFinal3.rds")
DFRMSEMedian <- readRDS(file="./DFRMSEMedianFinal3.rds")
DFBiasMean <- readRDS(file="./DFBiasMeanFinal3.rds")
DFBiasMedian <- readRDS(file="./DFBiasMedianFinal3.rds")


# D1 <- readRDS(file="./DFRMSEMean1.rds")
# D2 <- readRDS(file="./DFRMSEMedian1.rds")

apply(DFRMSEMean, FUN=mean, MARGIN=2)
apply(DFRMSEMean, FUN=median, MARGIN=2)
apply(DFRMSEMedian, FUN=mean, MARGIN=2)
apply(DFRMSEMedian, FUN=median, MARGIN=2)

# Example plots

### Intercept

ggplot() + 
  geom_line(data=as.data.frame(IndependentModelMatern$IndepedentModel$marginals.fixed$Intercept), mapping=aes(x=x,y=y, group=1, col="IM-EN")) +
  geom_line(data=as.data.frame(FeedbackIndependentModelMatern$IndepedentModel$marginals.fixed$Intercept), mapping=aes(x=x,y=y, group=1, col="IFM-EN")) +
  geom_vline(xintercept=betaModel[1], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-EN', 'IFM-EN', 'Real Value'),
    values=c('IM-EN'='black', 'IFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = "Intercept", y=expression(pi*"(Intercept)")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(IndependentModelPCMatern$IndepedentModel$marginals.fixed$Intercept), mapping=aes(x=x,y=y, group=1, col="IM-PC")) +
  geom_line(data=as.data.frame(FeedbackIndependentModelPCMatern$IndepedentModel$marginals.fixed$Intercept), mapping=aes(x=x,y=y, group=1, col="IFM-PC")) +
  geom_vline(xintercept=betaModel[1], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-PC', 'IFM-PC', 'Real Value'),
    values=c('IM-PC'='black', 'IFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = "Intercept", y=expression(pi*"(Intercept)")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(DependentModelMatern$DepedentModel$marginals.fixed$Interceptor), mapping=aes(x=x,y=y, group=1, col="PM-EN")) +
  geom_line(data=as.data.frame(FeedbackDependentModelMatern$DepedentModel$marginals.fixed$Interceptor), mapping=aes(x=x,y=y, group=1, col="PFM-EN")) +
  geom_vline(xintercept=betaModel[1], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-EN', 'PFM-EN', 'Real Value'),
    values=c('PM-EN'='black', 'PFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = "Intercept", y=expression(pi*"(Intercept)")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(DependentModelPCMatern$DepedentModel$marginals.fixed$Interceptor), mapping=aes(x=x,y=y, group=1, col="PM-PC")) +
  geom_line(data=as.data.frame(FeedbackDependentModelPCMatern$DepedentModel$marginals.fixed$Interceptor), mapping=aes(x=x,y=y, group=1, col="PFM-PC")) +
  geom_vline(xintercept=betaModel[1], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-PC', 'PFM-PC', 'Real Value'),
    values=c('PM-PC'='black', 'PFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = "Intercept", y=expression(pi*"(Intercept)")) + theme_bw()

### Covariate

ggplot() + 
  geom_line(data=as.data.frame(IndependentModelMatern$IndepedentModel$marginals.fixed$Cov), mapping=aes(x=x,y=y, group=1, col="IM-EN")) +
  geom_line(data=as.data.frame(FeedbackIndependentModelMatern$IndepedentModel$marginals.fixed$Cov), mapping=aes(x=x,y=y, group=1, col="IFM-EN")) +
  geom_vline(xintercept=betaModel[2], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-EN', 'IFM-EN', 'Real Value'),
    values=c('IM-EN'='black', 'IFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = "Cov", y=expression(pi*"(Cov)")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(IndependentModelPCMatern$IndepedentModel$marginals.fixed$Cov), mapping=aes(x=x,y=y, group=1, col="IM-PC")) +
  geom_line(data=as.data.frame(FeedbackIndependentModelPCMatern$IndepedentModel$marginals.fixed$Cov), mapping=aes(x=x,y=y, group=1, col="IFM-PC")) +
  geom_vline(xintercept=betaModel[2], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-PC', 'IFM-PC', 'Real Value'),
    values=c('IM-PC'='black', 'IFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = "Cov", y=expression(pi*"(Cov)")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(DependentModelMatern$DepedentModel$marginals.fixed$Cov), mapping=aes(x=x,y=y, group=1, col="PM-EN")) +
  geom_line(data=as.data.frame(FeedbackDependentModelMatern$DepedentModel$marginals.fixed$Cov), mapping=aes(x=x,y=y, group=1, col="PFM-EN")) +
  geom_vline(xintercept=betaModel[2], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-EN', 'PFM-EN', 'Real Value'),
    values=c('PM-EN'='black', 'PFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = "Cov", y=expression(pi*"(Cov)")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(DependentModelPCMatern$DepedentModel$marginals.fixed$Cov), mapping=aes(x=x,y=y, group=1, col="PM-PC")) +
  geom_line(data=as.data.frame(FeedbackDependentModelPCMatern$DepedentModel$marginals.fixed$Cov), mapping=aes(x=x,y=y, group=1, col="PFM-PC")) +
  geom_vline(xintercept=betaModel[2], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-PC', 'PFM-PC', 'Real Value'),
    values=c('PM-PC'='black', 'PFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = "Cov", y=expression(pi*"(Cov)")) + theme_bw()

### Precision

ggplot() + 
  geom_line(data=as.data.frame(IndependentModelMatern$IndepedentModel$marginals.hyperpar[[1]]), mapping=aes(x=x,y=y, group=1, col="IM-EN")) +
  geom_line(data=as.data.frame(FeedbackIndependentModelMatern$IndepedentModel$marginals.hyperpar[[1]]), mapping=aes(x=x,y=y, group=1, col="IFM-EN")) +
  geom_vline(xintercept=100, color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-EN', 'IFM-EN', 'Real Value'),
    values=c('IM-EN'='black', 'IFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(phi), y=expression(pi*"("*phi*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(IndependentModelPCMatern$IndepedentModel$marginals.hyperpar[[1]]), mapping=aes(x=x,y=y, group=1, col="IM-PC")) +
  geom_line(data=as.data.frame(FeedbackIndependentModelPCMatern$IndepedentModel$marginals.hyperpar[[1]]), mapping=aes(x=x,y=y, group=1, col="IFM-PC")) +
  geom_vline(xintercept=100, color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-PC', 'IFM-PC', 'Real Value'),
    values=c('IM-PC'='black', 'IFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(phi), y=expression(pi*"("*phi*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(DependentModelMatern$DepedentModel$marginals.hyperpar[[1]]), mapping=aes(x=x,y=y, group=1, col="PM-EN")) +
  geom_line(data=as.data.frame(FeedbackDependentModelMatern$DepedentModel$marginals.hyperpar[[1]]), mapping=aes(x=x,y=y, group=1, col="PFM-EN")) +
  geom_vline(xintercept=100, color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-EN', 'PFM-EN', 'Real Value'),
    values=c('PM-EN'='black', 'PFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(phi), y=expression(pi*"("*phi*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(DependentModelPCMatern$DepedentModel$marginals.hyperpar[[1]]), mapping=aes(x=x,y=y, group=1, col="PM-PC")) +
  geom_line(data=as.data.frame(FeedbackDependentModelPCMatern$DepedentModel$marginals.hyperpar[[1]]), mapping=aes(x=x,y=y, group=1, col="PFM-PC")) +
  geom_vline(xintercept=100, color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-PC', 'PFM-PC', 'Real Value'),
    values=c('PM-PC'='black', 'PFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(phi), y=expression(pi*"("*phi*")")) + theme_bw()

### Sigma

ggplot() + 
  geom_line(data=as.data.frame(inla.tmarginal(function(x) 1*exp(x), IndependentModelMatern$IndepedentModel$marginals.hyperpar[[2]])), mapping=aes(x=x,y=y, group=1, col="IM-EN")) +
  geom_line(data=as.data.frame(inla.tmarginal(function(x) 1*exp(x),FeedbackIndependentModelMatern$IndepedentModel$marginals.hyperpar[[2]])), mapping=aes(x=x,y=y, group=1, col="IFM-EN")) +
  geom_vline(xintercept=DF$sigma[i], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-EN', 'IFM-EN', 'Real Value'),
    values=c('IM-EN'='black', 'IFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(sigma), y=expression(pi*"("*sigma*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(IndependentModelPCMatern$IndepedentModel$marginals.hyperpar[[3]]), mapping=aes(x=x,y=y, group=1, col="IM-PC")) +
  geom_line(data=as.data.frame(FeedbackIndependentModelPCMatern$IndepedentModel$marginals.hyperpar[[3]]), mapping=aes(x=x,y=y, group=1, col="IFM-PC")) +
  geom_vline(xintercept=DF$sigma[i], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-PC', 'IFM-PC', 'Real Value'),
    values=c('IM-PC'='black', 'IFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(sigma), y=expression(pi*"("*sigma*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(inla.tmarginal(function(x) 1*exp(x), DependentModelMatern$DepedentModel$marginals.hyperpar[[2]])), mapping=aes(x=x,y=y, group=1, col="PM-EN")) +
  geom_line(data=as.data.frame(inla.tmarginal(function(x) 1*exp(x), FeedbackDependentModelMatern$DepedentModel$marginals.hyperpar[[2]])), mapping=aes(x=x,y=y, group=1, col="PFM-EN")) +
  geom_vline(xintercept=DF$sigma[i], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-EN', 'PFM-EN', 'Real Value'),
    values=c('PM-EN'='black', 'PFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(sigma), y=expression(pi*"("*sigma*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(DependentModelPCMatern$DepedentModel$marginals.hyperpar[[3]]), mapping=aes(x=x,y=y, group=1, col="PM-PC")) +
  geom_line(data=as.data.frame(FeedbackDependentModelPCMatern$DepedentModel$marginals.hyperpar[[3]]), mapping=aes(x=x,y=y, group=1, col="PFM-PC")) +
  geom_vline(xintercept=DF$sigma[i], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-PC', 'PFM-PC', 'Real Value'),
    values=c('PM-PC'='black', 'PFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(sigma), y=expression(pi*"("*sigma*")")) + theme_bw()

### Range

ggplot() + 
  geom_line(data=as.data.frame(inla.tmarginal(function(x) 0.2*exp(x), IndependentModelMatern$IndepedentModel$marginals.hyperpar[[3]])), mapping=aes(x=x,y=y, group=1, col="IM-EN")) +
  geom_line(data=as.data.frame(inla.tmarginal(function(x) 0.2*exp(x),FeedbackIndependentModelMatern$IndepedentModel$marginals.hyperpar[[3]])), mapping=aes(x=x,y=y, group=1, col="IFM-EN")) +
  geom_vline(xintercept=DF$range[i], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-EN', 'IFM-EN', 'Real Value'),
    values=c('IM-EN'='black', 'IFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(rho), y=expression(pi*"("*rho*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(IndependentModelPCMatern$IndepedentModel$marginals.hyperpar[[2]]), mapping=aes(x=x,y=y, group=1, col="IM-PC")) +
  geom_line(data=as.data.frame(FeedbackIndependentModelPCMatern$IndepedentModel$marginals.hyperpar[[2]]), mapping=aes(x=x,y=y, group=1, col="IFM-PC")) +
  geom_vline(xintercept=DF$range[i], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('IM-PC', 'IFM-PC', 'Real Value'),
    values=c('IM-PC'='black', 'IFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(rho), y=expression(pi*"("*rho*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(inla.tmarginal(function(x) 0.2*exp(x), DependentModelMatern$DepedentModel$marginals.hyperpar[[3]])), mapping=aes(x=x,y=y, group=1, col="PM-EN")) +
  geom_line(data=as.data.frame(inla.tmarginal(function(x) 0.2*exp(x), FeedbackDependentModelMatern$DepedentModel$marginals.hyperpar[[3]])), mapping=aes(x=x,y=y, group=1, col="PFM-EN")) +
  geom_vline(xintercept=DF$range[i], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-EN', 'PFM-EN', 'Real Value'),
    values=c('PM-EN'='black', 'PFM-EN'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(rho), y=expression(pi*"("*rho*")")) + theme_bw()

ggplot() + 
  geom_line(data=as.data.frame(DependentModelPCMatern$DepedentModel$marginals.hyperpar[[2]]), mapping=aes(x=x,y=y, group=1, col="PM-PC")) +
  geom_line(data=as.data.frame(FeedbackDependentModelPCMatern$DepedentModel$marginals.hyperpar[[2]]), mapping=aes(x=x,y=y, group=1, col="PFM-PC")) +
  geom_vline(xintercept=DF$range[i], color="red") +
  scale_color_manual(
    name='Model',
    breaks=c('PM-PC', 'PFM-PC', 'Real Value'),
    values=c('PM-PC'='black', 'PFM-PC'='blue', 'Real Value'='red')) +
  theme(legend.position="right") + labs(x = expression(rho), y=expression(pi*"("*rho*")")) + theme_bw()


# Bias

mean(Sim$DFSim$ysim-DependentModelMatern$DepedentModel$summary.fitted.values[DependentModelMatern$index.pred,"mean"])
mean(Sim$DFSim$ysim-FeedbackDependentModelMatern$DepedentModel$summary.fitted.values[FeedbackDependentModelMatern$index.pred,"mean"])
mean(Sim$DFSim$ysim-DependentModelPCMatern$DepedentModel$summary.fitted.values[DependentModelPCMatern$index.pred,"mean"])
mean(Sim$DFSim$ysim-FeedbackDependentModelPCMatern$DepedentModel$summary.fitted.values[FeedbackDependentModelPCMatern$index.pred,"mean"])

plot(inla.tmarginal(function(x) 1*exp(x),DependentModelMatern$DepedentModel$marginals.hyperpar$`Theta2 for spatial`), type="l")
lines(inla.tmarginal(function(x) 1*exp(x),data.frame(x=seq(0.,7, length.out=100), y=dnorm(seq(0,7, length.out=100), mean=2.112771, sd=0.6555))), col="red")

lines(DependentModelPCMatern$DepedentModel$marginals.hyperpar$`Stdev for spatial`, col="red")

plot(inla.tmarginal(function(x) 0.2*exp(x),DependentModelMatern$DepedentModel$marginals.hyperpar$`Theta2 for spatial`), type="l")
lines(DependentModelPCMatern$DepedentModel$marginals.hyperpar$`Range for spatial`, col="red")


plot(Sim$DFSim$ysim-DependentModelMatern$DepedentModel$summary.fitted.values[DependentModelPCMatern$index.pred,"mean"])
plot(Sim$DFSim$ysim-FeedbackDependentModelMatern$DepedentModel$summary.fitted.values[FeedbackDependentModelPCMatern$index.pred,"mean"])

plot(Sim$DFSim$ysim-IndependentModelPCMatern$IndepedentModel$summary.fitted.values[DependentModelPCMatern$index.pred,"mean"])
plot(Sim$DFSim$ysim-FeedbackDependentModelPCMatern$DepedentModel$summary.fitted.values[FeedbackDependentModelPCMatern$index.pred,"mean"])




apply(DFRMSEMean, FUN=mean, MARGIN=2)
apply(DFRMSEMean, FUN=median, MARGIN=2)
apply(DFRMSEMedian, FUN=mean, MARGIN=2)
apply(DFRMSEMedian, FUN=median, MARGIN=2)


saveRDS(object=DFRMSEMean, file="./DFRMSEMean1.rds")
saveRDS(object=DFRMSEMedian, file="./DFRMSEMedian1.rds")

saveRDS(object=DFBiasMean, file="./DFBiasMean1.rds")
saveRDS(object=DFBiasMedian, file="./DFBiasMedian1.rds")

DFMean <- readRDS("./DFRMSEMean1.rds")
DFMedian <- readRDS("./DFRMSEMedian2.rds")


DFMean <- DFMean[1:54,]
DFMedian <- DFMedian[1:54,]

apply(DFMedian, FUN=median, MARGIN=2)


# obMean <- readRDS("./DFRMSEMean.rds")

plot(inla.tmarginal(function(x) 1*exp(x), DependentModelMatern$DepedentModel$marginals.hyperpar$`Theta1 for spatial`), type="l")
lines(DependentModelPCMatern$DepedentModel$marginals.hyperpar$`Stdev for spatial`, col="red")

plot(inla.tmarginal(function(x) 0.2*exp(x), DependentModelMatern$DepedentModel$marginals.hyperpar$`Theta2 for spatial`), type="l")
lines(DependentModelPCMatern$DepedentModel$marginals.hyperpar$`Range for spatial`, col="red")

CorrectedDFmean <- DFRMSEMean[-as.vector(apply(X=DFRMSEMean, FUN=which.max, MARGIN=2)),]
CorrectedDFmedian <- DFRMSEMean[-as.vector(apply(X=DFRMSEMedian, FUN=which.max, MARGIN=2)),]
CorrectedDFmedian <- CorrectedDFmedian[order(CorrectedDFmedian$INDM),]
CorrectedDFmedian$indx <- 1:nrow(CorrectedDFmedian)

ggplot() + geom_tile(data=Sim$DFSim, mapping=aes(x=x, y=y, fill=ysim)) + 
  geom_point(data=Sim$DFSampleDep, mapping=aes(x=x, y=y)) +
  scale_fill_viridis_c(option="turbo") +
  theme_void() + theme(legend.position = "none")
