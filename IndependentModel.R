# Independent model

IndependentModel <- function(dataSample, dataSim, mesh, priormaterntype, priormatern, feedbacktype, prior.fixed, prior.hyper.design, prior.mode.location){
  if(priormaterntype=="pcmatern"){
    prior.range <- priormatern$range
    prior.sigma <- priormatern$sigma
    spde <- inla.spde2.pcmatern(mesh=mesh, prior.range=prior.range, prior.sigma=prior.sigma)
  } else if(priormaterntype=="matern"){
    range0 <- priormatern$range0
    sigma0 <- priormatern$sigma0
    theta.prior.mean <- priormatern$theta.prior.mean
    theta.prior.prec <- priormatern$theta.prior.prec
    kappa0 <- sqrt(8)/range0; tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
    spde <- inla.spde2.matern(mesh, B.tau=cbind(log(tau0),-1,+1), 
                              B.kappa=cbind(log(kappa0),0,-1),
                              theta.prior.mean=theta.prior.mean, theta.prior.prec=theta.prior.prec)
  }
  
  spde.index <- inla.spde.make.index(name="spatial", n.spde=spde$n.spde)

  A.inf <- inla.spde.make.A(mesh=mesh, loc=as.matrix(dataSample[,1:2]))
  A.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(dataSim[,1:2]))

  Inf.stack.effects <- list(list(spatial=spde.index$spatial),
                            list(Intercept = rep(1, nrow(dataSample)),
                                 Cov = dataSample[,4]))
  
  Pred.stack.effects <- list(list(spatial=spde.index$spatial),
                            list(Intercept = rep(1, nrow(dataSim)),
                                 Cov = dataSim[,4]))
  
  Inf.stack <- inla.stack(data=list(y=dataSample$ysim),
                          A=list(A.inf,1),
                          effects=Inf.stack.effects,
                          tag="Inference")
  
  Pred.stack <- inla.stack(data=list(y=NA),
                          A=list(A.pred,1),
                          effects=Pred.stack.effects,
                          tag="Prediction")
  
  Total.stack <- inla.stack(Inf.stack, Pred.stack)
  
  formula = y ~ -1 + Intercept + Cov + f(spatial, model=spde)
  
  control.fixed <- list(mean=list(Interceptor=prior.fixed$meanInterceptor, Cov=prior.fixed$meanCov),
                     prec=list(Interceptor=prior.fixed$precInterceptor, Cov=prior.fixed$precCov))
  
  if(feedbacktype=="moments"){
    if(missing(prior.mode.location)){
      IndependentModel <- inla(formula = formula, family="gamma",
                               data = inla.stack.data(Total.stack),
                               control.inla = list(strategy="auto", int.strategy="auto"),
                               control.predictor = list(A=inla.stack.A(Total.stack), compute=TRUE, link=1),
                               control.fixed = control.fixed,
                               control.compute = list(config=FALSE, dic=TRUE, waic=TRUE),
                               inla.mode="compact",
                               verbose=FALSE)
    } else if(!missing(prior.mode.location)){
      control.mode <- list(theta=prior.mode.location, restart=TRUE)
      IndependentModel <- inla(formula = formula, family="gamma",
                               data = inla.stack.data(Total.stack),
                               control.inla = list(strategy="auto", int.strategy="auto"),
                               control.predictor = list(A=inla.stack.A(Total.stack), compute=TRUE, link=1),
                               control.fixed = control.fixed,
                               control.mode = control.mode,
                               control.compute = list(config=FALSE, dic=TRUE, waic=TRUE),
                               inla.mode="compact",
                               verbose=FALSE)
    }
    
  } else if(feedbacktype=="fullhyperpar"){
    IndependentModel <- inla(formula = formula, family="gamma",
                             data = inla.stack.data(Total.stack),
                             control.inla = list(strategy="auto", int.strategy="user", int.design = prior.hyper.design),
                             control.predictor = list(A=inla.stack.A(Total.stack), compute=TRUE, link=1),
                             control.fixed = control.fixed,
                             control.mode = control.mode,
                             control.compute = list(config=FALSE, dic=TRUE, waic=TRUE),
                             inla.mode="compact",
                             verbose=FALSE)
  }
  
  index.pred <- inla.stack.index(stack=Total.stack, tag="Prediction")$data
  DFPredIM <- data.frame(x=dataSim$x, y=dataSim$y)
  DFPredIM$MeanDataPred <- IndependentModel$summary.fitted.values[index.pred, "mean"]
  DFPredIM$MedianDataPred <- IndependentModel$summary.fitted.values[index.pred, "0.5quant"]
  DFPredIM$SdDataPred <- IndependentModel$summary.fitted.values[index.pred, "sd"]
  
  return(list(IndepedentModel=IndependentModel, index.pred=index.pred, RMSEmean=rmse(Sim=dataSim$ysim, Pred=DFPredIM$MeanDataPred), RMSEmedian=rmse(Sim=dataSim$ysim, Pred=DFPredIM$MedianDataPred)))
}

# limlattice <- c(0,1)
# lengthlattice <- 25
# xlattice0 <- seq(limlattice[1], limlattice[2], length.out=lengthlattice)
# ylattice0 <- seq(limlattice[1], limlattice[2], length.out=lengthlattice)
# lattice0 <-  inla.mesh.lattice(xlattice0, ylattice0)
# meshInf <-  inla.mesh.create(lattice=lattice0, refine=list(max.edge=0.08*abs(diff(limlattice))))
# 
# 
# 
# priormatern <- list(range0=0.2, sigma0=1, theta.prior.mean=c(0,0), theta.prior.prec=c(0.01,0.01))
# prior.fixed <- list(meanInterceptor=0, meanCov=0, precInterceptor=0.001, precCov=0.001)
# IndependentModelMatern <- IndependentModel(dataSample=Sim$DFSampleInd, dataSim=Sim$DFSim, mesh=meshInf, 
#                                    priormaterntype="matern", priormatern=priormatern,
#                                    feedbacktype="moments", prior.fixed=prior.fixed,
#                                    prior.mode.location = c(NA,NA,NA))
# 
# priormatern <- list(range=c(0.2,0.5), sigma=c(1,0.5))
# IndependentModelPCMatern <- IndependentModel(dataSample=Sim$DFSampleInd, dataSim=Sim$DFSim, mesh=meshInf, 
#                                      priormaterntype="pcmatern", priormatern=priormatern,
#                                      feedbacktype="moments", prior.fixed=prior.fixed,
#                                      prior.mode.location = c(NA,NA,NA))
