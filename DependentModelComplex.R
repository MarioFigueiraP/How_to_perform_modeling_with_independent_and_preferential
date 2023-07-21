# Dependent model

DependentModelComplex <- function(dataSample, dataSim, mesh, priormaterntype, cov.formula=c("0.5*x + 1.2*y**2"), priormatern, feedbacktype, prior.fixed, prior.hyper.design, prior.mode.location, prior.scale.factor=c(0, 0.1)){
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
    spde <- inla.spde2.matern(meshInf, B.tau=cbind(log(tau0),-1,+1), 
                              B.kappa=cbind(log(kappa0),0,-1),
                              theta.prior.mean=theta.prior.mean, theta.prior.prec=theta.prior.prec)
  }
  
  spde.index <- inla.spde.make.index(name="spatial", n.spde=spde$n.spde)
  
  lim <- c(0,1)
  ldomain <- cbind(c(lim[1], lim[2], lim[2], lim[1], lim[1]), c(lim[1], lim[1], lim[2], lim[2], lim[1]))
  dmesh <- mesh.dual(mesh = mesh)
  domain.polys <- Polygons(list(Polygon(ldomain)), '0')
  domainSP <- SpatialPolygons(list(domain.polys))
  w <- sapply(1:length(dmesh), function(i) {
    if (gIntersects(dmesh[i, ], domainSP))
      return(gArea(gIntersection(dmesh[i, ], domainSP)))
    else return(0)
  })
  
  n <- nrow(dataSample)
  n.pp.ind <- nrow(dataSample[dataSample$kind=="ind",])
  n.pp.dep <- nrow(dataSample[dataSample$kind=="dep",])
  nv <- mesh$n
  # y.pp <- rep(0:1, c(nv, n))
  y.pp.ind <- rep(0:1, c(nv, n.pp.ind))
  y.pp.dep <- rep(0:1, c(nv, n.pp.dep))
  # e.pp <- c(w, rep(0, n))
  e.pp.ind <- c(w, rep(0, n.pp.ind))
  e.pp.dep <- c(w, rep(0, n.pp.dep))
  imat <- Diagonal(nv, rep(1, nv))
  # imat0 <- Diagonal(nv, rep(0, nv))
  lmat <- inla.spde.make.A(mesh=mesh, loc=as.matrix(dataSample[,1:2]))
  # lmat.pp.ind <- inla.spde.make.A(mesh=mesh, loc=as.matrix(dataSample[dataSample$kind=="ind",1:2]))
  lmat.pp.dep <- inla.spde.make.A(mesh=mesh, loc=as.matrix(dataSample[dataSample$kind=="dep",1:2]))
  # lmat.pp.ind <- lmat; lmat.pp.ind[dataSample$kind=="ind"] <- 0
  # lmat.pp.dep <- lmat; lmat.pp.dep[dataSample$kind=="dep"] <- 0
  
  # A.pp.lin.ind <- c(rep(1,nv), as.numeric(dataSample$kind=="ind"), 
  #                   rep(0,nv), rep(0, nrow(dataSample)))
  # A.pp.lin.dep <- c(rep(0,nv), rep(0, nrow(dataSample)),
  #                   rep(1,nv), as.numeric(dataSample$kind=="dep"))
  
  A.pp <- rbind(imat, lmat.pp.dep)
  # A.pp.ind <- rbind(imat, lmat.pp.ind)
  # A.pp.dep <- rbind(imat, lmat.pp.dep)
  
  A.pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(dataSim[,1:2]))
  
  cov.fun <- function(x,y,s){
    z <- eval(parse(text=s))
  }
  Cov.pp.mesh <- cov.fun(x=mesh$loc[,1], y=mesh$loc[,2], s=cov.formula)
  
  # Cov.stack.inf <- inla.stack(data = list(y = c(dataSample$cov)),
  #                             A = list(lmat, 1),
  #                             effects = list(list(spatial.cov = spde.index$spatial), 
  #                                            list(Interceptor.cov=rep(1, n))),
  #                             tag = 'Inference.cov')
  # 
  # A.cov.predmesh <- inla.spde.make.A(mesh=mesh, loc=mesh$loc[,1:2])
  # 
  # Cov.stack.pred <- inla.stack(data = list(y = rep(NA, nrow(A.cov.predmesh))),
  #                              A = list(A.cov.predmesh, 1),
  #                              effects = list(list(spatial.cov = spde.index$spatial), 
  #                                             list(Interceptor.cov=rep(1, nrow(A.cov.predmesh)))),
  #                              tag = 'Prediction.cov')
  # 
  # Total.stack.cov <- inla.stack(Cov.stack.inf, Cov.stack.pred)
  # 
  # Cov.formula <- y ~ -1 + (Interceptor.cov  +  f(spatial.cov, model = spde))
  # Cov.model <- inla(formula=Cov.formula, family = "gaussian", 
  #                   data = inla.stack.data(Total.stack.cov),
  #                   control.predictor = list(A = inla.stack.A(Total.stack.cov), compute=TRUE, link=1),
  #                   verbose=FALSE)
  # 
  # index.covpred <- inla.stack.index(stack=Total.stack.cov, tag="Prediction.cov")$data
  # Cov.pp.mesh <- Cov.model$summary.fitted.values[index.covpred, "mean"] 
  
  Inf.gs.stack <- inla.stack(data = list(y = cbind(dataSample$ysim, NA), e = rep(0, n)),
                             A = list(lmat, 1),
                             effects = list(list(spatial = spde.index$spatial), 
                                            list(Interceptor=rep(1, n), 
                                                 Cov=dataSample$cov)),
                             tag = 'Inference.gs')
  
  Inf.pp.stack.ind <- inla.stack(data = list(y = cbind(NA, y.pp.ind), e = e.pp.ind),
                             A = list(1),
                             effects = list(list(Interceptor.pp0=rep(1, nv+n.pp.ind))),
                             tag = 'Inference.pp.ind')
  
  Inf.pp.stack.dep <- inla.stack(data = list(y = cbind(NA, y.pp.dep), e = e.pp.dep),
                                 A = list(A.pp, 1),
                                 effects = list(list(spatial.pp=spde.index$spatial),
                                                list(
                                                  Interceptor.pp = rep(1, nv+n.pp.dep),
                                                  Cov.pp = c(Cov.pp.mesh, dataSample$cov[dataSample$kind=="dep"]))
                                 ),
                                 tag = 'Inference.pp.dep')
  
  
  Pred.gs.stack <- inla.stack(data=list(y=matrix(NA,nrow=nrow(dataSim), ncol=2)),
                              A = list(A.pred,1),
                              effects = list(list(spatial=spde.index$spatial),
                                             list(Interceptor=rep(1,nrow(dataSim)), 
                                                  Cov=dataSim$cov)),
                              tag = 'Prediction.gs')
  
  # Total.stack <- inla.stack(Inf.gs.stack, Pred.gs.stack, Inf.pp.stack)
  Total.stack <- inla.stack(Inf.gs.stack, Pred.gs.stack, Inf.pp.stack.ind, Inf.pp.stack.dep)

  control.fixed <- list(mean=list(Interceptor=prior.fixed$meanInterceptor, Cov=prior.fixed$meanCov, Interceptor.pp=prior.fixed$meanInterceptPP, Cov.pp=prior.fixed$meanCovPP),
                        prec=list(Interceptor=prior.fixed$precInterceptor, Cov=prior.fixed$precCov, Interceptor.pp=prior.fixed$precInterceptPP, Cov.pp=prior.fixed$precCovPP))
  
  # if(feedbacktype=="moments"){
    gaus.prior <- list(prior = 'gaussian', param = prior.scale.factor)
    formula <- y ~ -1 + (Interceptor.pp0 + Interceptor.pp + Cov.pp + f(spatial.pp, copy = 'spatial', fixed = FALSE, hyper = list(beta = gaus.prior))) + (Interceptor + Cov  +  f(spatial, model = spde))
    
    # if(missing(prior.mode.location)){
      DependentModel <- inla(formula = formula, family = c("gamma", "poisson"),
                             data = inla.stack.data(Total.stack),
                             E = inla.stack.data(Total.stack)$e,
                             control.inla = list(strategy="auto", int.strategy="auto"),
                             control.predictor = list(A=inla.stack.A(Total.stack), compute=TRUE, link=1),
                             # control.fixed = control.fixed,
                             control.compute = list(config=FALSE, dic=TRUE, waic=TRUE),
                             inla.mode="compact",
                             verbose=FALSE)
  #   } else if(!missing(prior.mode.location)){
  #     control.mode <- list(theta=prior.mode.location, restart=TRUE)
  #     DependentModel <- inla(formula = formula, family = c("gamma", "poisson"),
  #                            data = inla.stack.data(Total.stack),
  #                            E = inla.stack.data(Total.stack)$e,
  #                            control.inla = list(strategy="auto", int.strategy="auto"),
  #                            control.predictor = list(A=inla.stack.A(Total.stack), compute=TRUE, link=1),
  #                            control.fixed = control.fixed,
  #                            control.mode = control.mode,
  #                            control.compute = list(config=FALSE, dic=TRUE, waic=TRUE),
  #                            inla.mode="compact",
  #                            verbose=FALSE)
  #   }
  # } else if(feedbacktype=="fullhyperpar"){
  #   gaus.prior <- list(prior = 'gaussian', param = c(0, 0.001))
  #   formula <- y ~ -1 +  (Intercept.pp + Cov.pp + f(spatial.pp, copy = 'spatial', fixed = FALSE, hyper = list(beta = gaus.prior))) + (Interceptor + Cov  +  f(spatial, model = spde))
  #   DependentModel <- inla(formula = formula, family="gamma",
  #                          data = inla.stack.data(Total.stack),
  #                          E = inla.stack.data(Total.stack)$e,
  #                          control.inla = list(strategy="auto", int.strategy="user.expert", int.design = prior.hyper.design),
  #                          control.predictor = list(A=inla.stack.A(Total.stack), compute=TRUE, link=1),
  #                          control.fixed = control.fixed,
  #                          control.mode = control.mode,
  #                          control.compute = list(config=FALSE, dic=TRUE, waic=TRUE),
  #                          inla.mode="compact",
  #                          verbose=FALSE)
  # }
  
  index.pred <- inla.stack.index(stack=Total.stack, tag="Prediction.gs")$data
  DFPredDM <- data.frame(x=dataSim$x, y=dataSim$y)
  DFPredDM$MeanDataPred <- DependentModel$summary.fitted.values[index.pred, "mean"]
  DFPredDM$MedianDataPred <- DependentModel$summary.fitted.values[index.pred, "0.5quant"]
  DFPredDM$SdDataPred <- DependentModel$summary.fitted.values[index.pred, "sd"]
  
  return(list(DepedentModel=DependentModel, index.pred=index.pred, RMSEmean=rmse(Sim=dataSim$ysim, Pred=DFPredDM$MeanDataPred), RMSEmedian=rmse(Sim=dataSim$ysim, Pred=DFPredDM$MedianDataPred)))
}


# limlattice <- c(0,1)
# lengthlattice <- 25
# xlattice0 <- seq(limlattice[1], limlattice[2], length.out=lengthlattice)
# ylattice0 <- seq(limlattice[1], limlattice[2], length.out=lengthlattice)
# lattice0 <-  inla.mesh.lattice(xlattice0, ylattice0)
# meshInf <-  inla.mesh.create(lattice=lattice0, refine=list(max.edge=0.08*abs(diff(limlattice))))
# 
# 
# priormatern <- list(range0=0.2, sigma0=1, theta.prior.mean=c(0,0), theta.prior.prec=c(0.01,0.01))
# prior.fixed <- list(meanInterceptor=0, meanCov=0, precInterceptor=0.001, precCov=0.001)
# DependentModelMatern <- DependentModel(dataSample=Sim$DFSampleDep, dataSim=Sim$DFSim, mesh=meshInf,
#                                  priormaterntype="matern", priormatern=priormatern,
#                                  feedbacktype="moments", prior.fixed=prior.fixed,
#                                  cov.formula="0.25*(x-mean(x))**2 + 0.5*(y-mean(y))**2", prior.scale.factor=c(0, 0.001),
#                                  prior.mode.location = c(NA,NA,NA,NA))
# 
# priormatern <- list(range=c(0.2,0.5), sigma=c(1,0.5))
# DependentModelPCMatern <- DependentModel(dataSample=Sim$DFSampleDep, dataSim=Sim$DFSim, mesh=meshInf,
#                                  priormaterntype="pcmatern", priormatern=priormatern,
#                                  feedbacktype="moments", prior.fixed=prior.fixed,
#                                  cov.formula="0.25*(x-mean(x))**2 + 0.5*(y-mean(y))**2", prior.scale.factor=c(0, 0.001),
#                                  prior.mode.location = c(NA,NA,NA,NA))
