# Libraries and auxiliry functions

library(INLA)
library(rgeos)
library(ggplot2)

rmse <- function(Sim,Pred){
  return(sqrt(mean((Sim-Pred)**2)))
}

mesh.dual <- function(mesh){
  if (mesh$manifold=='R2'){
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0)
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1],
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2,
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

## Mesh-grid for simulation

Simulation <- function(globalseed, range0, sigma0, nsamplesInd, nsamplesDep, rscale, cov.formula="0.25*(x-mean(x))**2 + 0.5*(y-mean(y))**2", beta.vec=c(-1,1)){
  if(!missing(globalseed)){set.seed(globalseed)}
  limlattice <- c(0,1)
  lengthlattice <- 50
  xlattice0 <- seq(limlattice[1], limlattice[2], length.out=lengthlattice)
  ylattice0 <- seq(limlattice[1], limlattice[2], length.out=lengthlattice)
  lattice0 <-  inla.mesh.lattice(xlattice0, ylattice0)
  meshSim <-  inla.mesh.create(lattice=lattice0, refine=list(max.edge=0.08*abs(diff(limlattice))))
  lattice <- as.matrix(expand.grid(x=xlattice0,y=ylattice0))
  
  # sigma0 <- 0.6; range0 <- 0.2
  kappa0 <- sqrt(8)/range0; tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)
  spde <- inla.spde2.matern(meshSim, B.tau=cbind(log(tau0),-1,+1), B.kappa=cbind(log(kappa0),0,-1),
                            theta.prior.mean=c(0,0), theta.prior.prec=c(0.1,0.1))
  Qu <-  inla.spde.precision(spde, theta=c(0, 0))
  if(missing(globalseed)){
    u <- as.vector(inla.qsample(n=1, Q=Qu))
  } else{
    u <- as.vector(inla.qsample(n=1, Q=Qu, seed=globalseed))
  }
  
  grid <- as.matrix(expand.grid(
    seq(limlattice[1], limlattice[2], length.out=150),
    seq(limlattice[1], limlattice[2], length.out=150)))
  
  u.eff <- as.vector(inla.spde.make.A(mesh=meshSim, loc=grid)%*%u)
  
  cov.fun <- function(x,y,s){
    z <- eval(parse(text=s))
  }
  
  cov <- cov.fun(x=grid[,1], y=grid[,2], s=cov.formula)
  
  beta <- beta.vec
  predictor <- beta[1] + beta[2]*cov + u.eff
  
  var.gamma=0.001
  ysim <- rgamma(length(predictor), exp(predictor)**2/var.gamma, exp(predictor)/var.gamma)
  
  DFSim <- data.frame(x=grid[,1], y=grid[,2], ysim=ysim, cov=cov, u.effect=u.eff)
  
  ## Sample simulation
  
  DFSimSampleInd <- DFSim[sample(1:nrow(DFSim), size=nsamplesInd),]
  DFSimSampleDep <- DFSim[sample(1:nrow(DFSim), size=nsamplesDep, prob=exp(rscale*predictor)/sum(exp(rscale*predictor)/nrow(DFSim)**2)),]
  return(list(DFSim=DFSim, DFSampleInd=DFSimSampleInd, DFSampleDep=DFSimSampleDep, meshSim=meshSim))
}


# Sim <- Simulation(globalseed=1, range0=0.3, sigma0=1, nsamples=100, rscale=1.5, cov.formula="0.5*x + 1.2*y**2")
