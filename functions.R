library(Matrix)
library(ramps)
library(foreach)
library(doMC)
library(MASS)
library(xgboost)
library(caret)
library(data.table)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(magrittr)
library(ggpubr)
library(cowplot)

make_fixed_layers_covs <- function(x_matern_range = 1.5,
                                   x_matern_scale = 2,
                                   z_matern_range = 1,
                                   z_matern_scale = 1,
                                   noise_matern_range = 0.5,
                                   noise_matern_scale = 1,
                                   Scal = 50^2) {
  require(ramps); require(nlme)
  x_matern = corRMatern(value = c(x_matern_range,x_matern_scale), form = ~ x + y)
  z_matern = corRMatern(value = c(z_matern_range,z_matern_scale), form = ~ x + y)
  noise_matern = corRMatern(value = c(noise_matern_range,noise_matern_scale), form = ~ x + y)
  coords <- expand.grid(x = seq(-1,1,len=sqrt(Scal)), y = seq(-1,1,len=sqrt(Scal))) %>% data.table()
  cor_mat_x <- corMatrix(Initialize(x_matern,coords))
  cor_mat_z <- corMatrix(Initialize(z_matern,coords))
  cor_mat_noise <- corMatrix(Initialize(noise_matern,coords))
  list(cor_mat_x = cor_mat_x, cor_mat_z = cor_mat_z, cor_mat_noise=cor_mat_noise)
}

make_fixed_layers <- function(covs,
                              beta_X=c(1, 2, 0.5, -0.5, 0, 0, 0), Z_level=1,
                              rmvncores = 16) {
  Scal <- dim(covs$cor_mat_x)[1]
  coords <- expand.grid(x = seq(-1,1,len=sqrt(Scal)), y = seq(-1,1,len=sqrt(Scal))) %>% data.table()
  require(ramps); require(mvnfast)
  X <- rmvn(length(beta_X), mu=rep(1,Scal), sigma = covs$cor_mat_x, ncores = rmvncores) %>% t() %>% scale()
  colnames(X) <- paste("X",1:(length(beta_X)),sep="")
  Z <- as.numeric(rmvn(1, mu=rep(1,Scal), sigma = covs$cor_mat_z, ncores = rmvncores) %>% t() %>% scale())
  Y <- as.numeric((X %*% beta_X + Z_level * Z)  %>% scale() )
  dball <- cbind(coords, X,"Z" = Z,"Y" = Y) 
  
  dball_melt <- melt(dball, id.vars = c("x","y"), measure.vars = colnames(dball)[-(1:2)])
  plot <- ggplot(data = dball_melt, aes(x=x,y=y,z=value)) + 
    geom_raster(aes(fill = value)) + 
    geom_contour(color = "white", binwidth = 0.5) + 
    facet_wrap(~ variable, ncol = 3) +
    scale_fill_gradient2() +
    labs(x="",y="")
  
  list(coords = coords, X=X, Z=Z, Y=Y, dball=dball, plot=plot)
}

loss <- function(y,yhat) {(y-yhat)^2}

risk <- function(losses, iw=NULL) {
  if (is.null(iw)) {r <- mean(losses)}
  else {r <- mean(losses * iw)}
  return(r)
}

make_p_cov <- function(p_matern_range=1,p_matern_scale=1,Scal = 50^2) {
  require(ramps); require(nlme); 
  p_matern = corRMatern(value = c(p_matern_range,p_matern_scale), form = ~ x + y)
  coords <- expand.grid(x = seq(-1,1,len=sqrt(Scal)), y = seq(-1,1,len=sqrt(Scal))) %>% data.table()
  cor_mat_p <- corMatrix(Initialize(p_matern,coords))
  list(cor_mat_p=cor_mat_p)
}


make_distributions <- function(pcov,rmvncores=16) {
  Scal <- dim(pcov$cor_mat_p)[1]
  q0 <- rep(1,Scal)
  q <- rep(1/Scal,Scal)
  coords <- expand.grid(x = seq(-1,1,len=sqrt(Scal)), y = seq(-1,1,len=sqrt(Scal))) %>% data.table()
  
  require(MASS); require(scales); require(mvnfast)
  p0 <- as.numeric(rmvn(1, mu=rep(1,Scal), sigma = pcov$cor_mat_p, ncores = rmvncores)  %>% 
                     scales::rescale(to = c(0,1)))
  p0 <- p0
  p <- p0/sum(p0)
  
  list(q=q, q0=q0, p=p, p0=p0)
}


sample_over_space <- function(Scal, S=100, distribution) {
  
  p0 <- distribution$p0
  p <- distribution$p
  S_p <- sample(rep(1:Scal, round(S/10 * p0)), size = S, replace = T)   
  
  list(S_p=S_p,p=p)
}


estimate_dist <- function(db_p,n) {
  ker <- kde2d(db_p$x, db_p$y, n = n+0.2*n,
               h = 1,
               lims = c(-1.2, 1.2, -1.2, 1.2))
  xfixed <- ker$x[ker$x >= -1 & ker$x <= 1]
  yfixed <- ker$y[ker$y >= -1 & ker$y <= 1]
  allphat0 <- ker$z[ker$y %in% yfixed, ker$x %in% xfixed] #* (4/nrow(ker$z)^2)
  allphat <- allphat0/sum(allphat0)
  xindex <- sapply(db_p$x, function(s) which.min(abs(xfixed-s)))
  yindex <- sapply(db_p$y, function(s) which.min(abs(yfixed-s)))
  phat <- numeric(length(db_p$x))
  for (i in seq_along(phat)) {phat[i] <- allphat[xindex[i],yindex[i]]}
  list(phat=phat,allphat=allphat)
}

run_cv <- function(db_p,k=10,b_repeats=5,method="loo",Scal, trn_ratio= 0.75) {
  
  db_p$phat <- estimate_dist(db_p,n=sqrt(Scal))$phat
  db_p$q <- 1/Scal
  formula <- "Y~X1+X2+X3+X4+X5+X6+X7"
  
  if (method == "loo") {
    cv <- matrix(NA,nrow(db_p),2)
    
    for (i in 1:nrow(db_p)) {
      db_p_trn <- db_p[-i,]
      db_p_val <- db_p[i,]
      m_p <- lm(formula, db_p_trn)
      yhat_p <- predict(m_p, db_p_val, allow.new.levels=TRUE)
      cv[i,1] <- risk(loss(db_p_val$Y,yhat_p))
      cv[i,2] <- risk(loss(db_p_val$Y,yhat_p), iw = db_p_val$q/db_p_val$phat)
    }
  }
  
  if (method == "repeated") {
    cv <- matrix(NA,b_repeats,2)
    
    for (i in 1:b_repeats) {
      intrain <- sample(1:nrow(db_p), trn_ratio * nrow(db_p))
      db_p_trn <- db_p[intrain,]
      db_p_val <- db_p[-intrain,]
      m_p <- lm(formula, db_p_trn)
      yhat_p <- predict(m_p, db_p_val, allow.new.levels=TRUE)
      cv[i,1] <- risk(loss(db_p_val$Y,yhat_p))
      cv[i,2] <- risk(loss(db_p_val$Y,yhat_p), iw = db_p_val$q/db_p_val$phat)
    }
  }
  
  if (method == "kfold") {
    
    db_p <- db_p[sample(1:nrow(db_p)),]
    cv <- matrix(NA,k,2)
    folds <- cut(seq(1,nrow(db_p)),breaks=k,labels=FALSE)
    
    for (i in 1:k) {
      inval <- which(folds==i,arr.ind=TRUE)
      db_p_trn <- db_p[-inval,]
      db_p_val <- db_p[inval,]
      m_p <- lm(formula, db_p_trn)
      yhat_p <- predict(m_p, db_p_val, allow.new.levels=TRUE)
      cv[i,1] <- risk(loss(db_p_val$Y,yhat_p))
      cv[i,2] <- risk(loss(db_p_val$Y,yhat_p), iw = db_p_val$q/db_p_val$phat)
    }
  }
  
  colnames(cv) <- c("Naive_CV","IW_CV")
  colMeans(cv)
}

make_realizations <- function(covs, n_realizations,fixed_layers,beta_X,Z_level,noiselevel,rmvncores=1) {
  require(mvnfast)
  noise_r <- rmvn(n_realizations, mu=rep(0,dim(covs$cor_mat_noise)[1]),
                  sigma = covs$cor_mat_noise, ncores = rmvncores) %>% t() %>% scale()
  
  realizations <- lapply(split(noise_r, rep(1:ncol(noise_r), each = nrow(noise_r))), 
                         function (e) {
                           data.table(fixed_layers$coords,
                                      fixed_layers$X + noiselevel * as.numeric(e), 
                                      "Z" = fixed_layers$Z + noiselevel * as.numeric(e),
                                      "Y" = as.numeric((fixed_layers$X + noiselevel * as.numeric(e)) %*% beta_X + 
                                                         Z_level * (fixed_layers$Z + noiselevel * as.numeric(e)))
                           )})
  return(realizations)
}

compare_risks <- function(samples_S = 50, samples_target_risk = 50, S_p,
                          covs,fixed_layers,
                          beta_X=c(1, 2, 0.5, -0.5, 0, 0, 0),
                          noiselevel=1, Z_level=1,
                          k=10,b_repeats=5,trn_ratio=0.75,
                          rmvncores = 1, cores_samples = 10) {
  
  realizations_S <- make_realizations(covs = covs, fixed_layers = fixed_layers, n_realizations = samples_S,
                                      beta_X = beta_X, Z_level = Z_level, noiselevel = noiselevel,rmvncores = rmvncores)
  
  registerDoMC(cores = cores_samples)
  results <- 
    foreach(r = 1:samples_S, .combine = "rbind") %dopar% {
      
      db_all_r <- realizations_S[[r]]
      db_p <- db_all_r[S_p,]
      h_S <- lm(Y~X1+X2+X3+X4+X5+X6+X7, db_p)
      
      realizations_target_risk <- make_realizations(covs = covs, fixed_layers = fixed_layers,
                                                    n_realizations = samples_target_risk,
                                                    beta_X = beta_X, Z_level = Z_level, noiselevel = noiselevel)
      target_risk <- mean(
        sapply(realizations_target_risk,
               function(r) {mean((predict(h_S, r, allow.new.levels=TRUE) - r$Y)^2)}
        )
      )
      
      estimated_risks_LOO <- run_cv(db_p, method = "loo", Scal = dim(covs$cor_mat_x)[1])
      estimated_risks_RHO <- run_cv(db_p, method = "repeated", b_repeats = b_repeats, 
                                    Scal = dim(covs$cor_mat_x)[1],trn_ratio = trn_ratio)
      estimated_risks_KF <- run_cv(db_p, method = "kfold", k = k, Scal = dim(covs$cor_mat_x)[1])
      
      target_risk - c(estimated_risks_LOO, estimated_risks_RHO, estimated_risks_KF)
    } %>% data.table
  colnames(results) <- c("LOO","IWLOO","RHO","IWRHO","KF","IWKF")
  list(results=results, mean_results=colMeans(results))
}



gg_color_hue <- function(q) {
  hues = seq(15, 375, length = q + 1)
  hcl(h = hues, l = 65, c = 100)[1:q]
}


