# This code is 

source("functions.R")

# grid size
n <- 64 

# generating covariance matrices for mean fields
set.seed(103); covs <- make_fixed_layers_covs(Scal = n^2,
                                              x_matern_range = 1.5,
                                              x_matern_scale = 2,
                                              z_matern_range = 1,
                                              z_matern_scale = 2,
                                              noise_matern_range = 2,
                                              noise_matern_scale = 2)
# generating mean fields from GRFs
set.seed(103); fixed_layers <- make_fixed_layers(covs = covs, beta_X=c(1, 2, 0.5, -0.5, 0, 0, 0), rmvncores = 16) #15
# simulation parameters
noiselevel <- 0.2; beta_X=c(1, 2, 0.5, -0.5, 0, 0, 0); Z_level = 2

# generating realization
e_r <- as.numeric(mvnfast::rmvn(1, mu=rep(0,n^2), sigma = covs$cor_mat_noise, ncores = 1) %>% t() %>% scale())
X_r <- fixed_layers$X + noiselevel * e_r
Z_r <- fixed_layers$Z + noiselevel * e_r
Y_r <- as.numeric((X_r %*% beta_X + Z_level * Z_r)  %>% scale() )
coords <- expand.grid(x = seq(-1,1,len=sqrt(n^2)), y = seq(-1,1,len=sqrt(n^2))) %>% data.table()
dball_r <- cbind(coords,as.data.table(X_r),"Z"=Z_r,"Y" = Y_r)
dball_melt <- melt(dball_r,id.vars = colnames(dball_r)[1:2], measure.vars = colnames(dball_r)[-(1:2)])

########################################################################################################
breaks = -3:3 
mycolours = c("black","grey20","grey30","grey50","grey80","grey90","white") 
(Figure1_gs <- ggplot(data = dball_melt, aes(x=x,y=y,z=value)) +
    geom_raster(aes(fill = value)) +
    geom_contour(color = "white", binwidth = 0.5) +
    facet_wrap(~ variable, ncol = 3) +
      scale_fill_gradientn(colours = mycolours, breaks=breaks) +
    labs(x="",y="") + theme(strip.text.x = element_text(size = 16)))
########################################################################################################

# generating covariance for source gid cells distribution
set.seed(100); pcov <- make_p_cov(Scal = n^2, p_matern_range = 2, p_matern_scale = 1)

# generating source gid cells distribution from GRF
set.seed(101); distribution <- make_distributions(pcov=pcov,rmvncores=16) #1

# sampling one realization of training set grid cells
set.seed(103); S_p <- sample(rep(1:n^2, round(500/10 * (distribution$p0)^1.5)), size = 300, replace = T)

########################################################################################################
g2_gs <- ggplot(data = dball_r, aes(x=x,y=y,z=Y)) +
  geom_raster(aes(fill = Y)) +
  geom_contour(color = "white", binwidth = 0.5) +
  scale_fill_gradientn(colours = mycolours[2:6]) +
  geom_jitter(data = dball_r[S_p,], aes(x=x,y=y), size = 2, shape = 21) +
  labs(x="",y="") +
  theme(legend.position="bottom", legend.text = element_text(size=10))
dball_dens <- cbind(rbind(dball_r, dball_r[S_p,]),
                    "is_train" = c(rep("All Domain   \n(Target sample)",nrow(dball_r)),
                                   rep("Training set \n(Source sample)",length(S_p))))
g3_gs <- ggplot(data = dball_dens, aes(x=Y, linetype = is_train)) +
  geom_density(alpha=0.3, adjust = 2.75) +
  theme_classic() + xlim(-3,2.5) +
  theme(legend.position="bottom", legend.text = element_text(size=12),
        legend.title=element_blank())
(Figure2_gs <- grid.arrange(g2_gs,g3_gs,ncol=2))
########################################################################################################



# training set sizes
Nsizes <- c(100,300,1000)

# sampling training set grid cells from source gid cells distribution
# comatre cross-validation target risk estimates' bias.
# repeat for 3 Nsizes.
set.seed(85); S_p1 <- sample(rep(1:n^2, round(Nsizes[1] * (distribution$p0)^2)), size = Nsizes[1], replace = T)
registerDoSEQ()
example1_risks_bias1 <- compare_risks(samples_S = 100, samples_target_risk = 100, S_p=S_p1,
                                      covs=covs,fixed_layers=fixed_layers,
                                      beta_X=c(1, 2, 0.5, -0.5, 0, 0, 0),
                                      noiselevel=0.05, Z_level=2,
                                      k = 5, b_repeats = 10, trn_ratio = 0.5, 
                                      rmvncores = 1, cores_samples = 150)
ex1_gg <- cbind(melt(example1_risks_bias1$results),  
                approach = factor(rep(c("LOO","RHO","KF"), each = 200)),
                iw = factor(rep(rep(c("Naive","IW"),3), each = 100)))

set.seed(84); S_p2 <- sample(rep(1:n^2, round(Nsizes[2] * (distribution$p0)^2)), size = Nsizes[2], replace = T)
registerDoSEQ()
example1_risks_bias2 <- compare_risks(samples_S = 100, samples_target_risk = 100, S_p=S_p2,
                                      covs=covs,fixed_layers=fixed_layers,
                                      beta_X=c(1, 2, 0.5, -0.5, 0, 0, 0),
                                      noiselevel=0.05, Z_level=2,
                                      k = 5, b_repeats = 10, trn_ratio = 0.5, 
                                      rmvncores = 1, cores_samples = 150)
ex2_gg <- cbind(melt(example1_risks_bias2$results), 
                approach = factor(rep(c("LOO","RHO","KF"), each = 200)),
                iw = factor(rep(rep(c("Naive","IW"),3), each = 100)))

set.seed(86); S_p3 <- sample(rep(1:n^2, round(Nsizes[3] * (distribution$p0)^2)), size = Nsizes[3], replace = T)
registerDoSEQ()
example1_risks_bias3 <- compare_risks(samples_S = 100, samples_target_risk = 100, S_p=S_p3,
                                      covs=covs,fixed_layers=fixed_layers,
                                      beta_X=c(1, 2, 0.5, -0.5, 0, 0, 0),
                                      noiselevel=0.05, Z_level=2,
                                      k = 5, b_repeats = 10, trn_ratio = 0.5, 
                                      rmvncores = 1, cores_samples = 150)
ex3_gg <- cbind(melt(example1_risks_bias3$results), 
                approach = factor(rep(c("LOO","RHO","KF"), each = 200)),
                iw = factor(rep(rep(c("Naive","IW"),3), each = 100)))


########################################################################################################
library(ggridges)
gs1 <- ggplot(ex1_gg, aes(y = variable, x = value)) +
  geom_vline(xintercept = 0, color = "gray20", size = 1, linetype = "dashed") +
  geom_density_ridges(scale = 1.5) + 
  xlim(-0.3,0.4) + labs(title = paste("N =",Nsizes[1], sep = " ")) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
gs2 <- ggplot(ex2_gg, aes(y = variable, x = value)) +
  geom_vline(xintercept = 0, color = "gray20", size = 1, linetype = "dashed") +
  geom_density_ridges(scale = 1.5) + 
  xlim(-0.3,0.4) + labs(title = paste("N =",Nsizes[2], sep = " ")) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
gs3 <- ggplot(ex3_gg, aes(y = variable, x = value)) +
  geom_vline(xintercept = 0, color = "gray20", size = 1, linetype = "dashed") +
  geom_density_ridges(scale = 1.5) + 
  xlim(-0.3,0.4) + labs(title = paste("N =",Nsizes[3], sep = " ")) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
Figure3_gs <- grid.arrange(gs1,gs2,gs3, ncol = 3)
########################################################################################################


