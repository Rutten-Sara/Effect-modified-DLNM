
########################## Multiple time series ################################

rm(list = ls())
library(tidyverse) ; library(readxl); library(Rcpp); library(sf)
library(spdep); library(raster); library(exactextractr)
library(lubridate); library(mgcv)
library(dlnm)
library(tsModel); library(gnm)
library(R.utils);library(mixmeta); library(splines)

setwd("G:/My Drive/Onderzoek/DLNM/Spatial varying DLNM/Code/Covariate varying")

source('functions/DLNM_Laplace_covariate.R')
source('functions/DLNM_Laplace_covariate_NB.R')
source('functions/DLNM_Laplace_covariate_poisson.R')
source('functions/predRR_covariate.R')
source('functions/help_functions.R')
source('functions/BAM_predict.R')
source('functions/interpret_coefficients.R')



# Load data

load("data/daily_data.RData")
shapefile_bcn <- read_sf("data/shapefile_bcn.shp")
data$year <- year(data$date)


# Standardize the real temperature data
x <- data$temp
x <- (x-min(x))/diff(range(x))*10
data$x <- x

# Construct lagged exposures
group <- factor(paste(data$region, data$year, sep="-"))
Q <- Lag(x,0:8, group=group)

# Function to construct exposure-lag-response surface 
fflex <- function(x, coef) {
  as.numeric(outer(x,0:4,'^')%*%coef)
}

wdecay <- function(lag, denominator) exp(-lag/denominator)
ftemp <- function(x,lag, coef, denominator) 0.1 * (fflex(x, coef)-fflex(5,coef)) * 
  wdecay(lag, denominator)



# number of simulations
nsim <- 250

# nominal value
qn <- qnorm(0.975)



# True effect surface
trueeff <- function(fun, coef, denominator) {
  temp <- outer(seq(0,10,0.25),0:8,fun, coef, denominator)
  dimnames(temp) <- list(seq(0,10,0.25),paste("lag",0:8,sep=""))
  return(temp)
}


# Structured random effect
# Neighbourhood matrix
neig.map <- spdep::poly2nb(shapefile_bcn,row.names = shapefile_bcn$CBarri)
S <- length(unique(data$region))
Rn <- matrix(0, nrow = S, ncol = S)

for (s in 1:S) {
  # Diagonal elements (N_s)
  Rn[s, s] <- length(neig.map[[s]])
  
  # Off-diagonal elements (-1 for neighbors)
  for (u in neig.map[[s]]) {
    Rn[s, u] <- -1
  }
}
Rn <- Matrix::Matrix(Rn, sparse = TRUE)

# Random effect simulation
set.seed(1)
var_true_ind = 0.2
Q_spat_ind <- var_true_ind * solve(Matrix::Diagonal(n = S, 
                                                    x = 1 - 0.95) + 
                                     Matrix::Matrix(0.95*Rn, sparse = T)) 


eigen_spat_ind <- eigen(Q_spat_ind)
X <- matrix(rnorm(S*nsim),nsim)
spat_sim_ind <- eigen_spat_ind$vectors %*% diag(sqrt(eigen_spat_ind$values),S) %*% t(X)


# Simulate multiple spatially structured curves

#Coefficients
coefficients_all = lapply(1:5, function(i) matrix(NA, nrow = dim(shapefile_bcn)[1], ncol = nsim))

coef_base <- c(0.2118881,0.1406585,-0.0982663,0.0153671,-0.0006265) # base coefficients

coef_linear <- c(0.2,0.1,0.5,0.3,0.15) # Region-specific deviations

var_true = 0.4
z_sim = matrix(NA, nrow = dim(shapefile_bcn)[1], ncol = nsim) # Matrix with simulated covariate for each region
denominator_all = matrix(nrow = dim(shapefile_bcn)[1], ncol = nsim) # Matrix with simulated denominator (lag effect) for each region
for (k in 1:nsim){
  set.seed(10*k)
  cov_sim = rnorm(dim(shapefile_bcn)[1],0,var_true) # Simulate covariates
  
  z_sim[,k] = cov_sim - mean(cov_sim) # Center covariates
  for (i in 1:5){
    coefficients_all[[i]][,k] <- coef_base[i]*(1+coef_linear[i]*z_sim[,k]) # Region-specific coefficients
  } 
  denominator_all[,k] = rep(2, length(z_sim[,k])) # Common denominator
}

# Simulate population size (large regions)
set.seed(123)
offset_true = ceiling(rlnorm(n = 73, meanlog = log(6e6), sdlog = 0.9))



# Specifications DLNM
at_x = seq(0,10,0.25) # Values for prediction
L <- 8 # maximum lag
vx_pen <- 7 # number of basis for exposure var (penalized)
vl_pen <- 8 # number of basis for lag var (penalized)

# Knots for unpenalized DLNM
knots_vx <- quantile(x, probs = c(0.1,0.9))
knots_vl <- logknots(0:L, nk = 2)



# Store results
cov_RR <- rmse_RR <- bias_RR <- list("Smooth Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Linear Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Binary Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Smooth Unpen LPS" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Linear Unpen LPS" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Smooth Meta" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Linear Meta" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Smooth BAM" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "Linear BAM" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "CTS (Linear)" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)),
                                     "common" = matrix(0, nrow = 73, ncol = length(at_x)*(L+1)))



cov_all <- rmse_all <- bias_all <- list("Smooth Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Linear Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Binary Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Smooth Unpen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Linear Unpen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Smooth Meta" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Linear Meta" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Smooth BAM" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "Linear BAM" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "CTS (Linear)" = matrix(0, nrow = 73, ncol = length(at_x)),
                                        "common" = matrix(0, nrow = 73, ncol = length(at_x)))

predall_matrix   <-  list("Smooth Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "Linear Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "Binary Pen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "Smooth Unpen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "Linear Unpen LPS" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "Smooth Meta" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "Linear Meta" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "Smooth BAM" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "Linear BAM" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "CTS (Linear)" = matrix(0, nrow = 73, ncol = length(at_x)),
                          "common" = matrix(0, nrow = 73, ncol = length(at_x)))


time  <- rep(0,11)
DIC_percentage <- rep(0,4)
AIC_percentage <- rep(0,2)


which_fails <- rep(0,11)


tot_sim = nsim



# Helper function for meta analysis
fit_mixmeta <- function(formula, S, methods = c("reml", "ml", "mm", "vc")) {
  for (m in methods) {
    fit <- try(
      mixmeta(formula, S, method = m, control = list(maxiter = 10000)),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) {
      return(fit)
    }
  }
  stop("All methods failed for this model")
}

for (i in 1:nsim){
  
  succes_flags <- T
  succes_meta <- T
  
  print(i)
  
  #################################
  # Simulate true data            #
  #################################
  
  cumeff <- NULL
  
  # Extract (region-specific) coefficients 
  # Calculate DLNM effect based on exposure and coefficients
  for (j in 1:length(x)){
    coef_ij <- as.numeric(cbind((coefficients_all[[1]])[as.numeric(data$region[j]),i],
                                (coefficients_all[[2]])[as.numeric(data$region[j]),i],
                                (coefficients_all[[3]])[as.numeric(data$region[j]),i],
                                (coefficients_all[[4]])[as.numeric(data$region[j]),i],
                                (coefficients_all[[5]])[as.numeric(data$region[j]),i]))
    cumeff[j] = sum(do.call("ftemp", list(Q[j,], 0:8, coef_ij, denominator_all[as.numeric(data$region[j]),i])))
    
  }
  
  # True effect in each region for the variables seq(0,10,0.25)
  trueeff_sim = list()
  for (k in 1:73){
    coef_ik <- as.numeric(cbind((coefficients_all[[1]])[k,i],
                                (coefficients_all[[2]])[k,i],
                                (coefficients_all[[3]])[k,i],
                                (coefficients_all[[4]])[k,i],
                                (coefficients_all[[5]])[k,i]))
    trueeff_sim[[k]] = trueeff("ftemp", coef_ik, denominator_all[k,i])
  }
  
  
  
  
  trueeff_sim_mat <- do.call(rbind, lapply(trueeff_sim, function(x) as.numeric(x))) # true lag-specific RR
  trueeff_allsim_mat <- do.call(rbind, lapply(trueeff_sim, function(x) apply(x,1,sum))) #true overall RR
  
  # Simulate
  set.seed(12805+i)
  
  unique_area = data.frame(area = unique(data$region), offset = offset_true, z = z_sim[,i]) # Region specific size, covariate
  unique_area$random_effect = spat_sim_ind[,i] # Region-specific random effect
  
  random_area = inner_join(data.frame(area = data$region), unique_area,
                           by = "area")$random_effect # Merge random effect with data
  
  offset = inner_join(data.frame(area = data$region), unique_area,
                      by = "area")$offset # Merge population size with data
  
  z = inner_join(data.frame(area = data$region), unique_area,
                 by = "area")$z # Merge covariate with data
  

  suppressWarnings(y_all <- rnbinom(length(x),size = 5, mu = offset*exp(-10.5-z+cumeff+random_area))) # Simulate response
  data$y = y_all
  
  
  

  
  ##################################
  # Prepare model estimation       #
  ##################################
  
  #Prepare crossbasis matrix 
  crossbasis_pen <- crossbasis(x,
                               argvar=list(fun="ps",df = vx_pen, intercept=F),
                               arglag=list(fun="ps",df = vl_pen, intercept=T),
                               group = group, lag = L)
  
  crossbasis_unpen <- crossbasis(x,
                                 argvar=list(fun="ns",knots = knots_vx, intercept=F),
                                 arglag=list(fun="ns",knots = knots_vl, intercept=T),
                                 group = group, lag = L)
  
  
  
  # Prepare data
  z_lin = z
  data$z = z
  unique_area = unique(data$region)
  dlist <- split(data, data$region)[unique_area]
  
  z_binary = matrix(ifelse(z > median(z),1,0), ncol = 1)
  
  model_laplace = list()
  z_list = list(z_lin, z_lin, z_binary, z_lin, z_lin)
  
  type_pen = list(2,2,2,NULL,NULL)
  type_inter = c("non_linear","linear","factor","non_linear","linear")
  
  crossbasis_list = c(rep(list(crossbasis_pen),3), rep(list(crossbasis_unpen),2))
  y_all[which(is.na(crossbasis_pen[, 1]))] <- 0
  
  
  
  
  #########################################
  #               Fit CTS                 #
  #########################################
  tryCatch({
    data$stratum = factor(data$region)
    interaction_list <- lapply(seq_len(ncol(crossbasis_unpen)), function(j) {
      z_lin * crossbasis_unpen[, j]
    })
    interaction_matrix <- do.call(cbind, interaction_list)
    
    mtime <- proc.time()
    mod_cts <- gnm(y_all ~ z_lin + crossbasis_unpen + interaction_matrix, 
                   eliminate=stratum, data=data, family=nb(-1))
    time_CTS =  (proc.time()-mtime)[3]
    
    ID_coef_cts = grepl("crossbasis_unpen|interaction_matrix", names(coef(mod_cts)))
    selected_coefs_cts <- coef(mod_cts)[ID_coef_cts]
    selected_var_cts  <- vcov(mod_cts)[ID_coef_cts, ID_coef_cts]
    
  },error = function(e){
    which_fails[10]<<- which_fails[10]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  
  #########################################
  #         Fit Two-stage approach        #
  #########################################
  
  
  # Save data
  ymat <- matrix(NA,length(dlist),12,dimnames=list(unique_area,paste("b",seq(12),sep="")))
  
  
  Sall <- vector("list",length(dlist))
  names(Sall) <- unique_area
  # First stage: loop over different areas
  tryCatch({
    mtime = proc.time()
    for(k in seq(dlist)) {
      
      
      # data
      sub <- dlist[[k]]
      
      group_meta <- factor(sub$year)
      
      # cross-basis
      cb <- crossbasis(sub$x,lag=L,argvar=attributes(crossbasis_unpen)$argvar,
                       arglag=attributes(crossbasis_unpen)$arglag, group = group_meta)
      
      
      
      # Region-specific models
      #Slag2 <-  diag((0:(vl_pen-1))^2)
      #cbPen <- cbPen(cb,addSlag=list(Slag2)) 
      
      mfirst <- gam(y ~ cb,family=nb(-1),sub)
      
      
      # Store results
      ymat[k,] <- coef(mfirst)[-1]
      Sall[[k]] <- vcov(mfirst)[-1,-1]
      
    }
    
    # Prediction basis
    bvar <- do.call("onebasis",c(list(x=at_x),attr(cb,"argvar")))
    
  },error = function(e){
    which_fails[7]<<- which_fails[7]+1
    which_fails[6]<<- which_fails[6]+1
    succes_flags <<- F
    succes_meta <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  # Second stage (meta-analysis)
  time_common = (proc.time() - mtime)[3]
  
  
  if (succes_meta){
    tryCatch({
      # Cumulative effect
      mvall <- fit_mixmeta(ymat~unique(z_lin),Sall)
      
      
      time_meta_lin = (proc.time() - mtime)[3]
      
      predall_meta <- blup.mixmeta(mvall,vcov=T)
      
      
    },error = function(e){
      which_fails[7]<<- which_fails[7]+1
      succes_flags <<- F
      cat("ERROR :",conditionMessage(e), "\n")})
    
    time_lin = (proc.time() - mtime)[3] - time_common
    
    tryCatch({
      # Cumulative effect
      z_smooth = ns(unique(z_lin), knots = quantile(z_lin, probs = c(0.1, 0.5,0.9)), intercept = F)
      mvall_smooth <- fit_mixmeta(ymat~z_smooth,Sall)
      
      
      time_meta_smooth = ((proc.time()-mtime)[3] -time_lin)
      
      
      
      predall_meta_smooth <- blup.mixmeta(mvall_smooth,vcov=T)
      
      
      
      
      
    },error = function(e){
      which_fails[6]<<- which_fails[6]+1
      succes_flags <<- F
      cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  
  print("META")  
  
  #########################################
  #           Fit BAM                    #
  #########################################
  
  
  Slag2 <-  diag((0:(vl_pen-1))^2)
  cbPen <- cbPen(crossbasis_pen,addSlag=list(Slag2)) 
  
  AIC_list <- rep(0,2)
  time_bam = NULL
  model_bam = list()
  
  for (c in 1:2){
    if(c==1){
      crossbasis_by = ps(z_lin, df = 5,
                         intercept = F)
      
      interaction_list <- lapply(seq_len(ncol(crossbasis_pen)), function(j) {
        crossbasis_by * crossbasis_pen[, j]
      })
      
      cbPen_inter <- list(Svar = as.matrix(cbPen$Svar%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)),
                          Slag = as.matrix(cbPen$Slag%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)),
                          Slag2 = as.matrix(cbPen$Slag2%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)),
                          Svar2 = as.matrix(Matrix::Diagonal(n = vx_pen, x = NULL)%x%Matrix::Diagonal(n = vl_pen, x = NULL)%x%attributes(crossbasis_by)$S))
      
    }else{
      interaction_list <- lapply(seq_len(ncol(crossbasis_pen)), function(j) {
        z_lin * crossbasis_pen[, j]
      })
      
      cbPen_inter <- cbPen
    }
    
    
    interaction_matrix <- do.call(cbind, interaction_list)
    
    ID_gam =   as.factor(data$region)
    z_gam = z_list[[c]]
    
    mtime = proc.time()
    
    tryCatch({
      model_bam[[c]] = tryCatch(
        {
          model_bam <- withTimeout(
            bam(y_all ~  z_gam + crossbasis_pen + interaction_matrix + s(ID_gam, bs="re")+
                  offset(log(offset)),family=nb(-1),
                paraPen=list(crossbasis_pen=cbPen, interaction_matrix = cbPen_inter)),
            timeout = 20*60,
            onTimeout = "error"
          )
        },
        error = function(err){
          model_bam[[c]] <- withTimeout(
            bam(y_all ~  z_gam + crossbasis_pen + interaction_matrix + s(ID_gam, bs="re")+
                  offset(log(offset)) ,family=nb(-1),
                paraPen=list(crossbasis_pen=cbPen, interaction_matrix = cbPen_inter), method = "REML"),
            timeout = 20*60,
            onTimeout = "error"
          )
        }
      )
    },error = function(e){
      which_fails[c+7]<<- which_fails[c+7]+1
      succes_flags <<- F
      cat("ERROR :",conditionMessage(e), "\n")})
    
    time_bam[c] = (proc.time()-mtime)[3]

    
    
  }
  
  print("BAM")
  

  
  #########################################
  #               Fit LPS                 #
  #########################################
  
  time_laplace = NULL
  for (c in 1:5){
    tryCatch({
      mtime <- proc.time()
      model_laplace[[c]] <- DLNM_Laplace_cov(y_all ~ 1,
                                             crossbasis = crossbasis_list[[c]],
                                             ID = as.factor(data$region) ,
                                             covar.ri = "Leroux",
                                             offset = offset,
                                             data = data, map = shapefile_bcn,
                                             z = z_list[[c]],
                                             type = type_inter[c],
                                             df_z = 1,
                                             pen_crossbasis = type_pen[[c]],
                                             DIC = T,
                                             family = "NB") 
      
    },error = function(e){
      which_fails[c]<<- which_fails[c]+1
      succes_flags <<- F
      cat("ERROR :",conditionMessage(e), "\n")})
    
    time_laplace[c] = (proc.time()-mtime)[3]
  }
  
  
  print("Laplace")
  
  
  #########################################
  #              Fit common LPS           #
  #########################################
  
  
  mtime <- proc.time()
  tryCatch({
    model_laplace_common <- DLNM_Laplace_cov(y_all ~ z_lin,
                                             crossbasis = crossbasis_pen,
                                             ID = as.factor(data$region) ,
                                             covar.ri = "Leroux",
                                             offset = offset,
                                             data = data, map = shapefile_bcn,
                                             type = 'none',
                                             pen_crossbasis = 2,
                                             DIC = T,
                                             family = "NB") 
    
  },error = function(e){
    which_fails[11]<<- which_fails[11]+1
    succes_flags <<- F
    cat("ERROR :",conditionMessage(e), "\n")})
  
  time_laplace[6] = (proc.time()-mtime)[3]

  
  print("Laplace common")
  
  
  
  #########################################
  #              Predict CTS              #
  #########################################
  
  if(succes_flags){
    predvar_matrix_cts = predvar_matrix_cts_lower = predvar_matrix_cts_upper = NULL
    predall_matrix_cts_lower = predall_matrix_cts_upper = predall_matrix_cts =NULL
    
   for(k in 1:dim(shapefile_bcn)[1]){
      
      ID_z = which(data$region == unique(data$region)[k])[1]
      at_z=z_lin[ID_z]
      
      
      pred_all_cts = predRR_bam(selected_coefs_cts, selected_var_cts,crossbasis_unpen, at_x, cen = 5, L = 8, 
                                at_inter = at_z, by = 1)
      predvar_matrix_cts <- rbind(predvar_matrix_cts, pred_all_cts$logpredX)
      predvar_matrix_cts_lower <- rbind(predvar_matrix_cts_lower, pred_all_cts$Qlower_logpredX)
      predvar_matrix_cts_upper <- rbind(predvar_matrix_cts_upper, pred_all_cts$Qupper_logpredX)
      
      predall_matrix[[10]][k,] <-  log(pred_all_cts$pred_all)
      predall_matrix_cts_lower <- rbind(predall_matrix_cts_lower, log(pred_all_cts$Qlower_all))
      predall_matrix_cts_upper <- rbind(predall_matrix_cts_upper, log(pred_all_cts$Qupper_all))
      
    }
    cov_RR[[10]] = cov_RR[[10]] + (trueeff_sim_mat >= predvar_matrix_cts_lower &
                                     trueeff_sim_mat <= predvar_matrix_cts_upper)
    rmse_RR[[10]] = rmse_RR[[10]] + (trueeff_sim_mat - predvar_matrix_cts)^2
    
    bias_RR[[10]] = bias_RR[[10]] + (trueeff_sim_mat - predvar_matrix_cts)
    
    
    cov_all[[10]]<- cov_all[[10]] + (trueeff_allsim_mat >= predall_matrix_cts_lower &
                                      trueeff_allsim_mat <= predall_matrix_cts_upper)
    rmse_all[[10]] <- rmse_all[[10]] + (trueeff_allsim_mat - predall_matrix[[10]])^2
    
    bias_all[[10]] <- bias_all[[10]] + (trueeff_allsim_mat - predall_matrix[[10]])
    
    
    time[10] = time[10] + time_CTS
    #########################################
    #          Predict Two-stage            #
    #########################################
    
    predvar_matrix_meta = predvar_matrix_meta_lower = predvar_matrix_meta_upper = NULL
    predall_matrix_meta_lower = predall_matrix_meta_upper = NULL
    
    for(k in 1:dim(shapefile_bcn)[1]){
      cpall_smooth <- crosspred(crossbasis_unpen,coef=predall_meta_smooth[[k]]$blup,vcov=predall_meta_smooth[[k]]$vcov,
                                model.link="log",by=1,at = at_x,cen=5)
      
      
      predall_matrix[[6]][k,] <- log(cpall_smooth$allRRfit)
      predall_matrix_meta_lower <- rbind(predall_matrix_meta_lower, log(cpall_smooth$allRRlow))
      predall_matrix_meta_upper <- rbind(predall_matrix_meta_upper, log(cpall_smooth$allRRhigh))
      
      
      
      
      predvar_matrix_meta <- rbind(predvar_matrix_meta, as.numeric(cpall_smooth$matfit))
      predvar_matrix_meta_lower <- rbind(predvar_matrix_meta_lower, log(as.numeric(cpall_smooth$matRRlow)))
      predvar_matrix_meta_upper <- rbind(predvar_matrix_meta_upper, log(as.numeric(cpall_smooth$matRRhigh)))
      
    }
    cov_RR[[6]] = cov_RR[[6]] + (trueeff_sim_mat >= predvar_matrix_meta_lower &
                                   trueeff_sim_mat <= predvar_matrix_meta_upper)
    rmse_RR[[6]] = rmse_RR[[6]] + (trueeff_sim_mat - predvar_matrix_meta)^2
    
    bias_RR[[6]] = bias_RR[[6]] + (trueeff_sim_mat - predvar_matrix_meta)
    
    
    cov_all[[6]] <- cov_all[[6]] + (trueeff_allsim_mat >= predall_matrix_meta_lower &
                                      trueeff_allsim_mat <= predall_matrix_meta_upper)
    rmse_all[[6]] <- rmse_all[[6]] + (trueeff_allsim_mat - predall_matrix[[6]])^2
    
    bias_all[[6]] <- bias_all[[6]] + (trueeff_allsim_mat - predall_matrix[[6]])
    
    time[6] = time[6] + time_meta_smooth
    
    
    predvar_matrix_meta = predvar_matrix_meta_lower = predvar_matrix_meta_upper = NULL
    predall_matrix_meta_lower = predall_matrix_meta_upper = NULL
    
    for(k in 1:dim(shapefile_bcn)[1]){
      cpall <- crosspred(crossbasis_unpen,coef=predall_meta[[k]]$blup,vcov=predall_meta[[k]]$vcov,
                         model.link="log",by=1,at = at_x,cen=5)
      
      
      predall_matrix[[7]][k,] <- log(cpall$allRRfit)
      predall_matrix_meta_lower <- rbind(predall_matrix_meta_lower, log(cpall$allRRlow))
      predall_matrix_meta_upper <- rbind(predall_matrix_meta_upper, log(cpall$allRRhigh))
      
      
      predvar_matrix_meta <- rbind(predvar_matrix_meta, as.numeric(cpall$matfit))
      predvar_matrix_meta_lower <- rbind(predvar_matrix_meta_lower, log(as.numeric(cpall$matRRlow)))
      predvar_matrix_meta_upper <- rbind(predvar_matrix_meta_upper, log(as.numeric(cpall$matRRhigh)))
      
    }
    cov_RR[[7]] = cov_RR[[7]] + (trueeff_sim_mat >= predvar_matrix_meta_lower &
                                   trueeff_sim_mat <= predvar_matrix_meta_upper)
    rmse_RR[[7]] = rmse_RR[[7]] + (trueeff_sim_mat - predvar_matrix_meta)^2
    
    bias_RR[[7]] = bias_RR[[7]] + (trueeff_sim_mat - predvar_matrix_meta)
    
    
    cov_all[[7]] <- cov_all[[7]] + (trueeff_allsim_mat >= predall_matrix_meta_lower &
                                      trueeff_allsim_mat <= predall_matrix_meta_upper)
    rmse_all[[7]] <- rmse_all[[7]] + (trueeff_allsim_mat - predall_matrix[[7]])^2
    
    bias_all[[7]] <- bias_all[[7]] + (trueeff_allsim_mat - predall_matrix[[7]])
    
    time[7] = time[7] + time_meta_lin    
    
    #########################################
    #              Predict BAM              #
    #########################################
    
    for (c in 1:2){
      AIC_list[c] = AIC(model_bam[[c]])
      
      
      ID_coef = grepl("crossbasis_pen|interaction_matrix", names(coef(model_bam[[c]])))
      selected_coefs <- coef(model_bam[[c]])[ID_coef]
      selected_var  <- vcov(model_bam[[c]])[ID_coef, ID_coef]
      
      
      predvar_matrix_bam = predvar_matrix_bam_lower = predvar_matrix_bam_upper = NULL
      predall_matrix_bam_lower = predall_matrix_bam_upper = NULL
      
      for(k in 1:dim(shapefile_bcn)[1]){
        
        ID_z = which(data$region == unique(data$region)[k])[1]
        
        if(is.vector(z_list[[c]])){
          at_z=z_list[[c]][ID_z]
        }else{
          at_z = z_list[[c]][ID_z,]
        }
        
        if(c==1){
          at_z = ps(at_z, knots = attributes(crossbasis_by)$knots)
        }
        
        pred_all_bam = predRR_bam(selected_coefs, selected_var,crossbasis_pen, at_x, cen = 5, L = 8, 
                                  at_inter = at_z, by = 1)
        predvar_matrix_bam <- rbind(predvar_matrix_bam, pred_all_bam$logpredX)
        predvar_matrix_bam_lower <- rbind(predvar_matrix_bam_lower, pred_all_bam$Qlower_logpredX)
        predvar_matrix_bam_upper <- rbind(predvar_matrix_bam_upper, pred_all_bam$Qupper_logpredX)
        
        predall_matrix[[c+7]][k,] <- log(pred_all_bam$pred_all)
        predall_matrix_bam_lower <- rbind(predall_matrix_bam_lower, log(pred_all_bam$Qlower_all))
        predall_matrix_bam_upper <- rbind(predall_matrix_bam_upper, log(pred_all_bam$Qupper_all))
        
      }
      cov_RR[[c+7]] = cov_RR[[c+7]] + (trueeff_sim_mat >= predvar_matrix_bam_lower &
                                         trueeff_sim_mat <= predvar_matrix_bam_upper)
      rmse_RR[[c+7]] = rmse_RR[[c+7]] + (trueeff_sim_mat - predvar_matrix_bam)^2
      
      bias_RR[[c+7]] = bias_RR[[c+7]] + (trueeff_sim_mat - predvar_matrix_bam)
      
      
      cov_all[[c+7]]<- cov_all[[c+7]] + (trueeff_allsim_mat >= predall_matrix_bam_lower &
                                           trueeff_allsim_mat <= predall_matrix_bam_upper)
      rmse_all[[c+7]] <- rmse_all[[c+7]] + (trueeff_allsim_mat - predall_matrix[[c+7]])^2
      
      bias_all[[c+7]] <- bias_all[[c+7]] + (trueeff_allsim_mat - predall_matrix[[c+7]])
      
      time[c+7] = time[c+7] + time_bam[c]
    }
    
    
    ind_min_AIC = which.min(AIC_list)
    AIC_percentage[ind_min_AIC] = AIC_percentage[ind_min_AIC]+1
    
    
    #########################################
    #              Predict LPS              #
    #########################################
    
    for (c in 1:5){
      predvar_matrix = predvar_matrix_lower = predvar_matrix_upper = NULL
      predall_matrix_lower = predall_matrix_upper = NULL
      
 
      for(k in 1:dim(shapefile_bcn)[1]){
        
        ID_z = which(data$region == unique(data$region)[k])[1]
        
        if(is.vector(z_list[[c]])){
          at_z=z_list[[c]][ID_z]
        }else{
          at_z = z_list[[c]][ID_z,]
        }
        
        pred_all = predRR(model_laplace[[c]], at_x, cen = 5, L = 8, 
                          at_inter = at_z, by = 1)
        predvar_matrix <- rbind(predvar_matrix, pred_all$logpredX)
        predvar_matrix_lower <- rbind(predvar_matrix_lower, pred_all$Qlower_logpredX)
        predvar_matrix_upper <- rbind(predvar_matrix_upper, pred_all$Qupper_logpredX)
        
        predall_matrix[[c]][k,] <- log(pred_all$pred_all)
        predall_matrix_lower <- rbind(predall_matrix_lower, log(pred_all$Qlower_all))
        predall_matrix_upper <- rbind(predall_matrix_upper, log(pred_all$Qupper_all))
        
      }
      
      
      cov_RR[[c]] = cov_RR[[c]] + (trueeff_sim_mat >= predvar_matrix_lower &
                                     trueeff_sim_mat <= predvar_matrix_upper)
      rmse_RR[[c]] = rmse_RR[[c]] + (trueeff_sim_mat - predvar_matrix)^2
      
      bias_RR[[c]] = bias_RR[[c]] + (trueeff_sim_mat - predvar_matrix)
      
      
      cov_all[[c]]<- cov_all[[c]] + (trueeff_allsim_mat >= predall_matrix_lower &
                                       trueeff_allsim_mat <= predall_matrix_upper)
      rmse_all[[c]] <- rmse_all[[c]] + (trueeff_allsim_mat - predall_matrix[[c]])^2
      
      bias_all[[c]] <- bias_all[[c]] + (trueeff_allsim_mat - predall_matrix[[c]])
      
      
      
      time[c] = time[c] + time_laplace[c]
      
    }
    
    
    #########################################
    #        Predict common LPS             #
    #########################################
    
    
    pred_common_matrix = pred_var_common_matrix_lower = pred_var_common_matrix_upper = NULL
    predall_common_matrix_lower = predall_common_matrix_upper = NULL
    
    
    pred_all_common = predRR(model_laplace_common, at_x, cen = 5, L = 8, by = 1)
    predvar_common_matrix <- matrix(rep(pred_all_common$logpredX, 73), nrow = 73, byrow = T)
    predvar_common_matrix_lower <- matrix(rep(pred_all_common$Qlower_logpredX,  73), nrow = 73, byrow = T)
    predvar_common_matrix_upper <- matrix(rep(pred_all_common$Qupper_logpredX, 73), nrow = 73, byrow = T)
    
    predall_matrix[[11]] <- matrix(rep(log(pred_all_common$pred_all),  73), nrow = 73, byrow = T)
    predall_common_matrix_lower <- matrix(rep(log(pred_all_common$Qlower_all), 73), nrow = 73, byrow = T)
    predall_common_matrix_upper <- matrix(rep(log(pred_all_common$Qupper_all), 73), nrow = 73, byrow = T)
    
    
    cov_RR[[11]] = cov_RR[[11]] + (trueeff_sim_mat >= predvar_common_matrix_lower &
                                     trueeff_sim_mat <= predvar_common_matrix_upper)
    rmse_RR[[11]] = rmse_RR[[11]] + (trueeff_sim_mat - predvar_common_matrix)^2
    
    bias_RR[[11]] = bias_RR[[11]] + (trueeff_sim_mat - predvar_common_matrix)
    
    
    cov_all[[11]]<- cov_all[[11]] + (trueeff_allsim_mat >= predall_common_matrix_lower &
                                       trueeff_allsim_mat <= predall_common_matrix_upper)
    rmse_all[[11]] <- rmse_all[[11]] + (trueeff_allsim_mat - predall_matrix[[11]])^2
    
    bias_all[[11]] <- bias_all[[11]] + (trueeff_allsim_mat - predall_matrix[[11]])
    
    time[11] = time[11] + time_laplace[6]
    
    
    ind_min_DIC = which.min(c(model_laplace[[1]]$DIC, model_laplace[[2]]$DIC,
                              model_laplace[[3]]$DIC, model_laplace_common$DIC))
    DIC_percentage[ind_min_DIC] = DIC_percentage[ind_min_DIC]+1
    
    
  }else{
    tot_sim <- tot_sim - 1
  }
}


# Summary measures
unlist(lapply(cov_all, function(x) mean(x/tot_sim)))
unlist(lapply(cov_RR, function(x) mean(x/tot_sim)))
unlist(lapply(rmse_all, function(x) sqrt(mean(x/tot_sim))))
unlist(lapply(rmse_RR, function(x) sqrt(mean(x/tot_sim))))

which_fails/tot_sim
time/tot_sim
DIC_percentage/tot_sim
AIC_percentage/tot_sim


