################################################################################
# Function to calculate RR (for BAM)
################################################################################
#   - selected_coef: coefficients of cross-basis (including interaction)
#   - selected_var: variance-covariance matrix of cross-basis
#   - crossbasis: DLNM crossbasis
#   - at_x: vector of exposure values to predict in
#   - L: lag at which to calculate RR
#   - by: steps in lag vector
#   - at_inter: value of effect modifier (if effect modification)
#   - cen: centering value (for RR)
################################################################################

predRR_bam <-function(selected_coef, selected_var, crossbasis, at_x, L, by = 1, at_inter,
                      cen){
  predlag <- seq(0,L, by = by) # lag 0:L
  nrow <- length(at_x) # number of predictions
  at_x <- rep(at_x, length(predlag)) # Replicate values for every lag
  at_inter <- matrix(rep(at_inter, each = nrow*length(predlag)), nrow = nrow*length(predlag), byrow = FALSE) # Replicate values for every lag
  
  basisvar <- do.call("onebasis", c(list(x = at_x), attr(crossbasis,"argvar")))
  basislag <- do.call("onebasis", c(list(x = predlag), attr(crossbasis,"arglag")))
  
  # basis variables for center
  basiscen <- do.call("onebasis", c(list(x = cen), attr(crossbasis,"argvar")))
  # center data matrix
  Xpred_cen <- scale(basisvar,center = basiscen, scale = F)
  
  # Prediction matrix
  Xpred <- matrix(0, nrow=length(at_x), ncol = ncol(crossbasis))
  
  for (l in seq(length = length(predlag))){
    for (v in seq(length = ncol(basisvar))){
      for (k in 1:ncol(basislag)){
        Xpred[((l-1)*nrow+1):(nrow*l),(ncol(basislag)*(v-1)+k)] = Xpred_cen[((l-1)*nrow+1):(nrow*l),v]*basislag[l,k] # prediction for every variable at lag = l
      } 
    }
  }
  
  Mat_crossbasis_spat <- Matrix::Matrix(Xpred, sparse = TRUE) # cross basis
  
  interaction_list <- lapply(seq_len(ncol(Mat_crossbasis_spat)), function(j) {
    at_inter * Mat_crossbasis_spat[, j]
  })
  
  interaction_matrix <- do.call(cbind, interaction_list)
  
  # Global design mtarix
  Xpred <- Matrix::Matrix(cbind(Xpred,interaction_matrix), sparse = TRUE)
  
  ##########################
  logpredX <- as.numeric(Xpred %*% selected_coef)
  sd_logpredX <- as.numeric(sqrt(Matrix::diag(Xpred %*% selected_var %*% t(Xpred))))
  
  quantiles <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, logpredX, sd_logpredX)
  Qlower_logpredX <- quantiles[1,]
  Qupper_logpredX <- quantiles[2,]
  ###############
  
  # Overall RR
  Xpredall <- 0
  for (j in 1:length(predlag)) {
    if(predlag[j]==round(predlag[j])){
      ind_all <- seq(nrow) + nrow * (j - 1) 
      Xpredall <- Xpredall + Xpred[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
    }
  }
  pred_all <- exp(as.numeric(Xpredall %*% selected_coef))
  sd_all <- sqrt(Matrix::diag(Xpredall %*% selected_var %*% t(Xpredall)))
  quantiles_all <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, log(pred_all), sd_all)
  Qlower_all = exp(quantiles_all[1,])
  Qupper_all = exp(quantiles_all[2,])
  
  output <- list("Xpred" = Xpred,
                 "logpredX" = logpredX,
                 "sd_logpredX" = sd_logpredX,
                 "Qlower_logpredX" = Qlower_logpredX,
                 "Qupper_logpredX" = Qupper_logpredX,
                 "Xpredall" = Xpredall,
                 "Qlower_all" = Qlower_all,
                 "Qupper_all" = Qupper_all,
                 "pred_all" = pred_all,
                 "sd_all" = sd_all)
  
  
}
