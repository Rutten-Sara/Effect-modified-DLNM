
################################################################################
# Function to calculate RRR (changes in RR)
# - If interaction modification: value of the effect modifier and reference value are needed
# - Otherwise: calculate the RRR for an increase of 1 unit in the effect modifier
################################################################################
#   - model: result from DLNM_Laplace_cov
#   - at_x: value of exposure values to predict in
#   - cen: centering value (for RR)
#   - L: lag at which to calculate RR
#   - at_inter: value of effect modifier (if interaction modifying type)
#   - ref_inter: reference value of effect modifier (if interaction modifying type)
################################################################################

InterpretRR <- function(model,at_x, cen, L, at_inter = NULL, ref_inter = NULL){
  
  
  if(is.null(at_inter)& is.null(ref_inter)){

  crossbasis <- model$crossbasis
  ind.cb <- model$ind.cb
  ind.cb_inter <- model$ind.cb_inter
  xi_mode_cb <- model$xi_mode_cb
  xi_mode_cb_inter <- model$xi_mode_cb_inter
  sigma_xi <- model$Sigma
  predlag <- L # lag 0:L
  nrow <- length(at_x) # number of predictions
  at_x <- rep(at_x, length(predlag)) # Replicate values for every lag
 
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
  

    # Overall RR
    Xpredall <- 0
    for (j in 1:length(predlag)) {
      if(predlag[j]==round(predlag[j])){
        ind_all <- seq(nrow) + nrow * (j - 1) 
        Xpredall <- Xpredall + Xpred[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
      }
    }

      log_RR_change = Xpredall%*%xi_mode_cb_inter
      
      
      sd_all <- sqrt(Matrix::diag(Xpredall %*% sigma_xi[ind.cb_inter, ind.cb_inter] %*% t(Xpredall)))
      quantiles_all <- mapply(function(mean_val, sd_val) {
        qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
      }, log_RR_change, sd_all)
      Qlower_all = quantiles_all[1,]
      Qupper_all = quantiles_all[2,]
      


  }else{
    log_RR_change = Qlower_all = Qupper_all = NULL
    ind_k = 1
    
    crossbasis <- model$crossbasis
    ind.cb <- model$ind.cb
    ind.cb_inter <- model$ind.cb_inter
    xi_mode_cb <- model$xi_mode_cb
    xi_mode_cb_inter <- model$xi_mode_cb_inter
    sigma_xi <- model$Sigma
    predlag <- L # lag 0:L
    nrow <- length(at_x) # number of predictions
    at_x <- rep(at_x, length(predlag)) # Replicate values for every lag
    
    
    
    
    
    # Design matrix for centering z value
    if(!is.null(model$knots_inter)){
      at_inter_ref = as.matrix(ps(ref_inter, knots = model$knots_inter))
    }else{
      at_inter_ref = ref_inter
    }
    
    at_inter_ref <- matrix(rep(at_inter_ref, each = nrow*length(predlag)), nrow = nrow*length(predlag), byrow = FALSE) 
    
    basisvar <- do.call("onebasis", c(list(x = at_x), attr(crossbasis,"argvar")))
    basislag <- do.call("onebasis", c(list(x = predlag), attr(crossbasis,"arglag")))
    
    # basis variables for center
    basiscen <- do.call("onebasis", c(list(x = cen), attr(crossbasis,"argvar")))
    # center data matrix
    Xpred_cen <- scale(basisvar,center = basiscen, scale = F)
    
    # Prediction matrix
    Xpred_global <- matrix(0, nrow=length(at_x), ncol = ncol(crossbasis))
    
    for (l in seq(length = length(predlag))){
      for (v in seq(length = ncol(basisvar))){
        for (k in 1:ncol(basislag)){
          Xpred_global[((l-1)*nrow+1):(nrow*l),(ncol(basislag)*(v-1)+k)] = Xpred_cen[((l-1)*nrow+1):(nrow*l),v]*basislag[l,k] # prediction for every variable at lag = l
        } 
      }
    }
    
    
    Mat_crossbasis_spat <- Matrix::Matrix(Xpred_global, sparse = TRUE) # cross basis
    
    interaction_list_ref <- lapply(seq_len(ncol(Mat_crossbasis_spat)), function(j) {
      at_inter_ref * Mat_crossbasis_spat[, j]
    })
    
    interaction_matrix_ref <- do.call(cbind, interaction_list_ref)
    
    # Global design marix
    Xpred_ref <- Matrix::Matrix(cbind(Xpred_global,interaction_matrix_ref), sparse = TRUE)
    
    
    
    
    
    for (k in at_inter){
      
      if(!is.null(model$knots_inter)){
        at_k = as.matrix(ps(k, knots = model$knots_inter))
      }else{
        at_k = k
      }
      
    at_k <- matrix(rep(at_k, each = nrow*length(predlag)), nrow = nrow*length(predlag), byrow = FALSE) 

    interaction_list <- lapply(seq_len(ncol(Mat_crossbasis_spat)), function(j) {
      at_k * Mat_crossbasis_spat[, j]
    })
    
    interaction_matrix <- do.call(cbind, interaction_list)
    
    # Global design marix
    Xpred <- Matrix::Matrix(cbind(Xpred_global,interaction_matrix), sparse = TRUE) - Xpred_ref
    
    
    # Overall RR
    Xpredall <- 0
    for (j in 1:length(predlag)) {
      if(predlag[j]==round(predlag[j])){
        ind_all <- seq(nrow) + nrow * (j - 1) 
        Xpredall <- Xpredall + Xpred[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
      }
    }

    log_RR_change[ind_k] = as.numeric(Xpredall%*%c(xi_mode_cb,xi_mode_cb_inter))
    
    sd_all <- sqrt(Matrix::diag(Xpredall %*% sigma_xi[c(ind.cb, ind.cb_inter), c(ind.cb, ind.cb_inter)] %*% t(Xpredall)))
    quantiles_all <- mapply(function(mean_val, sd_val) {
      qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
    }, log_RR_change[ind_k], sd_all)
    Qlower_all[ind_k] = quantiles_all[1,]
    Qupper_all[ind_k] = quantiles_all[2,]
    
    
    
    ind_k = ind_k + 1
    
  }

  }
  
  
  return(list(log_RR_change = log_RR_change, log_lower = Qlower_all, log_upper = Qupper_all))
  
    
    
}

