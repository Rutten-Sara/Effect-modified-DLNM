
################################################################################
#   - x: AN EXPOSURE VECTOR OR (ONLY FOR dir="back") A MATRIX OF LAGGED EXPOSURES
#   - z: THE VALUE OF THE EFFECT MODIFIER
#   - cases: THE CASES VECTOR OR (ONLY FOR dir="forw") THE MATRIX OF FUTURE CASES
#   - model: THE FITTED MODEL
#   - ID: ID OF OBSERVATION (THAT MATCHES ID IN OUTCOME MODEL)
#   - type: EITHER "an" OR "af" FOR ATTRIBUTABLE NUMBER OR FRACTION
#   - dir: EITHER "back" OR "forw" FOR BACKWARD OR FORWARD PERSPECTIVES
#   - tot: IF TRUE, THE TOTAL ATTRIBUTABLE RISK IS COMPUTED
#   - cen: THE REFERENCE VALUE USED AS COUNTERFACTUAL SCENARIO
#   - range: THE RANGE OF EXPOSURE. IF NULL, THE WHOLE RANGE IS USED
#   - sim: IF SIMULATION SAMPLES SHOULD BE RETURNED. ONLY FOR tot=TRUE
#   - nsim: NUMBER OF SIMULATION SAMPLES
#   - group: GROUP ID TO CALCULATE LAGGED EFFECTS
#   - ex.prob: CALULCATE EXCEEDANCE PROBABILITY OF 0?
################################################################################
attrdl_Laplace <- function(x,model,cases, z = NULL, ID = NULL,
                           type="af",dir="back",tot=TRUE,cen,range=NULL,sim=FALSE,nsim=5000, group=NULL,
                           ex.prob = NULL) {
  ################################################################################
  
  crossbasis <- model$crossbasis
  Sigma <- model$Sigma
  ind.cb <- model$ind.cb
  xi_mode_cb <- model$xi_mode_cb
  
  if(!is.null(z)){
  if(is.null(group) & !is.null(ID)){
    group = ID
  }
    ind.cb_inter <- model$ind.cb_inter
    xi_mode_cb_inter <- model$xi_mode_cb_inter
    
    xi_mode = c(xi_mode_cb, xi_mode_cb_inter)
    sigma_xi = Sigma[c(ind.cb, ind.cb_inter),c(ind.cb, ind.cb_inter)]
    
  }else{
    xi_mode = xi_mode_cb
    sigma_xi = Sigma[ind.cb, ind.cb]
  }

  
  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))
  #
  # DEFINE CENTERING
  if(missing(cen) && is.null(cen <- attr(crossbasis,"argvar")$cen))
    stop("'cen' must be provided")
  if(!is.numeric(cen) && length(cen)>1L) stop("'cen' must be a numeric scalar")
  attributes(crossbasis)$argvar$cen <- NULL
  #  
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- cen
  #
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"

  lag <- attr(crossbasis,"lag")
  if(NCOL(x)==1L) {
    at <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2]),group=group) else 
      matrix(rep(x,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(at <- x)!=diff(lag)+1) 
      stop("dimension of 'x' not compatible with 'basis'")
  }
  #
  # NUMBER USED FOR THE CONTRIBUTION AT EACH TIME IN FORWARD TYPE
  #   - IF cases PROVIDED AS A MATRIX, TAKE THE ROW AVERAGE
  #   - IF PROVIDED AS A TIME SERIES, COMPUTE THE FORWARD MOVING AVERAGE
  #   - THIS EXCLUDES MISSING ACCORDINGLY
  # ALSO COMPUTE THE DENOMINATOR TO BE USED BELOW
  if(NROW(cases)!=NROW(at)) stop("'x' and 'cases' not consistent")
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    if(is.null(ID)){
    den <- sum(rowMeans(cases,na.rm=TRUE),na.rm=TRUE)
    cases <- rowMeans(cases)
    }else{
      den <- tapply(1:nrow(cases), ID,
                    function(idx) sum(rowMeans(cases[idx, , drop=FALSE], na.rm = TRUE), na.rm = TRUE))
      cases <- rowMeans(cases)
    }
  } else {
    if(!is.null(ID)){
    den <- tapply(1:length(cases), ID,
                  function(idx) sum(cases[idx], na.rm = TRUE))
    }else{
      den <- sum(cases,na.rm=TRUE) 
    }
      
    if(dir=="forw") 
      cases <- rowMeans(as.matrix(tsModel:::Lag(cases,-seq(lag[1],lag[2]),group=group)))
  }
  
  ################################################################################
  #
  
  if(!is.null(model$knots_inter)){
    z = as.matrix(ps(z, knots = model$knots_inter))
  }
  
  # PREPARE THE ARGUMENTS FOR TH BASIS TRANSFORMATION
  predvar <- nrow(at) # number of predictions
  predlag <- lag[1]:lag[2]
  #  
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES 
  at_x_predvar = as.numeric(at)
  basisvar_predvar <- do.call("onebasis", c(list(x = at_x_predvar), attr(crossbasis,"argvar")))
  basislag_predvar <- do.call("onebasis", c(list(x = predlag), attr(crossbasis,"arglag")))
  
  # basis variables for center
  basiscen <- do.call("onebasis", c(list(x = cen), attr(crossbasis,"argvar")))
  
  # center data matrix
  Xpred_cen <- scale(basisvar_predvar,center = basiscen, scale = F)
  
  # Prediction matrix
  Xpred_predvar <- matrix(0, nrow=length(at_x_predvar), ncol = ncol(crossbasis))
  
  for (l in seq(length = length(predlag))){
    for (v in seq(length = ncol(Xpred_cen))){
      for (k in 1:ncol(basislag_predvar)){
        Xpred_predvar[((l-1)*predvar+1):(predvar*l),(ncol(basislag_predvar)*(v-1)+k)] = Xpred_cen[((l-1)*predvar+1):(predvar*l),v]*basislag_predvar[l,k] # prediction for every variable at lag = l
      }
    }
  }
  
  if(!is.null(z)){
  at_inter <- kronecker(Matrix(1, nrow = length(predlag), ncol = 1), z)  
  
  Mat_crossbasis_spat <- Matrix::Matrix(Xpred_predvar, sparse = TRUE) # cross basis
  
  interaction_list <- lapply(seq_len(ncol(Mat_crossbasis_spat)), function(j) {
    at_inter * Mat_crossbasis_spat[, j]
  })
  
  interaction_matrix <- do.call(cbind, interaction_list)
  
  # Global design mtarix
  Xpred_predvar <- Matrix::Matrix(cbind(Xpred_predvar,interaction_matrix), sparse = TRUE)
  
  }
  
  # Coefficients xi belonging to crossbasis
  
  Xpredall <- 0
  for (j in seq(length = length(predlag))) {
    ind_all <- seq(predvar) + predvar * (j - 1) # first lag period for all observations
    Xpredall <- Xpredall + Xpred_predvar[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
  }
  
  
  #  
  # CHECK DIMENSIONS  
  if(length(xi_mode)!=ncol(Xpredall))
    stop("arguments 'basis' do not match 'xi_mode'")
  if(any(dim(sigma_xi)!=c(length(xi_mode),length(xi_mode)))) 
    stop("arguments 'xi_mode' and 'sigma_xi' do not match")
  
  #
  ################################################################################
  #
  # COMPUTE AF AND AN 
  af <- 1-exp(-drop(as.matrix(Xpredall%*%xi_mode)))
  an <- af*cases
  #
  # TOTAL
  #   - SELECT NON-MISSING OBS CONTRIBUTING TO COMPUTATION
  #   - DERIVE TOTAL AF
  #   - COMPUTE TOTAL AN WITH ADJUSTED DENOMINATOR (OBSERVED TOTAL NUMBER)
  if(tot) {
    if(!is.null(ID)){
      an_tot = af_tot = NULL
      isna <- is.na(an)
      for (k in 1:length(unique(ID))){
          af_tot[k] <- sum(an[ID==unique(ID)[k] & !isna])/sum(cases[ID==unique(ID)[k] & !isna])
          an_tot[k] <- af_tot[k]*den[k]
      }
      an = an_tot
      af = af_tot
    }else{
      isna <- is.na(an)
      af <- sum(an[!isna])/sum(cases[!isna])
      an <- af*den
  }
  }
  #
  ################################################################################
  #
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    k <- length(xi_mode)
    eigen <- eigen(sigma_xi)
    X <- matrix(rnorm(length(xi_mode)*nsim),nsim)
    coefsim <- xi_mode + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    
    
    # RUN THE LOOP
    if(!is.null(ID)){
      afsim <- apply(coefsim,2, function(coefi) {
            ani = (1-exp(-drop(Xpredall%*%coefi)))*cases
            ani_result = NULL
            for (r in 1:length(unique(ID))){
                  ani_result[r] = sum(ani[ID==unique(ID)[r] & !is.na(ani)])/sum(cases[ID==unique(ID)[r] & !is.na(ani)])
            }
            return(ani_result)
        })
    
      ansim = matrix(NA, nrow = nrow(afsim), ncol = ncol(afsim))
      for (r in 1:nsim){
        ansim[,r] <- as.numeric(afsim[,r])*den
      }

    }else{
      afsim <- apply(coefsim,2, function(coefi) {
        ani <- (1-exp(-drop(Xpredall%*%coefi)))*cases
        sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
      })
      
      ansim <- afsim*den
  }}
  #
  ################################################################################
  #
  res <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an else af    
  }
  
  if (!is.null(ex.prob) & sim){
    #threshold = -log(1-ex.prob)
    #sd_all = sqrt(diag(drop(Xpredall) %*% sigma_xi %*% t(drop(Xpredall))))
    #exceedance_1 <- mapply(function(mean_val, sd_val) {
    #  pnorm(q = threshold, mean = mean_val, sd = sd_val, lower.tail = F)
    #},drop(as.matrix(Xpredall%*%xi_mode)), sd_all)
    if(!is.null(ID)){
      exceedance_1 = NULL
      for (r in 1:length(unique(ID))){
          exceedance_1[r] = sum(afsim[r,]>ex.prob, na.rm=T)/ncol(afsim)
      }
    }else{
      res <- list(res = af, probs = exceedance_1)
    }

    
  }
  return(res)
}

#

