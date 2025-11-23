DLNM_Laplace_cov_NB <- function(model,
                                crossbasis,
                                ID,
                                covar.ri = NULL,
                                data,
                                map,
                                type = 'non_linear',
                                df_inter = 5,
                                df_z = 1,
                                offset = NULL,
                                z = NULL, 
                                pen_crossbasis = NULL, DIC = T){
  
  
  if(type == "non_linear"){
    crossbasis_by = ps(z, df = df_inter,
                       intercept = F)
    knots_inter = attr(crossbasis_by, "knots")
    crossbasis_by = as.matrix(crossbasis_by)
    if(df_z>1){
      main_z = ps(z, df = df_z,
                  intercept = F)
      knots_z = attr(main_z, "knots")
      main_z = as.matrix(main_z)
    }else{
      main_z = matrix(z, ncol = 1)
      knots_z = NULL
    }
    
  }else if (!type=="none"){
    df_inter = 1
    if(is.vector(z)){
      knots_inter = NULL
      crossbasis_by = matrix(z, ncol = 1)
      if(df_z>1){
        main_z = ps(z, df = df_z,
                    intercept = F)
        knots_z = attr(main_z, "knots")
        main_z = as.matrix(main_z)
      }else{
        main_z = matrix(z, ncol = 1)
        knots_z = NULL
      }
      type = "linear"
    }else{
      df_inter = 1
      crossbasis_by = z
      main_z = z
      type = "factor"
      knots_z = NULL
      knots_inter = NULL
    }
  }else{
    knots_z = NULL
    knots_inter= NULL
    crossbasis_by = NULL
  }
  
  # Required packages
  required_packages <- c("Matrix", "spdep","sf")
  for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message(paste("Package", package, "is required but not installed."))
    }
  }

  
  # If no population offset --> put equal to 1 for every area
  if (is.null(offset)){
    offset <- rep(1, dim(crossbasis)[1])
  }
  
  vx <- attributes(crossbasis)$df[1]
  vl <-  attributes(crossbasis)$df[2]
  
  offset <- offset[!is.na(crossbasis[,1])] 
  if(!type=="none"){
    crossbasis_by <- as.matrix(crossbasis_by[!is.na(crossbasis[,1]),])
    main_z = matrix(main_z[!is.na(crossbasis[,1]), ], ncol = ncol(main_z))
  }
  
  if(!missing(data)){
    mod <- stats::model.frame(model, data = data, na.action = na.pass)
    na_ind = (rowSums(is.na(mod))>0)
    
    Z.linear_part <- Matrix::sparse.model.matrix(mod, data = data, drop.unused.levels = T)
    Z.linear <- Matrix::sparseMatrix(
      i = rep(which(!na_ind), times = ncol(Z.linear_part)),
      j = rep(seq_len(ncol(Z.linear_part)), each = length(which(!na_ind))),
      x = as.numeric(Z.linear_part),
      dims = c(length(na_ind), ncol(Z.linear_part))
    )
    y <- as.numeric(stats::model.extract(mod, "response"))[!is.na(crossbasis[,1])] # response
    
    Z.linear = Z.linear[!is.na(crossbasis[,1]), ,drop = FALSE] # linear covariates
    
    p <- ncol(Z.linear)
    
  } else{
    mod <- stats::model.frame(model, na.action = na.pass)
    
    na_ind = (rowSums(is.na(mod))>0)
    
    Z.linear_part <- Matrix::sparse.model.matrix(mod, drop.unused.levels = T)
    Z.linear <- Matrix::sparseMatrix(
      i = rep(which(!na_ind), times = ncol(Z.linear_part)),
      j = rep(seq_len(ncol(Z.linear_part)), each = length(which(!na_ind))),
      x = as.numeric(Z.linear_part),
      dims = c(length(na_ind), ncol(Z.linear_part))
    )
    
    y <- as.numeric(stats::model.extract(mod, "response"))[!is.na(crossbasis[,1])] # response
    Z.linear = Z.linear[!is.na(crossbasis[,1]), ,drop = FALSE]
    p <- ncol(Z.linear)

  }
  
  
  Wcb <- Matrix::Matrix(na.omit(crossbasis), sparse = TRUE) # cross basis
  
  if(type == "none"){
    # Global design mtarix
    X <- Matrix::Matrix(cbind(Z.linear, Wcb), sparse = TRUE)
    interaction_matrix <- matrix(1, ncol = 1, nrow = 1)
    main_z <- matrix(1, ncol = 1, nrow = 1)
  }else{
    interaction_list <- lapply(seq_len(ncol(Wcb)), function(j) {
      crossbasis_by * Wcb[, j]
    })
    
    interaction_matrix <- do.call(cbind, interaction_list)
    
    # Global design mtarix
    X <- Matrix::Matrix(cbind(Z.linear,main_z, Wcb,interaction_matrix), sparse = TRUE)
  }
  
  
  logpv.fixed <- function(v) 0
  logpv.inter <- function(v) 0
  logpv.fixed_z <- function(v) 0
  
  pen_dlnm <- 0
  
  Pv1 <- function(v){
    as(Matrix::Diagonal(n = dim(Wcb)[2],
                        x = zeta), "generalMatrix")
  }
  
  Pv2 <- function(v){
    as(Matrix::Diagonal(n = dim(interaction_matrix)[2],
                        x = zeta), "generalMatrix")
  }
  
  Pv3 <- function(v){
    as(Matrix::Diagonal(n = dim(main_z)[2],
                        x = zeta), "generalMatrix")
  }
  
  
  # Difference order of the penalty
  if(!is.null(pen_crossbasis)){
    penorder <- pen_crossbasis
    
    # Penalty for exposure
    Dx <- Matrix::Diagonal(vx+1, x = NULL)
    for (k in 1:penorder) Dx <- Matrix::diff(Dx)
    Px <- Matrix::t(Dx) %*% Dx
    Px <- Px[-1,-1]
    Px <- Px + Matrix::Diagonal(vx, x = 1e-12)
    
    # Penalty for lag variable
    Dl <- Matrix::Diagonal(vl, x = NULL)
    for (k in 1:penorder) Dl <- Matrix::diff(Dl)
    Pl <- Matrix::t(Dl) %*% Dl
    Pl <- Pl + Matrix::Diagonal(vl, x = 1e-12)
    
    # Additional penalty for delay (force lag effect going to zero)
    Dl_add <- diag((0:(vl-1))^2)
    #Dl_add <- diag(rep(0:1,c(6,4)))
    Pl_add <- Dl_add + diag(1e-12, vl)
    
    
    
    Pv1 <- function(v) exp(v[2])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)) +
      exp(v[3])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl) +
      exp(v[4])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add)
    
    
    logpv.fixed <- function(v) 0.5 * nu * (v[2]+v[3]+v[4])- 
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[2]))+log(b + 0.5*nu*exp(v[3]))+
                      log(b + 0.5*nu*exp(v[4])))+ 
      0.5*determinant_Pv(Pv1(v))
    
    pen_dlnm = pen_dlnm + 3
    
    if(type == "non_linear"){
      
      D_inter <- Matrix::Diagonal(dim(crossbasis_by)[2], x = NULL)
      for (k in 1:penorder) D_inter <- Matrix::diff(D_inter)
      P_inter <- Matrix::t(D_inter) %*% D_inter
      P_inter <- P_inter + Matrix::Diagonal(dim(crossbasis_by)[2], x = 1e-12)
      
      
      if(df_z>1){
        D_z <- Matrix::Diagonal(dim(main_z)[2], x = NULL)
        for (k in 1:penorder) D_z <- Matrix::diff(D_z)
        P_z <- Matrix::t(D_z) %*% D_z
        P_z <- P_z + Matrix::Diagonal(dim(main_z)[2], x = 1e-12)
        Pv3 <- function(v) exp(v[9])*P_z
        
        
        logpv.fixed_z <- function(v) 0.5 * nu * v[9]- 
          (0.5*nu + a)*log(b + 0.5*nu*exp(v[9]))+ 
          0.5*determinant_Pv(Pv3(v))
        
        pen_dlnm = pen_dlnm + 1
      }
      
      Pv2 <- function(v) exp(v[5])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)) +
        exp(v[6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)) +
        exp(v[7])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL))+
        exp(v[8])*(Matrix::Diagonal(n = vx, x = NULL)%x%Matrix::Diagonal(n = vl, x = NULL)%x%P_inter)
      
      
      logpv.inter <- function(v) 0.5 * nu * (v[5]+v[6]+v[7]+v[8])- 
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[5]))+
                        log(b + 0.5*nu*exp(v[6]))+log(b + 0.5*nu*exp(v[7]))+log(b + 0.5*nu*exp(v[8])))+ 
        0.5*determinant_Pv(Pv2(v))

      
      
      pen_dlnm = pen_dlnm + 4
      
    }else if (type != "none"){
      
      if(df_z>1){
        D_z <- Matrix::Diagonal(dim(main_z)[2], x = NULL)
        for (k in 1:penorder) D_z <- Matrix::diff(D_z)
        P_z <- Matrix::t(D_z) %*% D_z
        P_z <- P_z + Matrix::Diagonal(dim(main_z)[2], x = 1e-12)
        Pv3 <- function(v) exp(v[8])*P_z
        
        
        logpv.fixed_z <- function(v) 0.5 * nu * v[8]- 
          (0.5*nu + a)*log(b + 0.5*nu*exp(v[8]))+ 
          0.5*determinant_Pv(Pv3(v))
        
        pen_dlnm = pen_dlnm + 1
      }
      
      Pv2 <- function(v) exp(v[5])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)) +
        exp(v[6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)) +
        exp(v[7])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL))
      
        
        logpv.inter <- function(v) 0.5 * nu * (v[7]+v[5]+v[6])- 
          (0.5*nu + a)*(log(b + 0.5*nu*exp(v[7]))+log(b + 0.5*nu*exp(v[5]))+
                          log(b + 0.5*nu*exp(v[6])))+ 
          0.5*determinant_Pv(Pv2(v))
        
        
        pen_dlnm = pen_dlnm + 3
    }
    
  }
  
  zeta <- 1e-05 # precision for linear effect coefficient
  # Hyperparameters for penalty parameters
  a <- b <- 10^(-5)
  nu <- 3
  
  # Prior for overdispersion parameter phi
  a.disp <- 1e-05
  b.disp <- 1e-05
  
  logpv.rand <- function(v) 0
  
  
  
  if(!is.null(covar.ri)){
    atau <- btau <- 10^(-5)
    a.rho <- b.rho <- 0.5
    if(!missing(data)){
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID) + 0, data)[!is.na(crossbasis[,1]),] # model matrix for random effects
    } else {
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID) + 0)[!is.na(crossbasis[,1]),] # model matrix for random effects
    }
    q.rand <- dim(Z.rand)[2]
    v.rand <- 0
    
    
    if(covar.ri == "ind"){
      
      Gv <- function(v) Matrix::Diagonal(n = q.rand, x = exp(v[1+pen_dlnm+1]))
      logpv.rand <- function(v) 0.5 * (nu + q.rand) * v[1+pen_dlnm+1] - 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[1+pen_dlnm+1]))
      
    } else if (covar.ri == "ICAR"){
      neig.map <- spdep::poly2nb(map)
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
      
      Gv <- function(v) exp(v[1+pen_dlnm+1])*Rn
      logpv.rand <- function(v) 0.5 * (nu + q.rand) * v[1+pen_dlnm+1] - 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[1+pen_dlnm+1]))
      
    } else if (covar.ri == "Convolution") {
      neig.map <- spdep::poly2nb(map)
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
      
      Gv <- function(v) Matrix::bdiag(exp(v[1+pen_dlnm+1])*Rn,
                                      Matrix::Diagonal(n = q.rand, x = exp(v[1+pen_dlnm+1])))
      logpv.rand <- function(v) sum(0.5 * (nu + q.rand) * v[(1+pen_dlnm+1):(1+pen_dlnm+2)]) - 
        sum((0.5*nu + a)*log(b + 0.5*nu*exp( v[(1+pen_dlnm+1):(1+pen_dlnm+2)])))
      Z.rand <- cbind(Z.rand, Z.rand)
      v.rand <- c(0,1)
      
    } else if (covar.ri == "Leroux"){
      neig.map <- spdep::poly2nb(map)
      Rn <- matrix(0, nrow = q.rand, ncol =q.rand)
      
      for (s in 1:q.rand) {
        # Diagonal elements (N_s)
        Rn[s, s] <- length(neig.map[[s]])
        
        # Off-diagonal elements (-1 for neighbors)
        for (u in neig.map[[s]]) {
          Rn[s, u] <- -1
        }
      }
      Rn <- Matrix::Matrix(Rn, sparse = TRUE)
      
      Gv <- function(v) exp(v[1+pen_dlnm+1])*(Matrix::Diagonal(n = q.rand, 
                                                             x = 1 - exp(v[1+pen_dlnm+2])/(1+exp(v[1+pen_dlnm+2]))) + 
                                              Matrix::Matrix(exp(v[1+pen_dlnm+2])/(1+exp(v[1+pen_dlnm+2]))*Rn, sparse = T))
      logpv.rand <- function(v)  {
        value <- 0.5 * nu * v[1+pen_dlnm+1] - (0.5*nu + a)*log(b + 0.5*nu*exp(v[1+pen_dlnm+1])) + 
          0.5*sum(sapply(eigen(Gv(v),only.values = T)$values,log)) + 
          a.rho*v[1+pen_dlnm+2] - (a.rho + b.rho)*log(1 + exp(v[1+pen_dlnm+2]))
        return(as.numeric(value))}
      v.rand <- c(0,1)
    }
    
    # Global design mtarix
    X <- Matrix::Matrix(cbind(X, Z.rand), sparse = TRUE)
    
  } else {
    v.rand <- NULL
    
  }
  
  
  
  #Precision matrix for parameter xi
  if(type == "none" & is.null(covar.ri)){
    Qv <- function(v){
      result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                         x = zeta),
                                        Pv1(v)), "generalMatrix")
    }} else if (type == "none"& !is.null(covar.ri)){
    Qv <- function(v){
      result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                         x = zeta),
                                        Pv1(v),
                                        Gv(v)), "generalMatrix")
    }}else if(is.null(covar.ri) ){
    Qv <- function(v){
      result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                         x = zeta),
                                        Pv3(v),
                                        Pv1(v),
                                        Pv2(v)), "generalMatrix")
    }} else{
      Qv <- function(v){
        result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                           x = zeta), 
                                          Pv3(v),
                                          Pv1(v),
                                          Pv2(v),
                                          Gv(v)), "generalMatrix")
      }
    }
  
  
  # Negative binomial GLM with log-link
  Cv_xi <- function(xi, Xv){
    Cvxi <- exp(as.numeric(Xv %*% xi)+log(offset))
    return(Cvxi)
  }
  
  var.nb <- function(xi, Cvxi,v) {
    res <- Cvxi + (1 / exp(v[1])) * (Cvxi ^ 2)
    return(res)
  }
  W.nb <- function(xi,Cvxi, v) {
    varval <- Cvxi + (1 / exp(v[1])) * (Cvxi ^ 2)
    res <- Matrix::Diagonal(x = ((Cvxi) ^ 2) * (1 / varval))
    return(res)
  }
  
  D.nb <- function(xi, Cvxi) Matrix::Diagonal(x = 1/Cvxi)
  M.nb <- function(xi, Cvxi) Matrix::Diagonal(x = y - Cvxi)
  V.nb <- function(xi, Cvxi, v){
    varval <-  Cvxi + (1 / exp(v[1])) * (Cvxi ^ 2)
    res <- Matrix::Diagonal(x = Cvxi * (1/varval - (Cvxi/(varval^2)) * 
                                           (1 + 2*Cvxi*(1/exp(v[1])))))
    return(res)
  }
  gamma.nb <- function(xi, Cvxi,v) {
    res <- exp(v[1]) * log(Cvxi / (Cvxi + exp(v[ 1])))
    return(res)
  }
  bgamma.nb <- function(xi, Cvxi, v) - (exp(v[ 1])^2) *
    log(exp(v[1])/(exp(v[1]) + Cvxi))
  
  # Log conditional posterior of xi given v
  log_pxi <- function(xi, Qv, Cvxi,v) {
    value <- (1/exp(v[1])) * sum((y * gamma.nb(xi, Cvxi,v)) - bgamma.nb(xi, Cvxi,v)) -
      .5 * Matrix::t(xi) %*% Qv %*% xi
    return(as.numeric(value))
  }
  
  Grad.logpxi <- function(xi,Xv,Qv, Cvxi,v){
    varval <- Cvxi + (1 / exp(v[1])) * (Cvxi ^ 2)
    W.nbval <- Matrix::Diagonal(x = ((Cvxi) ^ 2) * (1 / varval))
    D.nbval <- Matrix::Diagonal(x = 1/Cvxi)
    value <- as.numeric(Matrix::t(Xv)%*%W.nbval%*%D.nbval%*%(y - Cvxi) -
                          Qv%*%xi)
    return(value)
  }


  
  Grad.logpxi_poisson <- function(xi, Qv, Cvxi){
    value <- Matrix::t(X)%*%(y-Cvxi) - Qv%*%xi
    as.numeric(value)
  }
  
  # Laplace approximation to conditional posterior of xi
  # using Newton-Raphson algorithm
  
  NR_xi <- function(xi0, Qv0,v0){
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 1000   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      Cvxi0 <- Cv_xi(xi0, X)
      grad_est =Grad.logpxi(xi0,X,Qv0, Cvxi0,v0)
      
      W.nb0 = as(W.nb(xi0,Cvxi0,v0),"generalMatrix")
      MV.nb0 = as(M.nb(xi0, Cvxi0)%*%V.nb(xi0,Cvxi0,v0),"generalMatrix")
      
      #dxi <- as.numeric((-1)* solve_sparse_cholesky_NB(Qv0, W.nb0, MV.nb0, as(Matrix::Diagonal(x = Matrix::diag(X)), "generalMatrix"), grad_est))
      dxi <- as.numeric((-1)* solve_sparse_cholesky_NB(Qv0, W.nb0, MV.nb0, X, grad_est))

      
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi0, Qv0, Cvxi0,v0)
      while (max(log_pxi(xi.new, Qv0, Cv_xi(xi.new, X),v0),logpxi.current, na.rm=T)==logpxi.current |
             max(abs(step * dxi))>10) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((xi.new - xi0) ^ 2))
      iter <- iter + 1
      xi0 <- xi.new
      if(dist < epsilon) break
    }
    
    
    if(iter == maxiter){
      print("Warning: algorithm did not converge")
    }
    
    xistar <- xi0
    return(xistar)
  }
  
  
  # Initial values for log-penalty and log-overdispersion parameter

  dimxi <- dim(X)[2]
  xi0_init = rep(0,dimxi)
  
  if(is.null(pen_crossbasis) & is.null(covar.ri)){
    v_init = 0
  }else{
    v_init = c(0,rep(2,pen_dlnm),v.rand)
  }

  #if(approx == F){
  Qv_init <- Qv(v_init)
  # }else{
  #  Qv_init <- Qv_approx(v_init)
  #}

  xi_init <- NR_xi(xi0 = xi0_init, Qv0 = Qv_init, v0 = v_init)
  
  Cvxi_temp <- Cv_xi(xi_init, X)
  
  #X_low_rank = svd_k_irlba(X, ceiling(dim(X)[2]*0.2))
  
  
  # Log-posterior for the penalty- vector
  Mnb_v <- M.nb(xi_init, Cvxi_temp)

  log_pv <- function(v){
    
    result <- tryCatch({
      vdisp <- v[1]
      Qv_est <- Qv(v)

      # Loglikelihood
 
      varval <-  Cvxi_temp + (1 / exp(vdisp)) * (Cvxi_temp ^ 2)
      Vnb <- Matrix::Diagonal(x = Cvxi_temp * (1/varval - (Cvxi_temp/(varval^2)) *
                                             (1 + 2*Cvxi_temp*(1/exp(vdisp)))))
      Wnb <- Matrix::Diagonal(x = ((Cvxi_temp) ^ 2) * (1 / varval))
      gammanb <- exp(vdisp) * log(Cvxi_temp / (Cvxi_temp + exp(vdisp)))
      bgammanb <- (-1) * (exp(vdisp)^2) * log(exp(vdisp)/(exp(vdisp) + Cvxi_temp))
      
      a2 <- sum((1/exp(vdisp))*((y * gammanb) - bgammanb) +
                  lgamma(y + exp(vdisp)) - lgamma(exp(vdisp)))
      
      
      # Determinant of posterior hessian
      W_v_opt = Matrix::diag(-(Mnb_v%*%Vnb-Wnb))
      #XWX_est = XWX_func(W_v_opt, X)
      # XWX_est <-  as(approx_XWX_from_svd(X_low_rank$U, X_low_rank$S, X_low_rank$V, W_v_opt),
      #               "CsparseMatrix")
      
      XWX_est = XWX_func(W_v_opt, as(Matrix::Diagonal(x = Matrix::diag(X)), "generalMatrix"))
      
      a1 <- -0.5*determinant_Pv(XWX_est+Qv_est)
      
      
      # Prior Q calculation
      a4 <- 0.5 * sum((xi_init * Qv_est) %*% xi_init)
      
      # Gamma prior overdispersion
      #a5 <- (0.5 * nu) * vdisp - ((0.5 * nu) + a.disp) *  log(0.5*(nu * exp(vdisp)) + b.disp)
      a5 <- a.disp*vdisp - b.disp*exp(vdisp)
      
      as.numeric(value <- a1+a2-a4+a5+logpv.fixed(v)+logpv.rand(v) + logpv.inter(v)+logpv.fixed_z(v))
    }, error = function(e) {
      -1e10  # penalize the optimizer
    })

    return(result)
  }

  
  # Mode a posteriori estimate of v
  v_mode <- optim(par = v_init, 
                  fn = log_pv, 
                  method="Nelder-Mead", 
                  control = list(fnscale = -1, reltol = 1e-12))$par
  

  Qv_mode <- Qv(v_mode)

  xi_mode <- NR_xi(xi0 = xi_init, Qv0 = Qv_mode, v0 = v_mode)
  
  
  
  
  if(type == "none"){
    ind.cb <- seq(dim(Z.linear)[2]+1, by = 1, length.out = dim(Wcb)[2])  
    xi_mode_cb <- xi_mode[ind.cb]
    ind.cb_inter = NULL
    xi_mode_cb_inter = NULL
  }else{
    ind.cb <- seq(dim(Z.linear)[2]+dim(main_z)[2]+1, by = 1, length.out = dim(Wcb)[2])  
    xi_mode_cb <- xi_mode[ind.cb]
    
    ind.cb_inter <- seq(dim(Z.linear)[2]+dim(main_z)[2]+dim(Wcb)[2]+1, by = 1, length.out = dim(interaction_matrix)[2])  
    xi_mode_cb_inter <- xi_mode[ind.cb_inter]
  }

  
  
  
  Cvxi_mode = Cv_xi(xi_mode, X)
  W.nb_mode = as(W.nb(xi_mode,Cvxi_mode,v_mode),"generalMatrix")
  MV.nb_mode = as(M.nb(xi_mode, Cvxi_mode)%*%V.nb(xi_mode,Cvxi_mode,v_mode),"generalMatrix")


  Prec <- - Hess_logpxi_NB(Qv_mode,W.nb_mode, MV.nb_mode, X)
  Sigma = Matrix::solve(Prec)
  
  if(DIC == T){
    DIC_pd = calculate_DIC_tot(xi_mode, Prec, Qv_mode,X, y, offset, v_mode, model.family = "NB")
    DIC_result = DIC_pd$DIC
    pd_result = DIC_pd$pd
  }else{
    DIC_result = NULL
    pd_result = NULL
  }
 


  if(!is.null(covar.ri)){
    if(covar.ri == "Convolution"){
      str <- tail(xi_mode, 2*q.rand)[1:q.rand]
      unstr <- tail(xi_mode, q.rand)
      xispat <- (str + unstr) - mean(str + unstr)
    } else {
      xispat <- tail(xi_mode, q.rand) - mean(tail(xi_mode, q.rand))
    } 
  } else {
    xispat <- NULL
  }

  
  
  
  
  
  
  output <- list("xi_mode" = xi_mode,
                 "Sigma" =  Sigma,
                 "ind.cb" = ind.cb,
                 "xi_mode_cb" = xi_mode_cb,
                 "ind.cb_inter" = ind.cb_inter,
                 "xi_mode_cb_inter" = xi_mode_cb_inter,
                 "crossbasis" = crossbasis,
                 "crossbasis_by" = crossbasis_by,
                 "xispat" = xispat,
                 "v_mode" = v_mode,
                 "DIC" = DIC_result,
                 "pd" = pd_result,
                 "knots_z" = knots_z,
                 "knots_inter" = knots_inter)


  
  
  }
