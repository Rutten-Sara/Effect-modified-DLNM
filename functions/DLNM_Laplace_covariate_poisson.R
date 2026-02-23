DLNM_Laplace_cov_pois <- function(model,
                                  crossbasis,
                                  ID,
                                  covar.ri = NULL,
                                  data,
                                  map,
                                  type = 'non_linear',
                                  df_inter = 5,
                                  df_z = 1,
                                  smooth = NULL,
                                  df_smooth = 10,
                                  smooth2 = NULL,
                                  df_smooth2 = 10,
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
    if(is.vector(z)){
      knots_inter = NULL
      crossbasis_by = matrix(z, ncol = 1)
      df_inter = 1
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
      crossbasis_by = z
      main_z = z
      df_inter = 1
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
  required_packages <- c("Matrix", "spdep","sf","dlnm")
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


  Pv1 <- function(v) exp(v[1])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)) +
    exp(v[2])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl) +
    exp(v[3])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add)
  

  logpv.fixed <- function(v) 0.5 * nu * (v[1]+v[2]+v[3])- 
    (0.5*nu + a)*(log(b + 0.5*nu*exp(v[1]))+log(b + 0.5*nu*exp(v[2]))+
                    log(b + 0.5*nu*exp(v[3])))+ 
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
    Pv3 <- function(v) exp(v[8])*P_z
    
    
    logpv.fixed_z <- function(v) 0.5 * nu * v[8]- 
      (0.5*nu + a)*log(b + 0.5*nu*exp(v[8]))+ 
      0.5*determinant_Pv(Pv3(v))
    
    pen_dlnm = pen_dlnm + 1
  }
    
    Pv2 <- function(v) exp(v[4])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)) +
      exp(v[5])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)) +
      exp(v[6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL))+
      exp(v[7])*(Matrix::Diagonal(n = vx, x = NULL)%x%Matrix::Diagonal(n = vl, x = NULL)%x%P_inter)
    
    
    logpv.inter <- function(v) 0.5 * nu * (v[5]+v[6]+v[7]+v[4])- 
      (0.5*nu + a)*(log(b + 0.5*nu*exp(v[5]))+
                      log(b + 0.5*nu*exp(v[6]))+log(b + 0.5*nu*exp(v[7]))+log(b + 0.5*nu*exp(v[4])))+ 
      0.5*determinant_Pv(Pv2(v))

    
    
    pen_dlnm = pen_dlnm + 4
    
  }else if (type != "none"){
    if(df_z>1){
      D_z <- Matrix::Diagonal(dim(main_z)[2], x = NULL)
      for (k in 1:penorder) D_z <- Matrix::diff(D_z)
      P_z <- Matrix::t(D_z) %*% D_z
      P_z <- P_z + Matrix::Diagonal(dim(main_z)[2], x = 1e-12)
      Pv3 <- function(v) exp(v[7])*P_z
      
      
      logpv.fixed_z <- function(v) 0.5 * nu * v[7]- 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[7]))+ 
        0.5*determinant_Pv(Pv3(v))
      
      pen_dlnm = pen_dlnm + 1
    }
    
    
    Pv2 <- function(v) exp(v[4])*(Px%x%Matrix::Diagonal(n = vl, x = NULL)%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)) +
      exp(v[5])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL)) +
      exp(v[6])*(Matrix::Diagonal(n = vx, x = NULL)%x%Pl_add%x%Matrix::Diagonal(n = dim(crossbasis_by)[2], x = NULL))
      
      
      logpv.inter <- function(v) 0.5 * nu * (v[4]+v[5]+v[6])- 
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[4]))+log(b + 0.5*nu*exp(v[5]))+
                        log(b + 0.5*nu*exp(v[6])))+ 
        0.5*determinant_Pv(Pv2(v))
      
      
      pen_dlnm = pen_dlnm + 3
  }

  }
  
  zeta <- 1e-05 # precision for linear effect coefficient
  # Hyperparameters for penalty parameters
  a <- b <- 10^(-5)
  nu <- 3
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
      
      Gv <- function(v) Matrix::Diagonal(n = q.rand, x = exp(v[pen_dlnm+1]))
      logpv.rand <- function(v) 0.5 * (nu + q.rand) * v[pen_dlnm+1] - 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_dlnm+1]))
      
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
      
      Gv <- function(v) exp(v[pen_dlnm+1])*Rn
      logpv.rand <- function(v) 0.5 * (nu + q.rand) * v[pen_dlnm+1] - 
        (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_dlnm+1]))
      
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
      
      Gv <- function(v) Matrix::bdiag(exp(v[pen_dlnm+1])*Rn,
                                      Matrix::Diagonal(n = q.rand, x = exp(v[pen_dlnm+1])))
      logpv.rand <- function(v) sum(0.5 * (nu + q.rand) * v[(pen_dlnm+1):(pen_dlnm+2)]) - 
        sum((0.5*nu + a)*log(b + 0.5*nu*exp( v[(pen_dlnm+1):(pen_dlnm+2)])))
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
      
      Gv <- function(v) exp(v[pen_dlnm+1])*(Matrix::Diagonal(n = q.rand, 
                                                    x = 1 - exp(v[pen_dlnm+2])/(1+exp(v[pen_dlnm+2]))) + 
                                     Matrix::Matrix(exp(v[pen_dlnm+2])/(1+exp(v[pen_dlnm+2]))*Rn, sparse = T))
      logpv.rand <- function(v)  {
        value <- 0.5 * nu * v[pen_dlnm+1] - (0.5*nu + a)*log(b + 0.5*nu*exp(v[pen_dlnm+1])) + 
          0.5*sum(sapply(eigen(Gv(v),only.values = T)$values,log)) + 
          a.rho*v[pen_dlnm+2] - (a.rho + b.rho)*log(1 + exp(v[pen_dlnm+2]))
        return(as.numeric(value))}
      v.rand <- c(0,1)
    }
    
    # Global design mtarix
    X <- Matrix::Matrix(cbind(X, Z.rand), sparse = TRUE)

  } else {
      v.rand <- NULL

  }
  
  logpv.smooth <- function(v) 0

  
  if(!is.null(smooth)){
    smooth <- smooth[!is.na(crossbasis[,1])]
    smooth_cov <- dlnm::ps(smooth, df = df_smooth, intercept = F)
    
    D_smooth = Matrix::t(Matrix::Matrix(attr(smooth_cov,"S")))+ Matrix::Diagonal(n = df_smooth, x = 1e-12)
    X_smooth <- Matrix::Matrix(smooth_cov)

    
    if(!is.null(smooth2)){
      smooth2 <- smooth2[!is.na(crossbasis[,1])]
      smooth2_cov <- dlnm::ps(smooth2, df = df_smooth2, intercept = F)
      
      D_smooth2 = Matrix::t(Matrix::Matrix(attr(smooth2_cov,"S")))+ Matrix::Diagonal(n = df_smooth2, x = 1e-12)
      Pv_smooth <- function(v) as(Matrix::bdiag(exp(v[length(v)-1])*D_smooth,
                                                exp(v[length(v)])*D_smooth2), "generalMatrix")
      X_smooth <- Matrix::Matrix(cbind(X_smooth, Matrix::Matrix(smooth2_cov)), sparse = TRUE)
      
      
      logpv.smooth <- function(v) 0.5 * nu * (v[length(v)]+v[length(v)-1])- 
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[length(v)]+v[length(v)-1])))+ 
        0.5*determinant_Pv(Pv_smooth(v))
    }else{
      
      Pv_smooth <- function(v) exp(v[length(v)])*as(D_smooth, "generalMatrix")
      
      logpv.smooth <- function(v) 0.5 * nu * (v[length(v)])- 
        (0.5*nu + a)*(log(b + 0.5*nu*exp(v[length(v)])))+ 
        0.5*determinant_Pv(Pv_smooth(v))
    }

    X <- Matrix::Matrix(cbind(X, X_smooth), sparse = TRUE) 
  }
  
  
  #Precision matrix for parameter xi
  if(is.null(smooth)){
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
 
  }else{
    if(type == "none" & is.null(covar.ri)){
      Qv <- function(v){
        result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                           x = zeta),
                                          Pv1(v), Pv_smooth(v)), "generalMatrix")
      }} else if (type == "none"& !is.null(covar.ri)){
        Qv <- function(v){
          result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                             x = zeta),
                                            Pv1(v),
                                            Gv(v), Pv_smooth(v)), "generalMatrix")
        }}else if(is.null(covar.ri) ){
          Qv <- function(v){
            result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                               x = zeta),
                                              Pv3(v),
                                              Pv1(v),
                                              Pv2(v), Pv_smooth(v)), "generalMatrix")
          }} else{
            Qv <- function(v){
              result_matrix <- as(Matrix::bdiag(Matrix::Diagonal(n = dim(Z.linear)[2],
                                                                 x = zeta), 
                                                Pv3(v),
                                                Pv1(v),
                                                Pv2(v),
                                                Gv(v), Pv_smooth(v)), "generalMatrix")
            }
          }
    
    
  }
  # For Poisson GLM with log-link
  Cv_xi <- function(xi, Xv){
      Cvxi <- exp(as.numeric(Xv %*% xi)+log(offset))
      return(Cvxi)
  }
    

  log_pxi <- function(xi, Qv, Cvxi) {
    post <- sum(y * log(Cvxi) - Cvxi) - .5 * Matrix::t(xi) %*% Qv %*% xi
    return(as.numeric(post))
  }
  
  Grad.logpxi <- function(xi, Qv, Cvxi){
    value <- Matrix::t(X)%*%(y-Cvxi) - Qv%*%xi
    as.numeric(value)
  }
  

  # Laplace approximation to conditional posterior of xi
  # using Newton-Raphson algorithm

  NR_xi <- function(xi0, Qv0){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 300   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      Cvxi0 <- Cv_xi(xi0, X)
      grad_est =Grad.logpxi(xi0, Qv0, Cvxi0)
      dxi <- as.numeric((-1)* solve_sparse_cholesky(Qv0, Cvxi0, X, grad_est))

      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi0, Qv0, Cvxi0)
      while (max(log_pxi(xi.new, Qv0, Cv_xi(xi.new, X)),logpxi.current, na.rm=T)==logpxi.current|
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
    
    #if(iter == maxiter){
    #  print("Warning: algorithm did not converge")
    #}
    
    
    xistar <- xi0
    return(xistar)
  }

  
  # Initial values for log-penalty and log-overdispersion parameter

  dimxi <- dim(X)[2]
  xi0_init = rep(0,dimxi)
  

  if(is.null(pen_crossbasis) & is.null(covar.ri)){
  v_init = 1
  }else{
    v_init = c(rep(2,pen_dlnm),v.rand)
  }

  if(!is.null(smooth)){
    v_init <- c(v_init,5)
    if(!is.null(smooth2)){
      v_init <- c(v_init,2)
    }
  }

  
  #if(approx == F){
  Qv_init <- Qv(v_init)
  # }else{
  #  Qv_init <- Qv_approx(v_init)
  #}

  xi_init <- NR_xi(xi0 = xi0_init, Qv0 = Qv_init)
  
  Cvxi_temp <- Cv_xi(xi_init, X)
  # Log-posterior for the penalty- vector
  XWX <- XWX_func(Cvxi_temp, X)
  Xxi <- X%*%xi_init + log(offset)

 
  if(length(v_init)>0){
  log_pv <- function(v){
    
    result <- tryCatch({
      Qv_est <- Qv(v)
      a1 <- -0.5*determinant_Pv(XWX+Qv_est)
      
      a2 <- sum(y*(Xxi))
      a3 <- sum(exp(Xxi))
      
      a4 <- 0.5 * sum((xi_init * Qv_est) %*% xi_init)
      
      as.numeric(value <- a1+a2-a3-a4+logpv.fixed(v)+logpv.rand(v) + logpv.smooth(v) + logpv.inter(v)+logpv.fixed_z(v))
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
  


  #v_init = v_mode
  #xi0_init = xi_init
  # }
  
  
  # Estimate of xi (regression parameter)
  #  if(approx == F){
    Qv_mode <- Qv(v_mode)
    # }else{
    #   Qv_mode <- Qv_approx(v_mode)
    # }

  xi_mode <- NR_xi(xi0 = xi_init, Qv0 = Qv_mode)
  }else{
    xi_mode = xi_init
  }
  
  
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

  Prec <- - Hess_logpxi(Qv_mode,Cv_xi(xi_mode, X), X)
  Sigma = Matrix::solve(Prec)
  
  if(DIC == T){
    DIC_pd = calculate_DIC_tot(xi_mode, Prec, Qv_mode,X, y, offset, v_mode, model.family = "poisson")
    DIC_result = DIC_pd$DIC
    pd_result = DIC_pd$pd
  }else{
    DIC_result = NULL
    pd_result = NULL
  }
 


  if(!is.null(covar.ri)){
    if(is.null(smooth)){
    if(covar.ri == "Convolution"){
      str <- tail(xi_mode, 2*q.rand)[1:q.rand]
      unstr <- tail(xi_mode, q.rand)
      xispat <- (str + unstr) - mean(str + unstr)
    } else {
      xispat <- tail(xi_mode, q.rand) - mean(tail(xi_mode, q.rand))
    } 
    }else{
      if(covar.ri == "Convolution"){
        str <- tail(xi_mode[1:(length(xi_mode)-dim(X_smooth)[2])], 2*q.rand)[1:q.rand]
        unstr <- tail(xi_mode[1:(length(xi_mode)-dim(X_smooth)[2])], q.rand)
        xispat <- (str + unstr) - mean(str + unstr)
      } else {
        xispat <- tail(xi_mode[1:(length(xi_mode)-dim(X_smooth)[2])], q.rand) - mean(tail(xi_mode[1:(length(xi_mode)-dim(X_smooth)[2])], q.rand))
      } 
  }
    }else {
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
