
################################################################################
# Function to fit an effect-modified DLNM model
################################################################################
#   - model: model formula (without crossbasis or effect modifier)
#   - crossbasis: Cross-basis computed from x
#   - ID: vector of ID of each observation (for random effect)
#   - covar.ri : Random effect structure (ind, ICAR, Convolution, Leroux or NULL)
#   - data: dataset
#   - map: only needed if spatial structure on random effect
#   - type: none, non_linear, linear or factor modification
#   - df_inter: degrees of freedom for effect modifier in interaction effect
#   - df_z: degrees of freedom for main effect of effect modifier
#   - offset: offset of the model
#   - z: vector with values of effect modifier
#   - pen_crossbasis: difference penalty order on crossbasis (if penalty)
#   - DIC: calculate DIC criterium (true or false)
#   - family: "poisson" or "NB"
################################################################################

DLNM_Laplace_cov <- function(model,
                             crossbasis,
                             ID,
                             covar.ri = NULL,
                             data,
                             map = NULL,
                             type = "non_linear",
                             df_inter = 5,
                             df_z = 1,
                             offset = NULL,
                             z = NULL,
                             pen_crossbasis = 2, DIC = F, family = "poisson"){

  
  if(family == "poisson"){
    output = DLNM_Laplace_cov_pois(model,
                               crossbasis,
                               ID ,
                               covar.ri,
                               data,
                               map,
                               type,
                               df_inter,
                               df_z,
                               offset ,
                               z, pen_crossbasis , DIC)
  }else if (family =="NB"){
    output = DLNM_Laplace_cov_NB(model,
                               crossbasis,
                               ID ,
                               covar.ri ,
                               data,
                               map,
                               type,
                               df_inter,
                               df_z,
                               offset,
                               z, pen_crossbasis, DIC)
    
  }
}
