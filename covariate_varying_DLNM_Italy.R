library(tidyverse) ; library(data.table); library(splines); library(dlnm); library(sf)

setwd("G:/My Drive/Onderzoek/DLNM/Spatial varying DLNM/Code/Covariate varying")

source('functions/DLNM_Laplace_covariate.R')
source('functions/DLNM_Laplace_covariate_NB.R')
source('functions/DLNM_Laplace_covariate_poisson.R')
source('functions/predRR_covariate.R')
source('functions/help_functions.R')
source('functions/interpret_coefficients.R')
source('functions/af_Laplace.R')

#---------------------------
# Read time series data
#---------------------------

# Read mortality data
mortdata <- fread("application/data/mortality.csv.gz")

# Read temperature series and merge
tempdata <- fread("application/data/tmean.csv.gz")
tsdata <- merge(mortdata, tempdata, all.x = T)

# Order
setkey(tsdata, city_code, date)

# Create date-related variables
tsdata[, ":="(year = year(date), dow = weekdays(date))]

#---------------------------
# Read Metadata
#---------------------------

# Create all city-age combinations
agelabs <- grep("deaths_[[:digit:]]", names(tsdata), value = T) |> 
  gsub(pattern = "deaths_", replacement = "")
metadf <- expand.grid(agegroup = agelabs, 
                      city_code = unique(tsdata$city_code))

# Load data from EUcityTRM and merge
metadata_spatial <- read.csv("application/data/metadata_spatial.csv.gz")
metadf <- merge(metadf, metadata_spatial) 

# Load age-specific demographic data
metadata_age <- read.csv("application/data/metadata_age.csv.gz")




#######################
# Summarize data      #
#######################
mortdata_agg <- mortdata %>%
  mutate(deaths = rowSums(across(contains("death")), na.rm = TRUE))

datafull <- inner_join(inner_join(mortdata_agg, tempdata, by = c("city_code","date")),
                       metadata_spatial, by = "city_code")

datafull[, ":="(year = lubridate::year(date), dow = weekdays(date))]


########################
# Read map of Italy    #
########################
italymap <- st_read("application/data/italymap.shp")
cities_locations <- datafull%>%
  group_by(city_name) %>%
  summarize(lon = min(lon), lat = min(lat), depriv = min(depriv)) %>%
  ungroup()

cities_sf <- st_as_sf(cities_locations, coords = c("lon", "lat"), crs = 4326)

pdf("application/Figures/map_Italy.pdf")
ggplot() +
  geom_sf(data = italymap, fill = "antiquewhite", color = "black") +
  geom_sf(data = cities_sf, aes(color = depriv), size = 2) +
  theme_minimal() +
  labs(title = "Selected Cities in Italy")
dev.off()


################################################################################
#Prepare crossbasis matrix 

library(dlnm)
L <- 21 # maximum lag
vx <- 7 # number of basis for exposure var
vl <- 8 # number of basis for lag var
group <- factor(datafull$city_name)
crossbasis <- crossbasis(datafull$tmean.x, lag=L, # penalized
                         argvar=list(fun="ps",df = vx, intercept=F),
                         arglag=list(fun="ps",df = vl, intercept=T), group=group)


datafull <- datafull %>%
  group_by(city_name) %>%
  mutate(
    quantile_temp = ecdf(tmean.x)(tmean.x)
  ) %>%
  ungroup()
  
crossbasis_quantiles <- crossbasis(datafull$quantile_temp, lag=L, #crossbasis based on quantiles
                                  argvar=list(fun="ps",df = vx, intercept=F),
                                  arglag=list(fun="ps",df = vl, intercept=T), group=group)


y_all = datafull$deaths

datafull$dpv_scaled = scale(datafull$depriv) # scale deprivation index (effect modifier)


# Fit different models
model_laplace_common <- DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + 
                                           ns(dpv_scaled, df = 3, intercept = F) + prop_65p,
                                         crossbasis = crossbasis,
                                         ID = as.factor(datafull$city_name),
                                         covar.ri = "ind",
                                         offset = datafull$pop,
                                         data = datafull, 
                                         type = "none",
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "poisson") 

model_laplace_common_quantile <- DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + 
                                           ns(dpv_scaled, df = 3, intercept = F) + prop_65p,
                                         crossbasis = crossbasis_quantiles,
                                         ID = as.factor(datafull$city_name),
                                         covar.ri = "ind",
                                         offset = datafull$pop,
                                         data = datafull, 
                                         type = "none",
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "poisson") 

model_laplace_common_NB <- DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + 
                                           ns(dpv_scaled, df = 3, intercept = F) + prop_65p,
                                         crossbasis = crossbasis,
                                         ID = as.factor(datafull$city_name),
                                         covar.ri = "ind",
                                         offset = datafull$pop,
                                         data = datafull, 
                                         type = "none",
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "NB") 

model_laplace_common_quantile_NB <- DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + 
                                                    ns(dpv_scaled, df = 3, intercept = F) + prop_65p,
                                                  crossbasis = crossbasis_quantiles,
                                                  ID = as.factor(datafull$city_name),
                                                  covar.ri = "ind",
                                                  offset = datafull$pop,
                                                  data = datafull, 
                                                  type = "none",
                                                  pen_crossbasis = 2,
                                                  DIC = T,
                                                  family = "NB") 

model_laplace_inter <-  DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + prop_65p,
                                   crossbasis = crossbasis,
                                   ID = as.factor(datafull$city_name),
                                   covar.ri = "ind",
                                   offset = datafull$pop,
                                   data = datafull, 
                                   z = as.vector(datafull$dpv_scaled),
                                   type = "non_linear",
                                   df_z = 5,
                                   df_inter = 5,
                                   pen_crossbasis = 2,
                                   DIC = T,
                                   family = "poisson") 

model_laplace_inter_quantile <-  DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + prop_65p,
                                   crossbasis = crossbasis_quantiles,
                                   ID = as.factor(datafull$city_name),
                                   covar.ri = "ind",
                                   offset = datafull$pop,
                                   data = datafull, 
                                   z = as.vector(datafull$dpv_scaled),
                                   type = "non_linear",
                                   df_z = 5,
                                   df_inter = 5,
                                   pen_crossbasis = 2,
                                   DIC = T,
                                   family = "poisson") 


model_laplace_inter_NB <-  DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + prop_65p,
                                         crossbasis = crossbasis,
                                         ID = as.factor(datafull$city_name),
                                         covar.ri = "ind",
                                         offset = datafull$pop,
                                         data = datafull, 
                                         z = as.vector(datafull$dpv_scaled),
                                         type = "non_linear",
                                         df_z = 5,
                                         df_inter = 5,
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "NB") 

model_laplace_inter_quantile_NB <-  DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + prop_65p,
                                                  crossbasis = crossbasis_quantiles,
                                                  ID = as.factor(datafull$city_name),
                                                  covar.ri = "ind",
                                                  offset = datafull$pop,
                                                  data = datafull, 
                                                  z = as.vector(datafull$dpv_scaled),
                                                  type = "non_linear",
                                                  df_z = 5,
                                                  df_inter = 5,
                                                  pen_crossbasis = 2,
                                                  DIC = T,
                                                  family = "NB") 



model_laplace_linear <-  DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + prop_65p,
                                   crossbasis = crossbasis,
                                   ID = as.factor(datafull$city_name),
                                   covar.ri = "ind",
                                   offset = datafull$pop,
                                   data = datafull, 
                                   z = as.vector(datafull$dpv_scaled),
                                   type = "linear",
                                   df_z = 5,
                                   pen_crossbasis = 2,
                                   DIC = T,
                                   family = "poisson") 

model_laplace_linear_quantile <-  DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + prop_65p,
                                          crossbasis = crossbasis_quantiles,
                                          ID = as.factor(datafull$city_name),
                                          covar.ri = "ind",
                                          offset = datafull$pop,
                                          data = datafull, 
                                          z = as.vector(datafull$dpv_scaled),
                                          type = "linear",
                                          df_z = 5,
                                          pen_crossbasis = 2,
                                          DIC = T,
                                          family = "poisson") 


model_laplace_linear_NB <-  DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + prop_65p,
                                          crossbasis = crossbasis,
                                          ID = as.factor(datafull$city_name),
                                          covar.ri = "ind",
                                          offset = datafull$pop,
                                          data = datafull, 
                                          z = as.vector(datafull$dpv_scaled),
                                          type = "linear",
                                          df_z = 5,
                                          pen_crossbasis = 2,
                                          DIC = T,
                                          family = "NB") 

model_laplace_linear_quantile_NB <-  DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year))) + prop_65p,
                                                   crossbasis = crossbasis_quantiles,
                                                   ID = as.factor(datafull$city_name),
                                                   covar.ri = "ind",
                                                   offset = datafull$pop,
                                                   data = datafull, 
                                                   z = as.vector(datafull$dpv_scaled),
                                                   type = "linear",
                                                   df_z = 5,
                                                   pen_crossbasis = 2,
                                                   DIC = T,
                                                   family = "NB") 

dens_binary = matrix(ifelse(datafull$dpv_scaled > median(datafull$dpv_scaled),1,0), ncol = 1)
model_laplace_binary <- DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year)))+ prop_65p ,
                                         crossbasis = crossbasis,
                                         ID = as.factor(datafull$city_name),
                                         covar.ri = "ind",
                                         offset = datafull$pop,
                                         data = datafull, 
                                         z = dens_binary,
                                         type = "factor",
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "poisson") 

model_laplace_binary_quantile <- DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year)))+ prop_65p ,
                                         crossbasis = crossbasis_quantiles,
                                         ID = as.factor(datafull$city_name),
                                         covar.ri = "ind",
                                         offset = datafull$pop,
                                         data = datafull, 
                                         z = dens_binary,
                                         type = "factor",
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "poisson") 

model_laplace_binary_NB <- DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year)))+ prop_65p ,
                                         crossbasis = crossbasis,
                                         ID = as.factor(datafull$city_name),
                                         covar.ri = "ind",
                                         offset = datafull$pop,
                                         data = datafull, 
                                         z = dens_binary,
                                         type = "factor",
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "NB") 

model_laplace_binary_quantile_NB <- DLNM_Laplace_cov(y_all ~  dow + ns(date, df = 7 * length(unique(datafull$year)))+ prop_65p ,
                                                  crossbasis = crossbasis_quantiles,
                                                  ID = as.factor(datafull$city_name),
                                                  covar.ri = "ind",
                                                  offset = datafull$pop,
                                                  data = datafull, 
                                                  z = dens_binary,
                                                  type = "factor",
                                                  pen_crossbasis = 2,
                                                  DIC = T,
                                                  family = "NB") 


#Compare all models based on DIC
c(model_laplace_inter$DIC, model_laplace_linear$DIC, model_laplace_binary$DIC, model_laplace_common$DIC,
  model_laplace_inter_quantile$DIC, model_laplace_linear_quantile$DIC, model_laplace_binary_quantile$DIC,
  model_laplace_common_quantile$DIC,model_laplace_inter_NB$DIC, model_laplace_linear_NB$DIC, 
  model_laplace_binary_NB$DIC, model_laplace_common_NB$DIC,  model_laplace_inter_quantile_NB$DIC,
  model_laplace_linear_quantile_NB$DIC, model_laplace_binary_quantile_NB$DIC, model_laplace_common_quantile_NB$DIC)




################
# Calculate af #
################
data2021 <- datafull%>% filter(year==2021)
group_af = factor(data2021$city_name)
est_af = attrdl_Laplace(data2021$quantile_temp,model_laplace_inter_quantile_NB,data2021$deaths, data2021$dpv_scaled,
                        type="af",dir="back",tot=TRUE,cen = 0.6, sim=TRUE, nsim = 500, ID=group_af) 

est_af_fit = attrdl_Laplace(x = data2021$quantile_temp,model = model_laplace_inter_quantile_NB,cases = data2021$deaths, z = data2021$dpv_scaled,
                            type="af",dir="back",tot=TRUE,cen = 0.6, sim=FALSE, ID=group_af) 

result_af <- data.frame(ID = unique(group_af), AF = est_af_fit, depriv = data2021$depriv[seq(1,dim(data2021)[1], by = 365)],
                        lower = apply(est_af,1,quantile, probs = 0.025),
                        upper = apply(est_af,1,quantile, probs = 0.975)) %>%
  arrange(AF) %>%
  mutate(ID = factor(ID, levels = ID))



pdf("application/Figures/af_smooth.pdf", height = 16, width = 12)
ggplot(result_af, aes(x = AF, y = ID)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = depriv), height = 0.2, size = 1) +
  geom_point(aes(color = depriv), size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  scale_color_gradient(low = "pink", high = "purple") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Estimated attributable fraction",
    y = "City",
    color = "Deprivation"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14) 
  )


dev.off()


# Counterfactual scenarios

# Same deprivation
data2021$same_depriv = median(data2021$dpv_scaled)
est_af_depriv = attrdl_Laplace(data2021$quantile_temp,model_laplace_inter_quantile_NB,data2021$deaths, data2021$same_depriv,
                        type="af",dir="back",tot=TRUE,cen = 0.6, sim=TRUE, nsim = 500, ID=group_af )

est_af_fit_depriv = attrdl_Laplace(x = data2021$quantile_temp,model = model_laplace_inter_quantile_NB,cases = data2021$deaths, z = data2021$same_depriv,
                            type="af",dir="back",tot=TRUE,cen = 0.6, sim=FALSE, ID=group_af)

result_af_depriv <- data.frame(ID = unique(group_af), AF = est_af_fit_depriv,
                        lower = apply(est_af_depriv,1,quantile, probs = 0.025),
                        upper = apply(est_af_depriv,1,quantile, probs = 0.975),
                        depriv = data2021$depriv[seq(1,dim(data2021)[1], by = 365)]) %>%
  arrange(AF) %>%
  mutate(ID = factor(ID, levels = ID))



library(ggnewscale)
pdf("application/Figures/af_smooth_counter.pdf", height = 16, width = 12)

legend_df <- data.frame(
  type = c("True", "Counterfactual"),
  x = c(1, 1),
  y = c(0, 0)
)

ggplot(result_af, aes(x = AF, y = ID)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color =depriv), height = 0.2, size = 1,
                 show.legend = F) +
  geom_point(aes(color = depriv), size = 3,
             show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  scale_color_gradient(low = "lightcoral", high = "darkred") +
  new_scale_color() +
  geom_errorbarh(data = result_af_depriv, aes(xmin = lower, xmax = upper, color = depriv), 
                 height = 0.2, size = 1, alpha = 0.7, show.legend = F) +
  geom_point(data = result_af_depriv, aes(color =depriv), size = 3, shape = 17,
             show.legend = F) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Estimated attributable fraction",
    y = "City",
    color = "AF counterfactual"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14) 
  )+
  #  Dummy points for discrete legend
  new_scale_color() +
  geom_point(data = legend_df, 
             aes(x = x, y = y, color = type, shape = type), size = 4) +
  scale_color_manual(
    name = NULL,
    values = c("True" = "darkred", "Counterfactual" = "darkblue")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("True" = 16, "Counterfactual" = 17)
  ) +
  xlim (c(-0.02,0.1))

dev.off()


# Same temperature
data2021 <- data2021 %>%
  group_by(date) %>%
 mutate(quantile_common = mean(quantile_temp)) %>%
  ungroup()
est_af_temp = attrdl_Laplace(data2021$quantile_common,model_laplace_inter_quantile_NB,data2021$deaths, data2021$same_depriv,
                               type="af",dir="back",tot=TRUE,cen = 0.6, sim=TRUE, nsim = 500, ID=group_af) 

est_af_fit_temp = attrdl_Laplace(x = data2021$quantile_common,model = model_laplace_inter_quantile_NB,cases = data2021$deaths, z = data2021$same_depriv,
                                   type="af",dir="back",tot=TRUE,cen = 0.6, sim=FALSE, ID=group_af) 

result_af_temp <- data.frame(ID = unique(group_af), est = est_af_fit_temp,
                               lower = apply(est_af_temp,1,quantile, probs = 0.025),
                               upper = apply(est_af_temp,1,quantile, probs = 0.975)) %>%
  arrange(est) %>%
  mutate(ID = factor(ID, levels = ID))




######################
# Make predictions   #
######################
library(ggplot2)
library(plotly)
z_mean = mean(datafull$depriv)
z_sd = sd(datafull$depriv)

# Smooth Laplace quantiles
at_x_quantile = seq(0,1, length = 100)


# RR versus temperature for different deprivation quantiles

pred1_quantile = predRR(model_laplace_inter_quantile_NB, at_x_quantile, cen = 0.6, L = 21, quantile(datafull$dpv_scaled, 0.025), by = 0.5)
pred2_quantile = predRR(model_laplace_inter_quantile_NB, at_x_quantile, cen = 0.6, L = 21, quantile(datafull$dpv_scaled, 0.5), by = 0.5)
pred3_quantile = predRR(model_laplace_inter_quantile_NB, at_x_quantile, cen = 0.6, L = 21, quantile(datafull$dpv_scaled, 0.975), by = 0.5)


col <- c("darkgoldenrod3", "aquamarine3", "darkred")

pdf("application/Figures/RR_by_cov_smooth.pdf", width = 6, height = 6)
parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
# Plot gam
plot(at_x_quantile , pred3_quantile $pred_all, type = "l", ylim=c(0.8,1.5), ylab="RR", col=col[1], lwd=1.5,
     xlab="Temperature percentile")
lines(at_x_quantile , pred2_quantile $pred_all, col = col[2])
lines(at_x_quantile , pred1_quantile $pred_all, col = col[3])


abline(v = 0.025, lty = 2)
abline(v = 0.5, lty = 2)
abline(v = 0.975, lty = 2)
abline(h = 1, lty = 2)

legend("top", c("97.5% deprivation", "50% deprivation", "2.5% deprivation"), lty=1, lwd=1.5, col=col, bty="n",
       inset=0.05, y.intersp=2, cex=0.8)


fci <- function(x, high, low, ci.arg, plot.arg, noeff = NULL){
  polygon.arg <- modifyList(list(col = grey(0.9), border = NA), 
                            ci.arg)
  polygon.arg <- modifyList(polygon.arg, list(x = c(x, 
                                                    rev(x)), y = c(high, rev(low))))
  do.call(polygon, polygon.arg)
}

plot.arg1 <- list(type = "l",col=col[1],  lwd=1.5)
fci(x=at_x_quantile , high = pred3_quantile$Qupper_all,
    low = pred3_quantile$Qlower_all, ci.arg=list(col=alpha(col[1], 0.2)), plot.arg = plot.arg1)
plot.arg2 <- list(type = "l",col=col[2],  lwd=1.5)
fci(x=at_x_quantile , high = pred2_quantile$Qupper_all,
    low = pred2_quantile$Qlower_all, ci.arg=list(col=alpha(col[2], 0.2)), plot.arg = plot.arg2)
plot.arg3 <- list(type = "l",col=col[3],  lwd=1.5)
fci(x=at_x_quantile , high = pred1_quantile$Qupper_all,
    low = pred1_quantile$Qlower_all, ci.arg=list(col=alpha(col[3], 0.2)), plot.arg = plot.arg3)
dev.off()


z_pred = seq(quantile(datafull$dpv_scaled,0.025), quantile(datafull$dpv_scaled,0.975), length = 100)
by.z_pred_quantile = by.z_lower_quantile = by.z_upper_quantile = matrix(0, nrow = 3, ncol = length(z_pred))

for(pred_i in 1:length(z_pred)){
  by.z_all_quantile = predRR(model_laplace_inter_quantile_NB, c(0.025,0.5,0.975), 0.6,L,at_inter = z_pred[pred_i])
  by.z_pred_quantile[,pred_i] = by.z_all_quantile$pred_all
  by.z_lower_quantile[,pred_i] = by.z_all_quantile$Qlower_all
  by.z_upper_quantile[,pred_i] = by.z_all_quantile$Qupper_all
}



# RR versus deprivation for different temperature quantiles

pdf("application/Figures/RR_by_temp_smooth.pdf", width = 6, height = 6)
col <- c("darkgoldenrod3", "aquamarine3", "darkred")
parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
# Plot gam
plot(z_pred*z_sd+z_mean, by.z_pred_quantile[3,], type = "l", ylim=c(0.8,1.5), ylab="RR", col=col[1], lwd=1.5,
     xlab="Deprivation")
lines(z_pred*z_sd+z_mean, by.z_pred_quantile[2,], col = col[2])
lines(z_pred*z_sd+z_mean, by.z_pred_quantile[1,], col = col[3])


abline(v = quantile(datafull$depriv,0.025), lty = 2)
abline(v = quantile(datafull$depriv,0.5), lty = 2)
abline(v = quantile(datafull$depriv,0.975), lty = 2)
abline(h = 1, lty = 2)


legend("top", c("97.5% temperature", "50% temperature", "2.5% temperature"), lty=1, lwd=1.5, col=col, bty="n",
       inset=0.05, y.intersp=2, cex=0.8)


fci <- function(x, high, low, ci.arg, plot.arg, noeff = NULL){
  polygon.arg <- modifyList(list(col = grey(0.9), border = NA), 
                            ci.arg)
  polygon.arg <- modifyList(polygon.arg, list(x = c(x, 
                                                    rev(x)), y = c(high, rev(low))))
  do.call(polygon, polygon.arg)
}

plot.arg1 <- list(type = "l",col=col[1],  lwd=1.5)
fci(x=z_pred*z_sd+z_mean  , high = by.z_upper_quantile[3,],
    low = by.z_lower_quantile[3,], ci.arg=list(col=alpha(col[1], 0.2)), plot.arg = plot.arg1)
plot.arg2 <- list(type = "l",col=col[2],  lwd=1.5)
fci(x=z_pred*z_sd+z_mean  , high = by.z_upper_quantile[2,],
    low = by.z_lower_quantile[2,], ci.arg=list(col=alpha(col[2], 0.2)), plot.arg = plot.arg2)
plot.arg3 <- list(type = "l",col=col[3],  lwd=1.5)
fci(x=z_pred*z_sd+z_mean  , high = by.z_upper_quantile[1,],
    low = by.z_lower_quantile[1,], ci.arg=list(col=alpha(col[3], 0.2)), plot.arg = plot.arg3)
dev.off()


# Lag specific relationship for different deprivation quantiles

pred1_quantile_lag = predRR(model_laplace_inter_quantile_NB, 0.95, cen = 0.6, L = 21, quantile(datafull$dpv_scaled, 0.025), by = 1)
pred2_quantile_lag = predRR(model_laplace_inter_quantile_NB, 0.95, cen = 0.6, L = 21, quantile(datafull$dpv_scaled, 0.5), by = 1)
pred3_quantile_lag = predRR(model_laplace_inter_quantile_NB, 0.95, cen = 0.6, L = 21, quantile(datafull$dpv_scaled, 0.975), by = 1)


col <- c("darkgoldenrod3", "aquamarine3", "darkred")

pdf("application/Figures/RR_by_lag_smooth.pdf", width = 6, height = 6)
parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
# Plot gam
plot(0:21 , exp(pred3_quantile_lag$logpredX), type = "l", ylim=c(0.9,1.2), ylab="RR", col=col[1], lwd=1.5,
     xlab="Temperature percentile")
lines(0:21 , exp(pred2_quantile_lag$logpredX), col = col[2])
lines(0:21 , exp(pred1_quantile_lag$logpredX), col = col[3])

legend("top", c("97.5% deprivation", "50% deprivation", "2.5% deprivation"), lty=1, lwd=1.5, col=col, bty="n",
       inset=0.05, y.intersp=2, cex=0.8)


fci <- function(x, high, low, ci.arg, plot.arg, noeff = NULL){
  polygon.arg <- modifyList(list(col = grey(0.9), border = NA), 
                            ci.arg)
  polygon.arg <- modifyList(polygon.arg, list(x = c(x, 
                                                    rev(x)), y = c(high, rev(low))))
  do.call(polygon, polygon.arg)
}

plot.arg1 <- list(type = "l",col=col[1],  lwd=1.5)
fci(x=0:21 , high = exp(pred3_quantile_lag$Qupper_logpredX),
    low = exp(pred3_quantile_lag$Qlower_logpredX), ci.arg=list(col=alpha(col[1], 0.2)), plot.arg = plot.arg1)
plot.arg2 <- list(type = "l",col=col[2],  lwd=1.5)
fci(x=0:21 , high = exp(pred2_quantile_lag$Qupper_logpredX),
    low = exp(pred2_quantile_lag$Qlower_logpredX), ci.arg=list(col=alpha(col[2], 0.2)), plot.arg = plot.arg2)
plot.arg3 <- list(type = "l",col=col[3],  lwd=1.5)
fci(x=0:21 , high = exp(pred1_quantile_lag$Qupper_logpredX),
    low = exp(pred1_quantile_lag$Qlower_logpredX), ci.arg=list(col=alpha(col[3], 0.2)), plot.arg = plot.arg3)
dev.off()



# Overall RR for all areas
ind_area = 1
pred_overall_area = lower_overall_area = upper_overall_area = matrix(0, nrow = length(at_x_quantile), ncol = length(unique(datafull$depriv)))
for(pred_i in unique(datafull$dpv_scaled)){
  pred_all_area = predRR(model_laplace_inter_quantile_NB, at_x_quantile, 0.6,L,at_inter = pred_i)
  pred_overall_area[,ind_area] = pred_all_area$pred_all
  lower_overall_area[,ind_area] = pred_all_area$Qlower_all
  upper_overall_area[,ind_area] = pred_all_area$Qupper_all
  
  ind_area = ind_area+1
}


df_plot_area <- as.data.frame(pred_overall_area) %>%
  mutate(at_x = at_x_quantile) %>%
  pivot_longer(cols = -at_x,
               names_to = "variable",
               values_to = "RR") %>%
  mutate(depriv = rep(unique(datafull$depriv),100))

pdf("application/Figures/all_area.pdf", width = 6, height = 4)
ggplot(df_plot_area, aes(x = at_x, y = RR, group = variable, col = depriv)) +
  geom_line()+
  theme_minimal() +
  scale_color_gradientn(colors = c("#bdbdbd", "#ef3b2c", "#99000d"))  +
  labs(
    x = "Temperature percentile",
    y = "RR"
  )
dev.off()


# Grid of temperature and deprivation


smooth_pred_quantile = smooth_lower_quantile = smooth_upper_quantile = matrix(0, nrow = length(at_x_quantile), ncol = length(z_pred))

for(pred_i in 1:length(z_pred)){
  pred_smooth_all_quantile = predRR(model_laplace_inter_quantile_NB, at_x_quantile, 0.6,L,at_inter = z_pred[pred_i])
  smooth_pred_quantile[,pred_i] = pred_smooth_all_quantile$pred_all
  smooth_lower_quantile[,pred_i] = pred_smooth_all_quantile$Qlower_all
  smooth_upper_quantile[,pred_i] = pred_smooth_all_quantile$Qupper_all
}


df_plot_smooth_quantile <- expand.grid(at_x = at_x_quantile, z = z_pred*z_sd+z_mean) %>%
  mutate(
    RR = as.vector(smooth_pred_quantile),
    lower = as.vector(smooth_lower_quantile),
    upper = as.vector(smooth_upper_quantile)
  )



# Heat map

ggplot(df_plot_smooth_quantile, aes(x = at_x, y = z, fill = RR)) + # %>% mutate(RR = ifelse(RR>2.5,2.5,RR))
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()  +
  labs(x = "Temperature percentile", y = "Deprivation", fill = "RR")


heat_map1 <- ggplot(df_plot_smooth_quantile  %>% filter(at_x < 0.25), aes(x = at_x, y = z, fill = RR)) + # %>% mutate(RR = ifelse(RR>2.5,2.5,RR))
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(x = "Temperature percentile", y = "Deprivation", fill = "RR")

heat_map2 <- ggplot(df_plot_smooth_quantile  %>% filter(at_x > 0.75), aes(x = at_x, y = z, fill = RR)) + # %>% mutate(RR = ifelse(RR>2.5,2.5,RR))
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(x = "Temperature percentile", y = "Deprivation", fill = "RR")

library(gridExtra)
pdf("application/Figures/contour_smooth.pdf", width = 6, height = 4)
grid.arrange(heat_map1, heat_map2, ncol = 2)
dev.off()

# Interactive plot with slider for continuous z
fig_quantile <- plot_ly(
  df_plot_smooth_quantile,
  x = ~at_x,
  y = ~RR,
  frame = ~z,
  type = 'scatter',
  mode = 'lines',
  name = "RR"
) %>%
  add_ribbons(ymin = ~lower, ymax = ~upper, line = list(color = 'transparent'),
              fillcolor = 'rgba(0,100,80,0.2)', name = "95% CI")%>%
  layout(
    title = "RR vs temperature",
    xaxis = list(title = "temperature percentile"),
    yaxis = list(title = "RR")
  ) %>%
  animation_opts(frame = 0, redraw = TRUE) %>%
  animation_slider(currentvalue = list(prefix = "Deprivation = "))

fig_quantile

htmlwidgets::saveWidget(fig_quantile, "application/Figures/smooth_RR.html", selfcontained = TRUE)


# 3D plot
pdf("application/Figures/RR_3d_smooth.pdf", height = 6, width = 6)
plot_ly(
  x = ~z_pred*z_sd+z_mean,       
  y = ~at_x_quantile,    
  z = ~smooth_pred_quantile,  
  type = "surface",
  colorbar = list(title = "RR")
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Deprivation"),
      yaxis = list(title = "Temp Percentile"),
      zaxis = list(title = "RR")
    )
  )
dev.off()


# Red zone
pred_exc_smooth_quantile = matrix(0, nrow = length(at_x_quantile), ncol = length(z_pred))
for(pred_i in 1:length(z_pred)){
  pred_exc = predRR(model_laplace_inter_quantile_NB, at_x_quantile, 0.6,L,at_inter = z_pred[pred_i],
                                    exc.prob = T)
  pred_exc_smooth_quantile[,pred_i] = pred_exc$exc.prob
}



# Heat map
df_plot_smooth_heat <- expand.grid(at_x = at_x_quantile, z = z_pred*z_sd+z_mean) %>%
  mutate(
    exc = as.vector(pred_exc_smooth_quantile),
  )

df_plot_smooth_heat_low <- df_plot_smooth_heat %>% filter(at_x < 0.25)
df_plot_smooth_heat_high <- df_plot_smooth_heat %>% filter(at_x > 0.75)

pdf("application/Figures/exceedance_low_smooth.pdf", height = 4, width = 6)
ggplot(df_plot_smooth_heat_low, aes(x = at_x, y = z, fill = exc)) + # %>% mutate(RR = ifelse(RR>2.5,2.5,RR))
  geom_tile() +
  geom_hline(aes(yintercept=quantile(datafull$depriv,0.5)), linetype = "dashed")+
  scale_fill_gradient(low = "mistyrose", high = "darkred")+
  theme_minimal() +
  labs(x = "Temperature percentile", y = "Deprivation", fill = "exceedance")
dev.off()

pdf("application/Figures/exceedance_high_smooth.pdf", height = 4, width = 6)
ggplot(df_plot_smooth_heat_high, aes(x = at_x, y = z, fill = exc)) + # %>% mutate(RR = ifelse(RR>2.5,2.5,RR))
  geom_tile() +
  geom_hline(aes(yintercept=quantile(datafull$depriv,0.5)), linetype = "dashed")+
  scale_fill_gradient(low = "mistyrose", high = "darkred")+
  theme_minimal() +
  labs(x = "Temperature percentile", y = "Deprivation", fill = "exceedance")
dev.off()



# Independent random effects
cities_locations$random_effect = model_laplace_inter_quantile_NB$xispat
cities_sf <- st_as_sf(cities_locations, coords = c("lon", "lat"), crs = 4326)

pdf("application/Figures/map_Italy_re.pdf")
ggplot() +
  geom_sf(data = italymap, fill = "antiquewhite", color = "black") +
  geom_sf(data = cities_sf, aes(color = random_effect), size = 2) +
  theme_minimal() +
  labs(title = "Random effect of selected Cities in Italy")
dev.off()





# RRR
z_pred_RRR = quantile(datafull$dpv_scaled,0.9)
at_x_quantile_RRR = seq(0.025,0.975, length = 40)
pred_quantile_RRR = lower_quantile_RRR = upper_quantile_RRR = matrix(0, nrow = 40, ncol = length(z_pred_RRR))

for(pred_i in 1:length(at_x_quantile_RRR)){
  by.z_all_quantile_RRR = InterpretRR(model_laplace_inter_quantile_NB, at_x_quantile_RRR[pred_i], 0.6,0:L,at_inter = z_pred_RRR,
                                  ref_inter=quantile(datafull$dpv_scaled,0.1))
  pred_quantile_RRR[pred_i,] = by.z_all_quantile_RRR$log_RR_change
  lower_quantile_RRR[pred_i,] = by.z_all_quantile_RRR$log_lower
  upper_quantile_RRR[pred_i,] = by.z_all_quantile_RRR$log_upper
}


df_plot_RRR <- expand.grid(at_x = at_x_quantile_RRR, z = z_pred_RRR*z_sd+z_mean) %>%
  mutate(
    RRR = exp(as.vector(pred_quantile_RRR)),
    lower_RRR = exp(as.vector(lower_quantile_RRR)),
    upper_RRR = exp(as.vector(upper_quantile_RRR))
  ) %>%
  arrange(z)

pdf("application/Figures/RRR.pdf", height = 4, width = 6)
ggplot(df_plot_RRR, aes(x = at_x, y = RRR)) +
  geom_errorbar(aes(ymin = lower_RRR, ymax = upper_RRR), width = 0.01, size = 0.5) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    x = "temperature percentile",
    y = "RRR",
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14) 
  )

dev.off()

