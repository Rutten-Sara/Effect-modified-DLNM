library(tidyverse) ; library(data.table); library(splines); library(dlnm); library(sf)

setwd("G:/My Drive/Onderzoek/DLNM/Spatial varying DLNM/Code/Covariate varying")

source('functions/DLNM_Laplace_covariate.R')
source('functions/DLNM_Laplace_covariate_NB.R')
source('functions/DLNM_Laplace_covariate_poisson.R')
source('functions/predRR_covariate.R')
source('functions/help_functions.R')
source('functions/interpret_coefficients.R')
source('functions/af_Laplace.R')


# Data
italymap <- st_read("G:/My Drive/Onderzoek/DLNM/Spatial varying DLNM/Code/Unstructured spatially varying/Results Italy/data/Limiti2021/Com2021/Com2021.shp")
sicilymap <- italymap %>% filter(COD_REG == 19) %>%
  rename(COD_PROVCOM = PRO_COM) %>%
  arrange(COD_PROVCOM)


library(readxl)
full_data <- read_excel("C:/Users/lucp12805/Downloads/Composite fragility index - all municipalities (IT1,DF_COMP_FRA_IND_MUNICIPAL_01,1.0).xlsx")


full_data <- (full_data %>% rename("COMUNE" = "Indicator"))[ -7595,]

sicily_ind_data <- full_data %>% filter(COMUNE %in% sicilymap$COMUNE) %>%
  select(COMUNE, `Population dependency index`) %>%
  inner_join(sicilymap, by = "COMUNE")



data_temp <- read.csv("G:/My Drive/Onderzoek/DLNM/Spatial varying DLNM/Code/Unstructured spatially varying/Results Italy/data/Sicilia.csv", stringsAsFactors = FALSE)%>%
  arrange(COD_PROVCOM, date) %>%
  inner_join(sicily_ind_data[,c("COD_PROVCOM","POP21", "Population dependency index")], by = "COD_PROVCOM") 

sicilymap <- sicilymap %>% filter(COD_PROVCOM %in% unique(data_temp$COD_PROVCOM))


data_temp <- data_temp %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"))%>%
  mutate(dow = weekdays(date)) %>%
  mutate(dpv_scaled = scale(data_temp$`Population dependency index`))



# Fit model

library(dlnm)
L <- 21 # maximum lag
vx <- 7 # number of basis for exposure var
vl <- 8 # number of basis for lag var
group <- factor(data_temp$COD_PROVCOM)
crossbasis <- crossbasis(data_temp$temperature, lag=L, # penalized
                         argvar=list(fun="ps",df = vx, intercept=F),
                         arglag=list(fun="ps",df = vl, intercept=T), group=group)



y_all = data_temp$dtot



# Fit different models

model_laplace_common <- DLNM_Laplace_cov(y_all ~  dow,
                                         crossbasis = crossbasis,
                                         ID = as.factor(data_temp$COD_PROVCOM),
                                         covar.ri = "Leroux",
                                         smooth = data_temp$date,
                                         map= sicilymap,
                                         df_smooth = 7*length(unique(data_temp$year)),
                                         offset = data_temp$POP21,
                                         data = data_temp, 
                                         type = "none",
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "poisson") 




model_laplace_common_NB <- DLNM_Laplace_cov(y_all ~  dow ,
                                            crossbasis = crossbasis,
                                            ID = as.factor(data_temp$COD_PROVCOM),
                                            covar.ri = "Leroux",
                                            smooth = data_temp$date,
                                            map= sicilymap,
                                            df_smooth = 7*length(unique(data_temp$year)),
                                            offset = data_temp$POP21,
                                            data = data_temp, 
                                            type = "none",
                                            pen_crossbasis = 2,
                                            DIC = T,
                                            family = "NB") 




model_laplace_inter <-  DLNM_Laplace_cov(y_all ~  dow ,
                                         crossbasis = crossbasis,
                                         ID = as.factor(data_temp$COD_PROVCOM),
                                         covar.ri = "Leroux",
                                         smooth = data_temp$date,
                                         map= sicilymap,
                                         df_smooth = 7*length(unique(data_temp$year)),
                                         offset = data_temp$POP21,
                                         data = data_temp, 
                                         z = as.vector(data_temp$dpv_scaled),
                                         type = "non_linear",
                                         df_z = 10,
                                         df_inter = 5,
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "poisson") 



model_laplace_inter_NB <-  DLNM_Laplace_cov(y_all ~  dow,
                                            crossbasis = crossbasis,
                                            ID = as.factor(data_temp$COD_PROVCOM),
                                            covar.ri = "Leroux",
                                            smooth = data_temp$date,
                                            map= sicilymap,
                                            df_smooth = 7*length(unique(data_temp$year)),
                                            offset = data_temp$POP21,
                                            data = data_temp, 
                                            z = as.vector(data_temp$dpv_scaled),
                                            type = "non_linear",
                                            df_z = 10,
                                            df_inter = 5,
                                            pen_crossbasis = 2,
                                            DIC = T,
                                            family = "NB") 


model_laplace_linear <-  DLNM_Laplace_cov(y_all ~  dow,
                                          crossbasis = crossbasis,
                                          ID = as.factor(data_temp$COD_PROVCOM),
                                          covar.ri = "Leroux",
                                          smooth = data_temp$date,
                                          map= sicilymap,
                                          df_smooth = 7*length(unique(data_temp$year)),
                                          offset = data_temp$POP21,
                                          data = data_temp, 
                                          z = as.vector(data_temp$dpv_scaled),
                                          type = "linear",
                                          df_z = 10,
                                          pen_crossbasis = 2,
                                          DIC = T,
                                          family = "poisson") 



model_laplace_linear_NB <-  DLNM_Laplace_cov(y_all ~  dow ,
                                             crossbasis = crossbasis,
                                             ID = as.factor(data_temp$COD_PROVCOM),
                                             covar.ri = "Leroux",
                                             smooth = data_temp$date,
                                             map= sicilymap,
                                             df_smooth = 7*length(unique(data_temp$year)),
                                             offset = data_temp$POP21,
                                             data = data_temp, 
                                             z = as.vector(data_temp$dpv_scaled),
                                             type = "linear",
                                             df_z = 10,
                                             pen_crossbasis = 2,
                                             DIC = T,
                                             family = "NB") 




dens_binary = matrix(ifelse(data_temp$dpv_scaled > median(data_temp$dpv_scaled),1,0), ncol = 1)
model_laplace_binary <- DLNM_Laplace_cov(y_all ~  dow ,
                                         crossbasis = crossbasis,
                                         ID = as.factor(data_temp$COD_PROVCOM),
                                         covar.ri = "Leroux",
                                         smooth = data_temp$date,
                                         map= sicilymap,
                                         df_smooth = 7*length(unique(data_temp$year)),
                                         offset = data_temp$POP21,
                                         data = data_temp, 
                                         z = dens_binary,
                                         type = "factor",
                                         pen_crossbasis = 2,
                                         DIC = T,
                                         family = "poisson") 




model_laplace_binary_NB <- DLNM_Laplace_cov(y_all ~  dow,
                                            crossbasis = crossbasis,
                                            ID = as.factor(data_temp$COD_PROVCOM),
                                            covar.ri = "Leroux",
                                            smooth = data_temp$date,
                                            map= sicilymap,
                                            df_smooth = 7*length(unique(data_temp$year)),
                                            offset = data_temp$POP21,
                                            data = data_temp, 
                                            z = dens_binary,
                                            type = "factor",
                                            pen_crossbasis = 2,
                                            DIC = T,
                                            family = "NB") 




#Compare all models based on DIC
c(model_laplace_inter$DIC, model_laplace_linear$DIC, model_laplace_binary$DIC, model_laplace_common$DIC,
    model_laplace_inter_NB$DIC, model_laplace_linear_NB$DIC, model_laplace_binary_NB$DIC, model_laplace_common_NB$DIC)-
  min(c(model_laplace_inter$DIC, model_laplace_linear$DIC, model_laplace_binary$DIC, model_laplace_common$DIC,
       model_laplace_inter_NB$DIC, model_laplace_linear_NB$DIC,model_laplace_binary_NB$DIC, model_laplace_common_NB$DIC))



######################
# Make predictions   #
######################
library(ggplot2)
library(plotly)
z_mean = mean(data_temp$depriv)
z_sd = sd(data_temp$depriv)

# Smooth Laplace quantiles
at_x_quantile = seq(10,28, by = 0.5)


# RR versus temperature for different deprivation quantiles

pred1_quantile = predRR(model_laplace_inter_NB, at_x_quantile, cen = 20, L = 21, quantile(data_temp$dpv_scaled, 0.025), by = 0.5)
pred2_quantile = predRR(model_laplace_inter_NB, at_x_quantile, cen = 20, L = 21, quantile(data_temp$dpv_scaled, 0.5), by = 0.5)
pred3_quantile = predRR(model_laplace_inter_NB, at_x_quantile, cen = 20, L = 21, quantile(data_temp$dpv_scaled, 0.975), by = 0.5)


col <- c("darkgoldenrod3", "aquamarine3", "darkred")

parold <- par(no.readonly=T)
par(mar=c(4,4,1,0.5), las=1, mgp=c(2.5,1,0))
# Plot gam
plot(at_x_quantile , pred3_quantile $pred_all, type = "l", ylim=c(0.8,1.6), ylab="RR", col=col[1], lwd=1.5,
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




# Lag specific relationship for different deprivation quantiles

pred1_quantile_lag = predRR(model_laplace_inter_NB, 28, cen = 20, L = 21, quantile(data_temp$dpv_scaled, 0.025), by = 1)
pred2_quantile_lag = predRR(model_laplace_inter_NB, 28, cen = 20, L = 21, quantile(data_temp$dpv_scaled, 0.5), by = 1)
pred3_quantile_lag = predRR(model_laplace_inter_NB, 28, cen = 20, L = 21, quantile(data_temp$dpv_scaled, 0.975), by = 1)


col <- c("darkgoldenrod3", "aquamarine3", "darkred")

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







# RRR
z_pred_RRR = quantile(data_temp$dpv_scaled,0.9)
at_x_quantile_RRR = at_x_quantile
pred_quantile_RRR = lower_quantile_RRR = upper_quantile_RRR = matrix(0, nrow = 37, ncol = length(z_pred_RRR))

for(pred_i in 1:length(at_x_quantile_RRR)){
  by.z_all_quantile_RRR = InterpretRR(model_laplace_inter_NB, at_x_quantile_RRR[pred_i], 20,0:L,at_inter = z_pred_RRR,
                                      ref_inter=quantile(data_temp$dpv_scaled,0.1))
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




