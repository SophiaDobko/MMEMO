######################################################################+
# Externally calibrate Hydrus-1D using Particle Swarm Optimization ####
######################################################################+

library(hydrusR)
library(xts)
library(dplyr)
library(soilwater)
library(tidyverse)
library(lhs)

# Directories and basic information ####
# Adapt information in this first chapter 

source("D:/HydrusMMEMO/MMEMO/hydrus/hydpara_selector.R")

# define wd and other directories ####
comp <- Sys.getenv("COMPUTERNAME")

if(comp == "GK-NB-5"){
  setwd("D:/CosmicSense/Hydrus")
  hydrus_path <- "C:/hydrus"
  theta_observed <- "C:/Users/dobkowitz/data/crns/merged/rhi_pr2_FT.txt" # observed pr2 soil moisture timeseries, has columns datetime, theta_1, theta_2, theta_3, theta_4
  }

# Choose pr2 profile
pr_plot <- "f50" # or "center", "t25", "f25", "t50", "f50"

# Basic inputs
project_name = "Hydrus_rhi_pr_F50"
if (pr_plot == "t50"){
  project_name = "Hydrus_rhi_pr_T50"
}
project_path = paste0(getwd(), "/", project_name)
TimeUnit = "days"
SpaceUnit = "cm"
PrintTimes = 1
tz = "etc/GMT-1"
start_date = as.Date("2021-04-18")
end_date = as.POSIXct("2021-11-24")
nbefore = 107 # days in the year before first boundary conditions day

para <- data.frame(thr = c(0, 0, 0, 0), # initial soil hydraulic parameters for all materials
                   ths = c(0.797, 0.73, 0.741, 0.891),
                   Alpha = c(0.02, 0.012, 0.005, 0.003),
                   n = c(1.23, 1.12, 1.15, 1.16),
                   Ks = c(14, 6.69, 0.92, 104),
                   l = c(0.5, 0.5, 0.5, 0.5))

# initial estimates of parameters which will be calibrated in pso_
para_initial = cbind(c(para$ths, para$Alpha, para$n, para$Ks)) 


## Observed data ####
# adapt if observed soil moisture timeseries has not/more then 4 depths

th_obs <- read.table(theta_observed, sep = "\t", header = T)
th_obs <- th_obs[th_obs$plot==pr_plot,]
th_obs$date <- as.Date(th_obs$datetime)
th_obs <- th_obs[th_obs$date >= start_date & th_obs$date <= end_date,]
th_obs <- aggregate(x = th_obs[,c("theta_1","theta_2", "theta_3", "theta_4")], by = list(th_obs$date), FUN = mean)
colnames(th_obs)[1] <- "date"
th_obs_wide <- th_obs
th_obs$date <- (1+nbefore):(nrow(th_obs)+nbefore)
th_obs <- gather(th_obs, key = "probe", value = "theta","theta_1", "theta_2", "theta_3", "theta_4")

## Objective function ####
# adapt soil_para due to number of soil materials and parameters to be calibrated
# adapt if modeled soil moisture timeseries has not/more then 4 depths


zielfun = function(params) { # params is produced by ppso, dont worry about this
  soil_para <- data.frame("thr" = para$thr,
                     ths = c(params[1:4]),
                     Alpha = c(params[5:8]),
                     n = c(params[9:12]),
                     Ks = c(params[13:16]),
                     "l" = para$l)
  
  # update hydraulic parameters in selector.in for next run
  hydpara_selector(project.path = project_path, para = soil_para)
  
  # run hydrus
  call.H1D(project.path = project_path,
           hydrus.path = hydrus_path,
           show.output = TRUE)

  #wenn Fehler in Hydrus Berechnung, setze GOF=Inf
  if (file.exists(file.path(project_path, "Error.msg"))) { 
    file.remove(file.path(project_path, "Error.msg"))
    GOF = Inf
  }
  else {
    # read simulated soil moisture
    th_mod <- read.table(paste0(project_path, "/Obs_Node.out"), skip = 10, header = T, fill = T)
    th_mod <- th_mod[-nrow(th_mod),c(1,3,6,9,12)]
    colnames(th_mod) <- c("date","theta_1", "theta_2", "theta_3", "theta_4")
    th_mod$date <- as.numeric(th_mod$date)
    th_mod$date <- round(th_mod$date, digits = 0)
    th_mod <- aggregate(th_mod, by = list(th_mod$date), FUN = mean)
    th_mod <- th_mod[th_mod$date>nbefore,-1]
    th_mod_wide <- th_mod
    th_mod <- gather(th_mod, key = "probe", value = "theta","theta_1", "theta_2", "theta_3", "theta_4")
    
    
    #merge observed and modelled data into 1 dataframe
    df <- data.frame("obs" = th_obs$theta, "mod" = th_mod$theta)
    
    # Calculate Goodness Of Fit
    GOF = sum((df[,1]-df[,2])^2)/sum((df[,1]-mean(df[,1]))^2) #nash-sutcliff, without 1-..., so that minimum ist best
    #GOF = NSE(df[,1], df[,2]) # with this option optim = max_search, or GOF=-1*(NSE(...)-1)
  }
  return(GOF)
}


## Optimization ###########
####################
library(ppso)

start_time = Sys.time()
res=optim_pso(zielfun,number_of_parameters=16, 
              parameter_bounds = cbind(c(rep(0.7,4), rep(0.002,4), rep(1,4), rep(0.5,4)), #lower limit
                                       c(rep(0.9,4), rep(0.03,4), rep(1.3,4), rep(110,4))), #upper limit
              max_number_of_iterations=200, 
              logfile=paste0(project_path, "/ppso.log"), 
              projectfile=paste0(project_path, "/ppso.pro"),
              max_number_function_calls = 5000,
              initial_estimates = para_initial)

# look at results ####
res$value
res$function_calls
res$break_flag
GOF <- res$value
GOF

# save results ####
save(res, file = paste0(project_path, "/result_pso.RData"))

# read results ####
load(paste0(project_path, "/result_pso.RData"))
ppso <- read.table(paste0(project_path, "/ppso.log"), sep = "\t", header = T)
ppso$time <- as.POSIXct(ppso$time)

# choose best run ####
# if calibration has stopped because max number of function calls was reached but calibration has not converged 
# and you dont want to do more calibration runs

if (res$break_flag == "max number of function calls reached"){
  best <- ppso[ppso$objective_function == min(ppso$objective_function),]
  params <- as.numeric(best[1,2:17])
  GOF <- zielfun(params)
  }

th_mod <- read.table(paste0(project_path, "/Obs_Node.out"), skip = 10, header = T, fill = T)
th_mod <- th_mod[-nrow(th_mod),c(1,3,6,9,12)]
colnames(th_mod) <- c("date","theta_1", "theta_2", "theta_3", "theta_4")
th_mod$date <- as.numeric(th_mod$date)
th_mod$date <- round(th_mod$date, digits = 0)
th_mod <- aggregate(th_mod, by = list(th_mod$date), FUN = mean)
th_mod <- th_mod[th_mod$date>nbefore,-1]
th_mod_wide <- th_mod
th_mod <- gather(th_mod, key = "probe", value = "theta","theta_1", "theta_2", "theta_3", "theta_4")
#merge observed and modelled data into 1 dataframe
df <- data.frame("obs" = th_obs$theta, "mod" = th_mod$theta)

# Plot observed and modelled soil moisture ####
th_mod_wide$date <- th_obs_wide$date

png(paste0(project_path, "/pr_F50_obs_mod.png"), height=5, width=7, units = "in", res=500)
par(mfrow = c(2,1), mar = c(0,4,0,2), oma = c(4,1,2,0), las = 1)
plot(th_obs_wide$date, th_obs_wide$theta_1, type = "l", col = 1, ylim = c(0,1),
     ylab = "Soil moisture (m3/m3)", xaxt = "n")
lines(th_obs_wide$date, th_obs_wide$theta_2, type = "l", col = 2)
lines(th_obs_wide$date, th_obs_wide$theta_3, type = "l", col = 3)
lines(th_obs_wide$date, th_obs_wide$theta_4, type = "l", col = 4)

plot(th_mod_wide$date, th_mod_wide$theta_1, type = "l", col = 1, ylim = c(0,1),
     ylab = "Soil moisture (m3/m3)", xaxt = "n")
lines(th_mod_wide$date, th_mod_wide$theta_2, type = "l", col = 2)
lines(th_mod_wide$date, th_mod_wide$theta_3, type = "l", col = 3)
lines(th_mod_wide$date, th_mod_wide$theta_4, type = "l", col = 4)
legend("bottomright", legend = c("10cm", "20cm", "30cm", "40cm"),
       col = c(1:4), lty = 1, ncol = 2)
axis(1, at = pretty(th_mod_wide$date, n = 7), format(pretty(th_mod_wide$date, n=7), "%m/%y"))
dev.off()

# Plot mod vs obs and quality measure ####
png(paste0(project_path, "/pr_F50_NSE.png"), height=6, width=6, units = "in", res=500)
par(mfrow = c(1,1), las = 1, mar = c(4,4,2,2), oma = c(1,1,1,1))
plot(df$obs, df$mod, xlab = "Observed soil moisture [m3/m3]", ylab = "Modeled soil moisture [m3/m3]", xlim = c(0.3,0.8), ylim = c(0.3,0.8))
abline(0,1, col = 2)
legend("bottomright", legend = paste("NSE =", 1-round(GOF, digits = 2)))
dev.off()

# Plot evolution of parameter estimates ####

png(paste0(project_path, "/ppso_ths.png"), height=6, width=7, units = "in", res=500)
par(mfrow = c(2,1), mar = c(0,4,0,2), oma = c(4,1,2,0), las = 1)
plot(ppso$time, ppso$parameter_1, pch = ".", ylab = "ths", xaxt = "n")
points(ppso$time, ppso$parameter_2, pch = ".", col = 3)
points(ppso$time, ppso$parameter_3, pch = ".", col = 4)
points(ppso$time, ppso$parameter_4, pch = ".", col = 5)
plot(ppso$time, ppso$objective_function, pch = ".", ylab = "Objective function", col = 2)
dev.off()

png(paste0(project_path, "/ppso_alpha.png"), height=6, width=7, units = "in", res=500)
par(mfrow = c(2,1), mar = c(0,4,0,2), oma = c(4,1,2,0), las = 1)
plot(ppso$time, ppso$parameter_5, pch = ".", ylab = "Alpha", xaxt = "n")
points(ppso$time, ppso$parameter_6, pch = ".", col = 3)
points(ppso$time, ppso$parameter_7, pch = ".", col = 4)
points(ppso$time, ppso$parameter_8, pch = ".", col = 5)
plot(ppso$time, ppso$objective_function, pch = ".", ylab = "Objective function", col = 2)
dev.off()

png(paste0(project_path, "/ppso_n.png"), height=6, width=7, units = "in", res=500)
par(mfrow = c(2,1), mar = c(0,4,0,2), oma = c(4,1,2,0), las = 1)
plot(ppso$time, ppso$parameter_9, pch = ".", ylab = "n", xaxt = "n")
points(ppso$time, ppso$parameter_10, pch = ".", col = 3)
points(ppso$time, ppso$parameter_11, pch = ".", col = 4)
points(ppso$time, ppso$parameter_12, pch = ".", col = 5)
plot(ppso$time, ppso$objective_function, pch = ".", ylab = "Objective function", col = 2)
dev.off()

png(paste0(project_path, "/ppso_Ks.png"), height=6, width=7, units = "in", res=500)
par(mfrow = c(2,1), mar = c(0,4,0,2), oma = c(4,1,2,0), las = 1)
plot(ppso$time, ppso$parameter_13, pch = ".", ylab = "Ks", xaxt = "n")
points(ppso$time, ppso$parameter_14, pch = ".", col = 3)
points(ppso$time, ppso$parameter_15, pch = ".", col = 4)
points(ppso$time, ppso$parameter_16, pch = ".", col = 5)
plot(ppso$time, ppso$objective_function, pch = ".", ylab = "Objective function", col = 2)
dev.off()

