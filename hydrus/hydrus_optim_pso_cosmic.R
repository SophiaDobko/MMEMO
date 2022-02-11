######################################################################+
# Externally calibrate Hydrus-1D using Particle Swarm Optimization ####
######################################################################+

library(hydrusR)
library(xts)
library(dplyr)
library(soilwater)
library(tidyverse)
library(lhs)
library(hydroGOF)

# Directories and basic information ####
# Adapt information in this first chapter 

source("C:/users/bauers/data/MMEMO/hydrus/hydpara_selector.R")

# define wd and other directories ####
#comp <- Sys.getenv("COMPUTERNAME")
#if(comp == "GK-NB-5"){
#  setwd("D:/CosmicSense/Hydrus")
#  hydrus_path <- "C:/hydrus"
#  theta_observed <- "C:/Users/dobkowitz/data/crns/merged/rhi_pr2_FT.txt" # observed pr2 soil moisture timeseries, has columns datetime, theta_1, theta_2, theta_3, theta_4
#  }
setwd("C:/users/bauers/data/Hydrus-1D")

# Choose pr2 profile
#pr_plot <- "f50" # or "center", "t25", "f25", "t50", "f50"

# Basic inputs
project_name = "spo_cosmic_5"
#if (pr_plot == "t50"){
#  project_name = "Hydrus_rhi_pr_T50"
#}
project_path = paste0(getwd(), "/Projects/", project_name)
TimeUnit = "days"
SpaceUnit = "cm"
PrintTimes = 1
tz = "etc/GMT-1"
start_date = as.Date("2020-01-01")
end_date = as.POSIXct("2022-02-03")
nbefore = 0 # days in the year before first boundary conditions day

# number of soil horizons
nhor = 3

#para <- data.frame(thr = c(0, 0, 0, 0), # initial soil hydraulic parameters for all materials
#                   ths = c(0.797, 0.73, 0.741, 0.891),
#                   Alpha = c(0.02, 0.012, 0.005, 0.003),
#                   n = c(1.23, 1.12, 1.15, 1.16),
#                   Ks = c(14, 6.69, 0.92, 104),
#                   l = c(0.5, 0.5, 0.5, 0.5))

# initial estimates of parameters which will be calibrated in pso_
para_initial = cbind(c(para$ths, para$Alpha, para$n, para$Ks)) 


## Observed CRNS data ####
crns0 <- read.table("C:/Users/bauers/data/crns/Spo0_df_daily.txt", header=T, sep=",")
crns0$datetime <- as.POSIXct(crns0$datetime)
crns1 <- read.table("C:/Users/bauers/data/crns/Spo1_df_daily.txt", header=T, sep=",")
crns1$datetime <- as.POSIXct(crns1$datetime)
crns2 <- read.table("C:/Users/bauers/data/crns/Spo2_df_daily.txt", header=T, sep=",")
crns2$datetime <- as.POSIXct(crns2$datetime)

# rescale neutron counts of the different crns probes
crns1$cphc <- crns1$cphc * 3.23
crns0$cphc <- crns0$cphc * 3.23 / 0.9726

crns <- rbind(crns0[,c(1,20)], crns1[,c(1,17)], 
              crns2[c((which(crns2$datetime==max(crns1$datetime))+1):which(crns2$datetime=="2022-02-01")),c(1,20)])

## Objective function ####
# adapt soil_para due to number of soil materials and parameters to be calibrated
# adapt if modeled soil moisture timeseries has not/more then 4 depths


zielfun = function(params) { # params is produced by ppso, dont worry about this
  soil_para <- data.frame("thr" = para$thr,
                     ths = c(params[1:nhor]),
                     Alpha = c(params[(nhor+1):(nhor+3)]),
                     n = c(params[(2*nhor+1):(2*nhor+3)]),
                     Ks = c(params[(3*nhor+1):(3*nhor+3)]),
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
    # read simulated neutron counts
    mod_cosmic <- read.table(paste0(project_path,"/Cosmic.out"), skip=4)
    names(mod_cosmic) <- c("Time","NFlux")
    # set datetime for hydrus output
    mod_cosmic$Time <- seq.POSIXt(as.POSIXct("2020-01-01"),by="DSTday", length.out=length(mod_cosmic$Time))
    # take only modelled data for which exists observed data
    mod_cosmic_fit <- mod_cosmic[c(which(mod_cosmic$Time==min(crns0$datetime)):which(mod_cosmic$Time==max(crns0$datetime)),
                                   which(mod_cosmic$Time==min(crns1$datetime)):which(mod_cosmic$Time=="2022-02-01")),]
    
    #merge observed and modelled data into 1 dataframe
    df <- data.frame("obs" = crns$cphc, "mod" = mod_cosmic_fit$NFlux)
    
    # Calculate Goodness Of Fit
    #GOF = sum((df[,1]-df[,2])^2, na.rm=T)/sum((df[,1]-mean(df[,1], na.rm=T))^2) #nash-sutcliff, without 1-..., so that minimum ist best
    GOF = NSE(df[,1], df[,2]) # with this option optim = max_search, or GOF=-1*(NSE(...)-1)
  }
  return(GOF)
}


## Optimization ###########
####################
library(ppso)

hydrus_path = "C:/Users/bauers/Hydrus-1D 4.xx"
start_time = Sys.time()
res=optim_pso(zielfun, number_of_parameters=4*nhor, 
                                          # Qs            Alpha             n           Ks
              parameter_bounds = cbind(c(rep(0.35,nhor), rep(0.002,nhor), rep(0.9,nhor), rep(0.5,nhor)), #lower limit
                                       c(rep(0.42,nhor), rep(1.5,nhor), rep(4,nhor), rep(800,nhor))), #upper limit
              max_number_of_iterations=200, 
              logfile=paste0(project_path, "/ppso.log"), 
              projectfile=paste0(project_path, "/ppso.pro"),
              max_number_function_calls = 50,  #5000
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
  ppso <- ppso[is.finite(ppso$objective_function),]
  best <- ppso[ppso$objective_function == max(ppso$objective_function),]
  params <- as.numeric(best[1,2:(4*nhor+1)])
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

