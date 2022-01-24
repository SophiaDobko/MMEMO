#####################################
#### Calibrating parameter N in cosmic.in ####
#### (Cosmic operator HYDRUS 1D) ####
#####################################

# hydrusR from: https://github.com/shoebodh/hydrusR
# -> Use R to prepare hydrus input and analyse hydrus output

# install.packages("remotes")
# remotes::install_github("shoebodh/hydrusR")

library(hydrusR)
library(hydroGOF)

setwd("C:/Users/bauers/data/Hydrus-1D/Hydrus-R/Prepare_Hydrus_Input")

project.path = "C:/Users/bauers/data/Hydrus-1D/Projects/spo_cosmic_2"

# Set bulk densitiy and lattice water in cosmic.in file
bd <- 1.46 # buld density
lw <- 0.023 # lattice water
alpha <- 0.404 - 0.101*bd
L3 <- -31.65 + 99.29*bd
N <- 611  # initial value 

cosmic.in <- read.table(paste0(project.path,"/Cosmic.in"), sep="\t")
cosmic.in[,1] <- c(bd, lw, N, alpha, cosmic.in[5:6,1], L3, cosmic.in[8,1])
write.table(cosmic.in, paste0(project.path,"/Cosmic.in"), sep="\t", row.names=F,
            col.names=F, quote=F)


# Run hydrus model ####
call.H1D(project.path,
         hydrus.path = "C:/Users/bauers/Hydrus-1D 4.xx",
         show.output = TRUE)

##### Default hydrus path in Windows

# run.H1D.simulation(project.path = project_path,
# profile.depth = profile_depth,
# beginT = beginTime, endT = endTime, deltaT = time_step,
# bot.bc.type = bot_bc_type, bot.bc.value = const_botFlux,
# const.bot.bc = TRUE,atm.bc.data = atm_bc_data, TimeUnit = TimeUnit,
# show.output = T)


##### Read CRNS data #####
crns0 <- read.table("C:/Users/bauers/data/crns/Spo0_df24_hourly.txt", header=T, sep=",")
crns0$datetime <- as.POSIXct(crns0$datetime)
crns1 <- read.table("C:/Users/bauers/data/crns/Spo1_df24_hourly.txt", header=T, sep=",")
crns1$datetime <- as.POSIXct(crns1$datetime)
crns2 <- read.table("C:/Users/bauers/data/crns/Spo2_df24_hourly.txt", header=T, sep=",")
crns2$datetime <- as.POSIXct(crns2$datetime)

# rescale neutron counts of the different crns probes
crns1$cphc <- crns1$cphc * 3.23
crns0$cphc <- crns0$cphc * 3.23 / 0.9726

crns <- rbind(crns0[,c(1,20,21)], crns1[,c(1,17,18)], 
              crns2[c((which(crns2$datetime==max(crns1$datetime))+1):which(crns2$datetime=="2021-08-28")),c(1,20,21)])

##### Display model output #####
mod_cosmic <- read.table(paste0(project.path,"/Cosmic.out"), skip=4)
mod_flux <- read.table(paste0(project.path,"/T_Level.out"), skip=9, 
                       nrow=length(readLines(paste0(project.path,"/T_Level.out"))) - 10)
mod_obsnode <- read.table(paste0(project.path,"/Obs_Node.out"), skip=11, 
                          nrow=length(readLines(paste0(project.path,"/Obs_Node.out"))) - 12)

names(mod_cosmic) <- c("Time","NFlux")
names(mod_flux) <- c("Time","rTop","rRoot","vTop","vRoot","vBot","sum_rRoot","sum_rRoot","sum_vTop","sum_vRoot","sum_Bot","hTop","hRoot","hBot",
                     "RunOff","sum_RunOff","Volume","sum_Infil","sum_Evap","TLevel","Cum_WTrans","SnowLayer","01","02","03","04")
names(mod_obsnode) <- c("Time",paste0(c("h_","theta_","Temp_"),c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6)))

# set datetime for hydrus output
mod_cosmic$Time <- seq.POSIXt(as.POSIXct("2020-01-01"),by="hours", length.out=length(mod_cosmic$Time))
mod_flux$Time <- seq.POSIXt(as.POSIXct("2020-01-01"),by="hours", length.out=length(mod_cosmic$Time))
mod_obsnode$Time <- seq.POSIXt(as.POSIXct("2020-01-01"),by="hours", length.out=length(mod_cosmic$Time))


# Gütemaße für Modell
# Modellergebnisse für die es auch Daten gibt:
mod_cosmic_fit <- mod_cosmic[c(which(mod_cosmic$Time==min(crns0$datetime)):which(mod_cosmic$Time==max(crns0$datetime)),
                               which(mod_cosmic$Time==min(crns1$datetime)):which(mod_cosmic$Time=="2021-08-28")),]

(rmse(crns$cphc, mod_cosmic_fit$NFlux))

lm.daten.modell <- lm(mod_cosmic_fit$NFlux ~ crns$cphc)
(r_2 <- summary(lm.daten.modell)$r.squared)

(nse <- NSE(mod_cosmic_fit$NFlux, crns$cphc))

# mean of data
mean_mod <- mean(mod_cosmic_fit$NFlux)
mean_crns <- mean(crns$cphc, na.rm=T)


#### Loop for calibrating N in cosmic.in ####
nse_old <- nse-0.2
docu <- data.frame(N=cosmic.in[3,1],
                   NSE=nse,
                   R_2=r_2)
while #(abs(nse - nse_old)>0.005) {
(nse > nse_old) {
  nse_old <- nse
  
  # Read cosmic.in file and set parameter N
  cosmic.in <- read.table(paste0(project.path,"/Cosmic.in"), sep="\t")
  ifelse(mean_mod > mean_crns,   # condition
         cosmic.in[3,1] <- cosmic.in[3,1]-1,  # yes
         cosmic.in[3,1] <- cosmic.in[3,1]+1)  # no
  write.table(cosmic.in, paste0(project.path,"/Cosmic.in"), sep="\t", row.names=F,
              col.names=F, quote=F)
  
  # Run hydrus model ####
  call.H1D(project.path,
           hydrus.path = "C:/Users/bauers/Hydrus-1D 4.xx",
           show.output = TRUE)
  
  # Read model output
  mod_cosmic <- read.table(paste0(project.path,"/Cosmic.out"), skip=4)
  names(mod_cosmic) <- c("Time","NFlux")
  mod_cosmic$Time <- seq.POSIXt(as.POSIXct("2020-01-01"),by="hours", length.out=length(mod_cosmic$Time))
  
  mod_cosmic_fit <- mod_cosmic[c(which(mod_cosmic$Time==min(crns0$datetime)):which(mod_cosmic$Time==max(crns0$datetime)),
                                 which(mod_cosmic$Time==min(crns1$datetime)):which(mod_cosmic$Time=="2021-08-28")),]
  
  # mean of data
  mean_mod <- mean(mod_cosmic_fit$NFlux)
  mean_crns <- mean(crns$cphc, na.rm=T)
  
  # Gütemaße
  lm.daten.modell <- lm(mod_cosmic_fit$NFlux ~ crns$cphc)
  (r_2 <- summary(lm.daten.modell)$r.squared)
  
  (nse <- NSE(mod_cosmic_fit$NFlux, crns$cphc))
  
  # write docu
  docu <- rbind(docu, c(cosmic.in[3,1], nse, r_2))
  
}
# write the last N in cosmic.in
cosmic.in[3,1] <- docu[nrow(docu)-1, 1]
write.table(cosmic.in, paste0(project.path,"/Cosmic.in"), sep="\t", row.names=F,
            col.names=F, quote=F)
#write.table(docu, "docu_N_fit.txt", sep="\t", row.names=F,col.names=T, quote=F)



##### plot model output against data ####
#oldpar <- par()
par(oldpar)
plot(mod_cosmic, type="l", ylim=c(2200,3100))
lines(crns0$datetime, crns0$cphc, col=3)
lines(crns1$datetime, crns1$cphc, col=4)
lines(crns2$datetime, crns2$cphc, col=2)

# Plot model fluxes 
par(las=1, mgp=c(2,0.7,0), mar=c(3,3,1,1))
plot(mod_flux$Time, mod_flux$sum_rRoot, type="l", ylim=c(-30,30))
lines(mod_flux$Time, mod_flux$sum_vRoot, col=2)
lines(mod_flux$Time, mod_flux$sum_vTop, col=4)
lines(mod_flux$Time, mod_flux$sum_Bot, col=6)
legend("topleft",c())

# Plot theta of observation nodes
plot(mod_obsnode$Time, mod_obsnode[,3], type="l", ylim=c(0.1,0.4))
grid()
for (i in 2:6) {
  lines(mod_obsnode$Time, mod_obsnode[,3*i], col=i)
}

