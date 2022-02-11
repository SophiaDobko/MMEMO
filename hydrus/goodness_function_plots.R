# Goodness function and plots for Hydrus model
#
#end time for crns data has to be updated!!

library(hydroGOF)

##### Read CRNS data #####
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

## Function for computing goodness
goodness <- function(project.path) {
  
  # Display model output #####
  mod_cosmic <- read.table(paste0(project.path,"/Cosmic.out"), skip=4)
  mod_flux <- read.table(paste0(project.path,"/T_Level.out"), skip=9, 
                         nrow=length(readLines(paste0(project.path,"/T_Level.out"))) - 10)
  mod_obsnode <- read.table(paste0(project.path,"/Obs_Node.out"), skip=11, 
                            nrow=length(readLines(paste0(project.path,"/Obs_Node.out"))) - 12)
  
  names(mod_cosmic) <- c("Time","NFlux")
  names(mod_flux) <- c("Time","rTop","rRoot","vTop","vRoot","vBot","sum_rTop","sum_rRoot","sum_vTop","sum_vRoot","sum_Bot","hTop","hRoot","hBot",
                       "RunOff","sum_RunOff","Volume","sum_Infil","sum_Evap","TLevel","Cum_WTrans","SnowLayer","01","02","03","04")
  names(mod_obsnode) <- c("Time",paste0(c("h_","theta_","Temp_"),c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6)))
  
  # set datetime for hydrus output
  mod_cosmic$Time <- seq.POSIXt(as.POSIXct("2020-01-01"),by="DSTday", length.out=length(mod_cosmic$Time))
  mod_flux$Time <- seq.POSIXt(as.POSIXct("2020-01-01"),by="DSTday", length.out=length(mod_cosmic$Time))
  mod_obsnode$Time <- seq.POSIXt(as.POSIXct("2020-01-01"),by="DSTday", length.out=length(mod_cosmic$Time))
  
  
  # Gütemaße für Modell
  # Modellergebnisse für die es auch Daten gibt:
  mod_cosmic_fit <- mod_cosmic[c(which(mod_cosmic$Time==min(crns0$datetime)):which(mod_cosmic$Time==max(crns0$datetime)),
                                 which(mod_cosmic$Time==min(crns1$datetime)):which(mod_cosmic$Time=="2022-02-01")),]
  
  rmse <- rmse(crns$cphc, mod_cosmic_fit$NFlux)
  
  lm.daten.modell <- lm(mod_cosmic_fit$NFlux ~ crns$cphc)
  r_2 <- summary(lm.daten.modell)$r.squared
  
  nse <- NSE(mod_cosmic_fit$NFlux, crns$cphc)
  
  # mean of data
  mean_mod <- mean(mod_cosmic_fit$NFlux)
  mean_crns <- mean(crns$cphc, na.rm=T)
  
  
  ##### plot model output against data ####
  #oldpar <- par()
  par(las=0, mfrow=c(3,1), mgp=c(2.5,0.7,0), mar=c(1.8,3.5,0.5,1))
  plot(mod_cosmic, type="l", ylim=c(2100,3200))
  lines(crns0$datetime, crns0$cphc, col=3)
  lines(crns1$datetime, crns1$cphc, col=4)
  lines(crns2$datetime, crns2$cphc, col=2)
  legend("topright", c("model data", "CRNS 0", "CRNS 1", "CRNS 2"), col=c(1,3,4,2),
         lwd=2, cex=0.9)
  
  # Plot model fluxes 
  #par(las=1, mgp=c(2.5,0.7,0), mar=c(2,3.5,1,1))
  plot(mod_flux$Time, mod_flux$sum_vRoot, type="l", ylim=c(-35,30),
       xlab = "Time", ylab = "Flux [cm]")
  lines(mod_flux$Time, mod_flux$sum_vTop, col=2)
  lines(mod_flux$Time, mod_flux$sum_Bot, col=4)
  legend("topleft",c("Root Water Uptake", "Top Flux", "Bottom Flux"),
         col=c(1,2,4), lwd=2, cex=0.9)
  
  # Plot theta of observation nodes
  plot(mod_obsnode$Time, mod_obsnode[,3], type="l", ylim=c(0.09,0.4),
       xlab = "Time", ylab = "Theta")
  grid(nx=0, ny=NULL)
  for (i in 2:6) {
    lines(mod_obsnode$Time, mod_obsnode[,3*i], col=i)
  }
  abline(v=range(crns0$datetime), lty="dashed")
  abline(v=c(min(crns1$datetime), max(crns$datetime)),  lty="dashed")
  legend("bottomright", c("10cm","20cm","30cm","40cm","60cm","100cm","calibrating period"),
         col=c(1:6), cex=0.8, bty="o", ncol=3, lty=c(rep("solid",6),"dashed"))
  
  goods <- c(rmse, r_2, nse)
  
  return(goods)

}
