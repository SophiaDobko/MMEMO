#Optimierung einer einfachen Simulation mit HYDRUS-1D
# Marie-Therese Schmehl

#notwendige Pakete laden/installieren
  devtools::install_github("shoebodh/hydrusR")
  library(hydrusR)

if (!require(dplyr))
{  
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
} 
if (!require(xts))
{  
  install.packages("xts", dependencies = TRUE)
  library(xts)
}  


#set working directory - BITTE ANPASSEN
setwd("/home/marie/ownCloud/Numerik/Optimierung_nieplitz")   #setwd("P:/nieplitz/hydrus")
source("functions/write_atmosph_template_path.R")
source("functions/write_ini_cond.R")
#set timezone for equal setting in every data.frame/data
tz = "etc/GMT-1" 

## BEOBACHTUNGSDATEN #########
th_obs= read.table("theta_mm.txt", header=T)
th_obs$Index = trunc.POSIXt(as.POSIXct(th_obs$Index, format="%Y-%m-%d %H:%M:%S", tz = tz),
                            units = "hours") #setze alle Zeitstempel auf HH:00:00; für Vergleich mit th_mod notwendig
th_obs = th_obs[,1:2] #lösche unbenutzte Spalten
th_obs$Sonde.50.04 = th_obs$Sonde.50.04/10 #in cm umrechnen
th_obs = th_obs[which(!is.na(th_obs$Index)),]
th_obs = th_obs[which(!is.na(th_obs$Sonde.50.04)),]


#### HYDRUS - unveränderter Teil##########
###############################

## HYDRUS-Projekt mit "Basic" und "Profile" information manuell in HYDRUS-GUI anlegen
# "Initial conditions" und "variable boundary conditions" werden hier im Skript angepasst

# hier Pfad zum angelegten Hydrusprojekt angeben
project_path =  "P:/nieplitz/hydrus/niep18_2"

#Start- und Enddatum der zu simulierende Zeitspanne festlegen
#Länge muss mit Zeitspanne in HydrusProjekt übereinstimmen
start_date = as.POSIXct("2018-01-01", tz = tz)
end_date = as.POSIXct("2018-12-31", tz = tz)


### Atmospheric top boundary conditions ####
# hier Niederschlag und potential Evapotranspiration (pET)

## potential Evapotranspiration (pET) input vorbereiten
et_rate = read.csv("verdunstung/ET_HargreavesSamani.csv", nrows = 911) 
et_rate$date = as.POSIXct(row.names(et_rate), tz = tz) 
et_rate = et_rate[which(et_rate$date >= start_date  & et_rate$date <= end_date),] #auf zu simulierende Zeitspanne

## Precipitation input vorbereiten 
prec = read.table("n_day_frohnsdorf_2019-11-07", header = T, sep=";")
prec$date = as.POSIXct(prec$date, tz = tz)
prec = prec[which(prec$date >= start_date  & prec$date <= end_date),]
prec$x[is.na(prec$x)] = 0 #NAs in Niederschlag, geben Fehler in Hydrus -> auf numerischen Wert setzen

## numerischen Zeitvektor aus input data für Hydrus anlegen
time_step_et = unique(diff(et_rate$date)) 
  if(length(time_step_et)!=1) print("timestep is not unique. Check timesteps of pET data")
time_step_prec = unique(diff(prec$date)) 
  if(length(time_step_prec)!=1) print("timestep is not unique. Check timesteps of precipitation data")
  if(time_step_et != time_step_prec) print("timesteps of pET and prec must be equal, please check")
time_step = time_step_et; rm(time_step_et); rm(time_step_prec)
endTime = time_step+as.numeric(difftime(end_date, start_date))
tAtm = seq(time_step, endTime, time_step)
ntimesteps = length(tAtm)

## atmospheric boundary condition data festlegen
atm_bc_data = data.frame(tAtm = tAtm,
                         Prec = prec$x/10, #in cm umrechnen
                         rSoil = et_rate$ET.Daily/10, #in cm umrechnen
                         rRoot = numeric(ntimesteps),
                         hCritA = rep(10000, ntimesteps),
                         rB = numeric(ntimesteps),
                         hB = numeric(ntimesteps),
                         ht = numeric(ntimesteps),
                         RootDepth = numeric(ntimesteps))

write.atmosph.in(project.path = project_path, maxAL = endTime, deltaT = time_step,
                 atm.bc.data = atm_bc_data)

## Initial condition ####
#Anzahl der Profilknoten. Muss mit Angabe im HydrusProjekt übereinstimmen
profile_nodes = 11

#ein WErt für ganzes Profil aus Beobachtungszeitreihe ablesen
th.ini = th_obs[which(th_obs$Index >= start_date),] 
th.ini = mean(th.ini$Sonde.50.04[1:24])/100 #th_obs ist in cm angegeben, daher /100 -> %
pr.vec = rep(th.ini,length(profile_nodes)) 

write.ini.cond(project.path = project_path, pr.vec = pr.vec)



## ZIELFUNKTION ##########
######################
# bei Änderung der Anzahl der Paramter, diese im Skript in soil_para anpassen
zielfun = function(params) {
  #wandle vector in list() wie für run.H1D.simulation verwendbar
  soil_para = list(thr = 0.036,
                ths = 0.43, 
                Alfa = params[1] , #
                n = params[2] , #
                Ks = params[3],
                l = 0.45) #params[4]
  
  #passe soil_para für neuen Durchlauf an
  write.hydraulic.para(project.path = project_path, para = soil_para) #schreibt in selector.in
  
  #führe Hydrus aus
  call.H1D(project_path)
  
  #wenn Fehler in Hydrus Berechnung, setze GOF=Inf
  if (file.exists(file.path(project_path, "Error.msg"))) { 
    file.remove(file.path(project_path, "Error.msg"))
    GOF = Inf
  }
  else {
    #lese modellierte Werte aus
    th_mod = read.tlevel.out(project_path)
    th_mod$DateTime = trunc.POSIXt(as.POSIXct(th_mod$Time*24*60*60, tz = tz, origin=start_date),
                                   units = "hours") #sets all timestamps to same hour
    th_mod = th_mod[,c(23,17)]
    
    #merge observed and modelled data into 1 dataframe
    df = merge.xts(xts(th_obs$Sonde.50.04, order.by = th_obs$Index),
                   period.apply(xts(th_mod$Volume, order.by = th_mod$DateTime),
                                INDEX = endpoints(th_mod$DateTime, "hours"), mean),
                   join = "right")
    colnames(df) = c("obs", "mod")

    #Berechne Gütemaß (Goodness-of-fit)
    GOF = sum((df[,1]-df[,2])^2)/sum((df[,1]-mean(df[,1]))^2) #nash-sutcliff, ohne 1-..., damit minimum gegeben
    #GOF = NSE(df[,1], df[,2]) # hierbei optim auf max-suche stellen, oder GOF=-1*(NSE(...)-1)
  }
  return(GOF)
}


## OPTIMIERUNG ###########
####################
library(ppso)

#Funktion zum plotten des Ergebnisses mit den besten Paramtern
#Anzahl parameter ggf. im Skript anpassen
source("functions/plot_best_par.R") 

#Standard paramtersätze für versch. Bodentypen in HYDRUS implementiert
#soil_para_clay = 0,068	0,38	0,008	1,09	4,8	0,5	
#soil_para_loam = 0,078	0,43	0,036	1,56	24,96	0,5
#soil_para_sand = 0,045	0,43	0,145	2,68	712,8	0,5																																																																																														

start_time = Sys.time()
res=optim_pso(zielfun,number_of_parameters=3, 
              parameter_bounds = cbind(c(0.08, 2.5, 300),c(.15, 2.75, 600)), 
              max_number_of_iterations=200, 
              logfile=paste0(project_path, "ppso.log"), 
              projectfile=paste0(project_path, "ppso.pro"),
              max_number_function_calls = 5000)
print(paste("Parameter:",paste(format(res$par),collapse=" / "),"Güte:", format(res$value),"Dauer:", format(Sys.time()-start_time)))

plot.best(res$par)


# Nur ERGEBNISSE einlesen ############
#Pfade der Ergebisordner einlesen
project_path =  grep(dir(path=getwd(), pattern="niep18", full.names = TRUE),
                     pattern=".h1d|.RData", invert=TRUE, value=TRUE)
start_date = as.POSIXct("2018-01-01", tz = tz)

#für alle Ergebnisse ausführen
for(i in 1:length(project_path)) {
  #lese modellierte Werte aus
  th_mod = read.tlevel.out(project_path[i])
  th_mod$DateTime = trunc.POSIXt(as.POSIXct(th_mod$Time*24*60*60, tz = tz, origin=start_date),
                                 units = "hours") #sets all timestamps to same hour
  th_mod = th_mod[,c(23,17,4)]
  
  #merge observed and modelled data into 1 dataframe
  df = merge.xts(xts(th_obs$Sonde.50.04, order.by = th_obs$Index),
                 period.apply(xts(th_mod[,-1], order.by = th_mod$DateTime),
                              INDEX = endpoints(th_mod$DateTime, "hours"), mean),
                 join = "right")
  colnames(df) = c("theta_obs", "theta_mod", "actual top Flux")
  #verwendete pET
  pET = sub(project_path[i], pattern=paste0(getwd(), "/niep18_"), replacement ="")
  #plot(df)
  dy_graph = list(
    d_th = dygraph(df[,1:2], main = paste("potentielle ET =", pET), ylab = "Bodenfeuchte [cm]", group = "res") %>%
      dyRangeSelector() %>% dyLegend(show = "always"),
    d_vTop = dygraph(df[,3], ylab = "Flux [cm/Tag]", group = "res") %>% dyLegend(show = "always")
  )
  show(htmltools::browsable(htmltools::tagList(dy_graph)))
  #summe real ET in mm
  rET = round(sum(df[which(df[,3]>=0),3])*10)
  print(paste("Summe der realen ET ist", rET, "mm im Modellierungszeitraum. ET_pot =", pET))
}
