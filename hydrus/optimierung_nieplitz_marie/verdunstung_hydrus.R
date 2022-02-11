### Marie-Therese Schmehl

# Nieplitz evaporation calculation, from dwd data 
#Langenlipsdorf und Wiesenburg

###

setwd("P:/nieplitz/Visualisierung")
#setwd("/home/marie/Projekte/nieplitz/Visualisierung")

if (!require(lubridate))
{  
  install.packages("lubridate", dependencies = TRUE)
  library(lubridate)
}



library(xts)

#Evapotranspiration package ######

if (!require(Evapotranspiration))
{  
  install.packages("Evapotranspiration", dependencies = TRUE)
  library(Evapotranspiration)
}

library(xts)
#library(ggplot2)

#daten einlesen######
tempmin = read.csv("../externe_Daten/dwd/temp/TMIN/data/data_TNK_MN004.csv", header=T)
tempmin$datum = as.POSIXct(strptime(tempmin$Zeitstempel, format="%Y%m%d"), tz="UTC")
#unique(tempmin$SDO_ID) #2856 = Langenlipsdorf, 5546=Wiesenburg
tempmin = subset(tempmin[,-c(3,5,6)], SDO_ID == 5546)

tempmax = read.csv("../externe_Daten/dwd/temp/TMAX/data/data_TXK_MN004.csv", header=T)
tempmax$datum = as.POSIXct(strptime(tempmax$Zeitstempel, format="%Y%m%d"), tz="UTC")
#unique(tempmax$SDO_ID) #2856 = Langenlipsdorf, 5546=Wiesenburg
tempmax = subset(tempmax[,-c(3,5,6)], SDO_ID == 5546)

sonne = read.csv("../externe_Daten/dwd/sonnenscheindauer/cdc_download_2020-02-26_14_00/data/data_SDK_MN004.csv", header=T)
sonne$datum = as.POSIXct(strptime(sonne$Zeitstempel, format="%Y%m%d"), tz="UTC")
#unique(sonne$SDO_ID) #2856 = Langenlipsdorf, 5546=Wiesenburg
sonne = subset(sonne[,-c(3,5,6)], SDO_ID == 5546)

relLF = read.csv("../externe_Daten/dwd/humidity/cdc_download_2020-02-27_12_40/data/data_RF_TU_MN009.csv", header=T)
relLF$datum = as.POSIXct(strptime(relLF$Zeitstempel, format="%Y%m%d%H%M"), tz="UTC")
#unique(relLF$SDO_ID) #2856 = Langenlipsdorf, 5546=Wiesenburg
relLF = subset(relLF[,-c(3,5,6)], SDO_ID == 5546)
RHmax = apply.daily(xts(relLF$Wert, order.by = relLF$datum), max, na.rm=T)
RHmin = apply.daily(xts(relLF$Wert, order.by = relLF$datum), min, na.rm=T)

cloud = read.csv("../externe_Daten/dwd/cloud/cdc_download_2020-02-27_14_52/data/data_NM_MN004.csv", header=T)
cloud$datum = as.POSIXct(strptime(cloud$Zeitstempel, format="%Y%m%d"), tz="UTC")
#unique(cloud$SDO_ID) #2856 = Langenlipsdorf, 5546=Wiesenburg
cloud = subset(cloud[,-c(3,5,6)], SDO_ID == 5546)

wind = read.csv("../externe_Daten/dwd/wind/cdc_download_2020-02-27_14_48/data/data_FM_MN003.csv", header=T)
wind$datum = as.POSIXct(strptime(wind$Zeitstempel, format="%Y%m%d"), tz="UTC")
#unique(wind$SDO_ID) #2856 = Langenlipsdorf, 5546=Wiesenburg
wind = subset(wind[,-c(3,5,6)], SDO_ID == 5546)


#Daten prozessieren######
data("constants") #load constants provided by R Package
#change constants applying for Wiesenburg
constants$lat = 52.1207 #latitude Wiesenburg
constants$lat_rad = (constants$lat * pi) / (180)
constants$as = 0.25 # Empfehlung FAO bei nciht vorhanden Werten
constants$Elev = 187.0 #Höhe über Meeresspiegel
#hydiff = 60*60*24*61 

#Datum
date <- setNames(data.frame(do.call("rbind",strsplit(gsub("\\(|\\)|,","",tempmax$datum),split="-"))),c("Year","Month","Day"))
date = apply(date,2,as.numeric)

#Klimadaten  
clim_xts = merge.xts(
  Tmax = xts(tempmax$Wert, order.by = tempmax$datum), 
  Tmin = xts(tempmin$Wert, order.by = tempmin$datum), 
  RHmax , 
  RHmin ,  
  n = xts(sonne$Wert, order.by = sonne$datum),
  #Cd = xts(cloud$Wert, order.by = cloud$datum),
  uz = xts(wind$Wert, order.by = wind$datum)
)
clim_xts = apply.daily(clim_xts, mean, na.rm=T)  #all values into 1 line
clim = ReadInputs(colnames(clim_xts), cbind(date, as.data.frame(clim_xts)), 
                  constants , stopmissing = c(20, 20, 20))
  
#Berechnung von ET

#clim$RHmin=NULL
#clim$RHmax=NULL

ETmethods = c("ET.Turc", "ET.Makkink", "ET.Abtew", "ET.Hamon", "ET.HargreavesSamani",
              "ET.Romanenko")
ETmethods = t(read.table("../externe_Daten/ET.methods.txt", sep=","))
ETmethods = ETmethods[-c(1:3,5,9,11:14,18,19:21)]

results = list()
ET_year = list()

for (i in 1:length(ETmethods)) {
  results[[i]] = do.call(ETmethods[i], list(clim, constants, ts="daily", 
                                            solar="sunshine hours",
                                            AdditionalStats=F, save.csv = "yes"))
  plot(results[[i]]$ET.Daily, main = ETmethods[i])
  abline(h=0, col=2)
  ET_year[[i]] = apply.yearly(xts(results[[i]]$ET.Daily, order.by = index(results[[i]]$ET.Daily) ), sum)

}

ET_year = as.data.frame(ET_year, row.names = index(ET_year[[1]]))
colnames(ET_year) = ETmethods
barplot(as.numeric(ET_year[2,]))
barplot(as.numeric(ET_year[3,]), width=0.5, space=1.4, add=T, col=2)

results[[14]] = ET.Turc(clim, constants, ts="daily", 
                         solar="sunshine hours", humid = T,
                         AdditionalStats=F, save.csv = F)
plot(results[[14]]$ET.Daily, main = ETmethods[i])

