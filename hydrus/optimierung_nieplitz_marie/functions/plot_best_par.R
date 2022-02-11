#Funktion zum plotten von Hydrusergebnis mit bestem Parameterset
# Marie-Therese Schmehl
if (!require(dygraphs))
{  
  install.packages("dygraphs", dependencies = TRUE)
  library(dygraphs)
}  

plot.best = function(params) {
  #wandle vector in list() wie f?r run.H1D.simulation verwendbar
  soil_para = list(thr = 0.036,
                   ths = 0.43, 
                   Alfa = params[1] , #
                   n = params[2] , #
                   Ks = params[3],
                   l = 0.45) #params[4]
  
  #passe soil_para an
  write.hydraulic.para(project.path = project_path, para = soil_para) #schreibt in selector.in
  
  #rufe Hydrus auf
  call.H1D(project_path)
  
  #lese modellierte Werte aus
  th_mod = read.tlevel.out(project_path)
  th_mod$DateTime = trunc.POSIXt(as.POSIXct(th_mod$Time*24*60*60, tz = tz, origin=start_date),
                                 units = "hours") #sets all timestamps to same hour
  th_mod = th_mod[,c(23,17,4)]
  
  #merge observed and modelled data into 1 dataframe
  df = merge.xts(xts(th_obs$Sonde.50.04, order.by = th_obs$Index),
                 period.apply(xts(th_mod[,-1], order.by = th_mod$DateTime),
                              INDEX = endpoints(th_mod$DateTime, "hours"), mean),
                 join = "right")
  colnames(df) = c("theta_obs", "theta_mod", "actual top Flux")
  #plot(df)
  dy_graph = list(
    d_th = dygraph(df[,1:2]) %>%
      dyRangeSelector(),
    d_vTop = dygraph(df[,3])
  )
  d = htmltools::browsable(htmltools::tagList(dy_graph))
  #summe real ET in mm
  rET = round(sum(df[which(df[,3]>=0),3])*10)
  print(paste("Summe der realen ET ist", rET, "mm im Modellierungszeitraum"))
  
  return(d)
}