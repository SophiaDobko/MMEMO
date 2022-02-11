#' Write profile probe observations in FIT.IN
#'
#' This function writes profile probe observations in fit.in
#' 
#' @param project.path directory of the folder containing the hydrus input files
#' @param para dataframe containing the six columns ```thr```, ```ths```, ```Alfa```, ```n```, ```Ks``` and ```l```, containing the parameters for material 1,2...n
#' @return ```selector.in``` with the hydraulic soil parameters for the next calibration run
#' @importFrom data.table fwrite
#' @export


library(hydrusR)
#install.packages("tidyverse")
library(tidyverse)


comp <- Sys.getenv("COMPUTERNAME")
if(comp == "GK-NB-5"){
  setwd("D:/CosmicSense/Hydrus")
  crns_dir <- "C:/Users/Dobkowitz/data/crns/"}

tz = "etc/GMT-1"
start_date <- as.POSIXct("2021-04-18", tz)
end_date <- as.POSIXct("2021-11-24", tz)
nbefore = 107 # number of days in the year before start_date
project_name = "Hydrus_inv_pr1"
project_path = paste0(getwd(), "/", project_name)
products_dir <- "C:/Users/dobkowitz/git/cosmicsense-notebooks/notebooks/sponheim_rhinluch/timeseries_product/"
pr_plot <- "center" # or "t25", "f25", "t50", "f50"

# Read in and formatation of pr data
pr_prepare <- function(products_dir, project_path, pr_plot, start_date, end_date, tz = "etc/GMT-1"){
  pr <- read.table(paste0(products_dir,"Rhi1_pr_product.txt"), sep = ",", header = T)
  pr <- pr[pr$plot==pr_plot,]
  pr$datetime <- as.POSIXct(pr$datetime, format = "%Y-%m-%d %H:%M:%S", tz)
  pr$date <- as.Date(pr$datetime)
  pr <- pr[pr$datetime>= start_date & pr$datetime<= end_date,]
  pr <- aggregate(x = pr[,4:7], by = list(pr$date), FUN = mean)
  colnames(pr)[1] <- "date"
  pr$date <- 1:nrow(pr)
  
  prtidy <- gather(pr, key = "probe", value = "FOS","theta_1", "theta_2", "theta_3", "theta_4")
  #prtidy$POS <- paste0(substr(prtidy$probe,7,7),"0")
  prtidy$POS[prtidy$probe=="theta_1"] <- 2
  prtidy$POS[prtidy$probe=="theta_2"] <- 3
  prtidy$POS[prtidy$probe=="theta_3"] <- 4
  prtidy$POS[prtidy$probe=="theta_4"] <- 5
  prtidy$'ITYPE(N)' <- rep(2,nrow(prtidy))
  prtidy$WTS <- rep(1, nrow(prtidy))
  colnames(prtidy)[1] <- "HO(N)"
  prtidy <- prtidy[,c(1,3,5,4,6)]
  prtidy$`HO(N)` <- prtidy$`HO(N)` + nbefore # add number of days in the year before start_date
  return(prtidy)
}


# write pr soil moisture observations to fitin
pr_fitin <- function(para, products_dir, project_path, pr_plot, start_date, end_date, tz = "etc/GMT-1"){
 
  fitin = readLines(con = paste0(project_path, "/FIT.IN"),
                    n = -1L,
                    encoding = "unknown")
  cat(fitin[1:(10+5*nrow(para))],
      file = paste0(project_path, "/FIT.IN"),
      fill = F,
      sep = "\n")
  
  fwrite(x = prtidy,
         file = paste0(project_path, "/FIT.IN"),
         sep = "\t",
         col.names=T,
         append=T)
  
  cat(fitin[length(fitin)],
      file = paste0(project_path, "/FIT.IN"),
      append = T)
  
}






