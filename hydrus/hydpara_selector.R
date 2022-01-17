#' Write new hydraulic parameters in SELECTOR.IN
#'
#' This function writes during hydrus calibration runs new hydraulic soil parameters in selector.in
#' @param project.path directory of the folder containing the hydrus input file selector.in
#' @param para dataframe containing the six columns ```thr```, ```ths```, ```Alfa```, ```n```, ```Ks``` and ```l```, containing the parameters for material 1,2...n
#' @return ```selector.in``` with the hydraulic soil parameters for the next calibration run
#' @importFrom data.table fwrite
#' @export

# hydraulic parameter
#para <- data.frame(thr = c(1,2),
#                   ths = c(0,2),
#                   Alfa = c(1,2),
#                   n = c(1,2),
#                   Ks = c(1,2),
#                   l = c(1,2))

hydpara <- function(project.path, para){
  
  # write soil hydraulic parameters in selector.in
  selector = readLines(con = paste0(project.path, "/SELECTOR.IN"),
                       n = -1L,
                       encoding = "unknown")
  cat(selector[1:25],
      file = paste0(project.path, "/SELECTOR.IN"),
      fill = F,
      sep = "\n")
  
  fwrite(x = as.list(para),
         file = paste0(project.path, "/SELECTOR.IN"),
         sep = "\t",
         col.names=T,
         append=T)
  
  cat(selector[(27+nrow(para)):length(selector)],
      file = paste0(project.path, "/SELECTOR.IN"),
      fill = F,
      sep = "\n",
      append = T)
}