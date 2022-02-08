#' Write new hydraulic parameters in FIT.IN
#'
#' This function writes during hydrus calibration runs new hydraulic soil parameters in fit.in
#' @param project.path directory of the folder containing the hydrus input files
#' @param para dataframe containing the six columns ```thr```, ```ths```, ```Alfa```, ```n```, ```Ks``` and ```l```, containing the parameters for material 1,2...n
#' @return ```selector.in``` with the hydraulic soil parameters for the next calibration run
#' @importFrom data.table fwrite
#' @export




writefit <- function(project.path, para, mat, nmat, m, n, poptim){  # works only for model with 4 materials!
  
  # write optimized parameter of the last iteration into para
  if(!is.null(fitout$WCR))  para[m,]$thr <- fitout$WCR[nrow(fitout)]
  if(!is.null(fitout$WCS))  para[m,]$ths <- fitout$WCS[nrow(fitout)]
  if(!is.null(fitout$ALPHA)) para[m,]$Alpha <- fitout$ALPHA[nrow(fitout)]
  if(!is.null(fitout$N))   para[m,]$n <- fitout$N[nrow(fitout)]
  if(!is.null(fitout$CONDS)) para[m,]$Ks <- fitout$CONDS[nrow(fitout)]
  if(!is.null(fitout$L))   para[m,]$l <- fitout$L[nrow(fitout)]
  

  # materials
  mat[(1+4*seq(0, length.out=nmat)),] <- para
               
  mat[(2+4*seq(0, length.out=nmat)),] <- rep(0, 6)
  mat[(4*n-2),] <- poptim
  
  # write soil hydraulic parameters in fit.in
  fitin = readLines(con = paste0(project.path, "/FIT.IN"),
                       n = -1L,
                       encoding = "unknown")
  cat(fitin[1:10],
      file = paste0(project.path, "/FIT.IN"),
      fill = F,
      sep = "\n")
  
  for (i in seq(0, length.out=nmat)) {
    fwrite(x = as.list(mat[seq((i*4+1), length.out=4),]),
         file = paste0(project.path, "/FIT.IN"),
         sep = "\t",
         col.names=T,
         append=T) }
  
    cat(fitin[(11+5*nrow(para)):length(fitin)],
      file = paste0(project.path, "/FIT.IN"),
      fill = F,
      sep = "\n",
      append = T)
  
    
    # para docu
    para_g <- data.frame(para, 
                         RMSE = c(rep(NA, nmat-1),goods[1]),
                         R2 = c(rep(NA, nmat-1),goods[2]),
                         NSE = c(rep(NA, nmat-1),goods[3]))
    fwrite(x = as.list(para_g), file = paste0(project.path,"/para_calibration.txt"),
           sep = "\t", col.names = T, append = T)
    
    
  return(para)
}
