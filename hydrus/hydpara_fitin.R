#' Write new hydraulic parameters in FIT.IN
#'
#' This function writes during hydrus calibration runs new hydraulic soil parameters in fit.in
#' @param project.path directory of the folder containing the hydrus input files
#' @param para dataframe containing the six columns ```thr```, ```ths```, ```Alfa```, ```n```, ```Ks``` and ```l```, containing the parameters for material 1,2...n
#' @return ```selector.in``` with the hydraulic soil parameters for the next calibration run
#' @importFrom data.table fwrite
#' @export



# materials
mat <- rbind(para[1,],
              c(0,0,0,0,0,0),   # parameter optimized? yes/no
              c(0,0,0,0,5,0),   # min values
              c(0,0,0,0,1000,0), # max values
              para[2,],
              c(0,0,0,0,0,0),
              c(0,0,0,0,5,0),
              c(0,0,0,0,1000,0),
             para[3,],
             c(0,0,0,0,0,0),
             c(0,0,0,0,5,0),
             c(0,0,0,0,1000,0),
             para[4,],
             c(0,0,0,0,0,0),
             c(0,0,0,0,5,0),
             c(0,0,0,0,1000,0))


# number of material to be optimized next
n <- 1

# number of material that was optimized last
m <- 1

# parms to be optimized next (yes/no)
#           thr   ths   Alpha  n     Ks   l
poptim <- c( 0 ,   0 ,   1 ,   1 ,   0 ,  0 )



writefit <- function(project.path, para, mat, m, n, poptim){  # works only for model with 4 materials!
  
  # write optimized parameter of the last iteration into para
  # material number must be specified (rownumber of para)!!!
  if(!is.null(fitout$WCR))  para[m,]$thr <- fitout$WCR[nrow(fitout)]
  if(!is.null(fitout$WCS))  para[m,]$ths <- fitout$WCS[nrow(fitout)]
  if(!is.null(fitout$ALPHA)) para[m,]$Alpha <- fitout$ALPHA[nrow(fitout)]
  if(!is.null(fitout$N))   para[m,]$n <- fitout$N[nrow(fitout)]
  if(!is.null(fitout$CONDS)) para[m,]$Ks <- fitout$CONDS[nrow(fitout)]
  if(!is.null(fitout$L))   para[m,]$l <- fitout$L[nrow(fitout)]
  

  # materials
  mat[1,] <- para[1,]
  mat[5,] <- para[2,]
  mat[9,] <- para[3,]
  mat[13,] <- para[4,]
               
  mat[(4*n-2),] <- poptim
  
  # write soil hydraulic parameters in fit.in
  fitin = readLines(con = paste0(project.path, "/FIT.IN"),
                       n = -1L,
                       encoding = "unknown")
  cat(fitin[1:10],
      file = paste0(project.path, "/FIT.IN"),
      fill = F,
      sep = "\n")
  
  fwrite(x = as.list(mat[1:4,]),
         file = paste0(project.path, "/FIT.IN"),
         sep = "\t",
         col.names=T,
         append=T)
  fwrite(x = as.list(mat[5:8,]),
         file = paste0(project.path, "/FIT.IN"),
         sep = "\t",
         col.names=T,
         append=T)
  fwrite(x = as.list(mat[9:12,]),
         file = paste0(project.path, "/FIT.IN"),
         sep = "\t",
         col.names=T,
         append=T)
  fwrite(x = as.list(mat[13:16,]),
         file = paste0(project.path, "/FIT.IN"),
         sep = "\t",
         col.names=T,
         append=T)
  
  cat(fitin[(11+5*nrow(para)):length(fitin)],
      file = paste0(project.path, "/FIT.IN"),
      fill = F,
      sep = "\n",
      append = T)
  
  para <- return(para)
}
