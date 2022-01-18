# Read optimized hydraulic parameter from Fit.out
#
#

project.path = "C:/Users/bauers/data/Hydrus-1D/Projects/spo_cosmic_2_invers"


# Number of soil materials
nmat = 4

# Number of observed data
nobsdata = 178

# Number of iterations
niter = 11



readoptim <- function(project.path){
  
  # read Fit.out file
  options(warn = -1)
  fitout = data.table::fread(
            input = paste0(project.path, "/Fit.out"),
            skip = " Iteration     SSQ         ALPHA       N     ",
            header = T,
            fill = F,
            blank.lines.skip = F,
            nrows = 21)
  
  return(fitout)
}


