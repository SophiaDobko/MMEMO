# Read optimized hydraulic parameter from Fit.out
#
#

project.path = "C:/Users/bauers/data/Hydrus-1D/Projects/spo_cosmic_2_invers"


# Number of soil materials
nmat = 4

# Number of observed data
nobsdata = 178

# Number of iterations
niter = 9

readoptim <- function(project.path, nmat, nobsdata, niter){
  
  # read Fit.out file
  fitout = read.table(paste0(project.path, "/Fit.out"),
                     skip = (15+    # headings
                             10*nmat+  # number soil materials
                             3+
                             nobsdata+   # number observed data
                             6),
                     nrow = niter,       # number of iterations
                     header = T)
  
}
