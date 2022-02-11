#####################################
#### Calibrating hydraulic parameters in Hydrus ####
#####################################


# 1. Take measured or estimated hydraulic parameters from other model or do a quick
#    calibration with data of the profile probe.
# 2. Take the exe. file of the cosmic operator
# 3. Calibrate the cosmic parameter N with the script calibrating_cosmic_N.R
#    and a forward modelling
# 4. Copy the input files into the folder of the inverse project
# 5. Run this script
# 6. Start Hydrus model in the Hydrus interface
# 7. Run this skript again from "read optimized parameter"
#    n, m and poptim must be specified manually
# 8. Start Hydrus model again in the Hydrus interface and so on...
#    never save the project inside the Hydrus interface!


library(hydrusR)
library(hydroGOF)

setwd("C:/Users/bauers/data/Hydrus-1D/Hydrus-R/Prepare_Hydrus_Input")
#project.path = "C:/Users/bauers/data/Hydrus-1D/Projects/spo_cosmic_2"

source("C:/Users/bauers/data/MMEMO/hydrus/hydpara_fitin.R")
source("C:/Users/bauers/data/MMEMO/hydrus/hydpara_selector.R")
source("C:/Users/bauers/data/MMEMO/hydrus/read_optimized_parms.R")
source("C:/Users/bauers/data/MMEMO/hydrus/goodness_function_plots.R")


# hydraulic parameter
#para <- data.frame(thr = c(0.078, 0.078 ,0.078, 0.078),
#                   ths = c(0.37, 0.37, 0.37, 0.37),
#                   Alpha = c(0.024, 0.002, 0.063, 0.052),
#                   n = c(1.272, 1.683, 1.093, 1.05),
#                   Ks = c(500.5, 322.3, 118.6, 630.7),
#                   l = c(0.5, 0.5, 0.5, 0.5))

# number of materials
nmat <- 3

# materials
mat <- rbind(para[1,],
             c(0,0,0,0,0,0),   # parameter optimized? yes/no
             c(0,0,0,0,1,0),   # min values
             c(0,.5,0,0,1000,0), # max values
             para[2,],
             c(0,0,0,0,0,0),
             c(0,0,0,0,1,0),
             c(0,.45,0,0,1000,0),
             para[3,],
             c(0,0,0,0,0,0),
             c(0,0,0,0,1,0),
             c(0,.45,0,0,1000,0) ) #,
#             para[4,],
#             c(0,0,0,0,0,0),
#             c(0,0,0,0,1,0),
#             c(0,0,0,0,1000,0))


# read optimized parameter
(fitout <- readoptim(project.path = project.path))

# goodness function and plot model results
(goods <- goodness(project.path = project.path))

# number of material to be optimized next
n <- 2
# number of material that was optimized last
m <- 1
# parms to be optimized next (yes/no)
#           thr   ths   Alpha  n     Ks   l
poptim <- c( 0 ,   0 ,   1 ,   1 ,   1 ,  0 )

#fitout <- NULL  # take this to set parameter to initial values or if para should not be changed
#goods <- rep(NA, nmat) # take this to set parameter to initial values

# write new input files
para <- writefit(project.path=project.path, para=para, mat=mat, nmat=nmat, m=m, n=n, poptim=poptim)
hydpara(project.path=project.path, para=para)


