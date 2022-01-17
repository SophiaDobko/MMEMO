#####################################
#### Calibrating hydraulic parameters in Hydrus ####
#####################################


# 1. Take measured or estimated hydraulic parameters from other model or do a quick
#    calibration with data of the profile probe.
# 2. Take the exe. file of the cosmic operator
# 3. Calibrate the cosmic parameter N with the script calibrating_cosmic_N.R
# 4. Copy the input files into the folder of the inverse project

# ...
# write selector
# never save the project inside the Hydrus interface
# docu


setwd("C:/Users/bauers/data/Hydrus-1D/Hydrus-R/Prepare_Hydrus_Input")
project.path = "C:/Users/bauers/data/Hydrus-1D/Projects/spo_cosmic_2_invers"

source("C:/Users/bauers/data/MMEMO/hydrus/hydpara_fitin.R")
source("C:/Users/bauers/data/MMEMO/hydrus/hydpara_selector.R")
source("C:/Users/bauers/data/MMEMO/hydrus/read_optimized_parms.R")


# hydraulic parameter
para <- data.frame(thr = c(0.078, 0.078 ,0.078, 0.078),
                   ths = c(0.37, 0.37, 0.37, 0.37),
                   Alpha = c(0.024, 0.002, 0.063, 0.052),
                   n = c(1.272, 1.683, 1.093, 1.05),
                   Ks = c(500.5, 322.3, 118.6, 630.7),
                   l = c(0.5, 0.5, 0.5, 0.5))


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



writefit(project.path=project.path, para=para, mat=mat, m=m, n=n, poptim=poptim)

hydpara(project.path=project.path, para=para)

readoptim()
