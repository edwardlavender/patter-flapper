###########################
###########################
#### prepare-patter-param-plot.R

#### Aims
# 1) Visualise observation and movement models

#### Prerequisites
# 1) Define models


###########################
###########################
#### Set up 

#### Wipe workspace 
rm(list = ls())
try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
library(dv)

#### Load data
src()


###########################
###########################
#### Plot models 

#### Plot observational model (detection)
# TO DO

#### Plot observational model (depths) 
plot(0, xlim = c(-50, 25), ylim = c(0, 1), 
     axes = FALSE, xlab = "", ylab = "",
     type = "n")
add_depth_error_model(150)
axis(side = 1, seq(-50, 25, by = 25), pos = 0)
axis(side = 2, pos = -50, las = TRUE)
mtext(side = 1, "Depth (m)", line = 2)
mtext(side = 2, "Weight", line = 2)

#### Plot movement model
# TO DO


#### End of code. 
###########################
###########################