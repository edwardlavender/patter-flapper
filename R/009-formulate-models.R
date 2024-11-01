###########################
###########################
#### formulate-models.R

#### Aims
# 1) Formulate movement and observation models

#### Prerequisites
# 1) Current speed analysis (process-data-spatial-fvcom.R)
# 2) Parameter choices are informed by run-patter-trials.R


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())
# try(pacman::p_unload("all"), silent = TRUE)
dv::clear()

#### Essential packages
dv::src()
library(prettyGraphics)
library(fvcom.tbx)
library(truncdist)

#### Load data
map <- terra::rast(here_data("spatial", "bathy.tif"))

#### Local pars
# Sensitivity thresholds
# * Consider inflating/deflating parameters by a constant amount for comparability 
# * Uncertainty in the movement and the detection probability model is not equal
# * Different degrees of uncertainty are required
# * These degrees are for the movement model 
psensitivity <- 0.5
pinflate     <- 1 + psensitivity
pdeflate     <- 1 - psensitivity
tsensitivity <- 0.1
tinflate     <- 1 + tsensitivity
tdeflate     <- 1 - tsensitivity 
# Graphics
yat         <- c(0, 0.25, 0.5, 0.75, 1)
lwds        <- c(1.5, 1, 1)
cols        <- c("black", "red", "blue")


###########################
###########################
#### Movement model

###########################
#### Resting

#### Define parameters
pmovement_rest <- data.table(k1 = c(0, 0, 0), 
                             theta1 = c(5, 5, 5), 
                             mobility = c(1095, 990, 1125))
p       <- copy(pmovement_rest)
p$label <- c("Best-guess", "Restrictive", "Flexible")

png(here_fig("model-move-step-length-resting.png"), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
x <- seq(0, 100, length.out = 1e5)
y <- dtrunc(x, "cauchy", a = 0, b = p$mobility[1], location = p$k1[1], scale = p$theta1[1])
y <- y / max(y)
xlim <- range(x)
pretty_plot(x, y,
            xlab = "", ylab = "",
            pretty_axis_args = list(pretty = list(n = 4), 
                                    axis = list(x = list(), 
                                                y = list(at = yat))
                                    ),
            xlim = xlim,
            type = "l")
if (FALSE) {
  # Compare new values to current model 
  x <- seq(0, 100, length.out = 1e5)
  y <- dtrunc(x, "cauchy", a = 0, b = 1095, location = 0, scale = 5)
  y <- y / max(y)
  lines(x, y, col = "red")
}
xat <- pretty_axis(lim = list(xlim / 120 / 1.75))[[1]]$axis$at
axis(1, at = xlim, lwd.ticks = 0, labels = FALSE, pos = -0.2)
axis(side = 1, at = xat * 120 * 1.75, labels = xat, pos = -0.2)
dev.off()


###########################
#### Active

#### Define parameters
# pmovement_active <-  data.table(k2 = c(5, 1, 10), 
#                                 theta2 = c(100, 60, 150), 
#                                 mobility = c(1095, 990, 1125))
pmovement_active <-  data.table(k2 = c(5, 5 * pdeflate, 5 * pinflate), 
                                theta2 = c(100, 100 * pdeflate, 100 * pinflate), 
                                mobility = c(1095, 1095 * tdeflate, 1095 * tinflate))
p       <- copy(pmovement_active)
p$label <- c("Best-guess", "Restrictive", "Flexible")
p$lwd   <- lwds
p$col   <- cols

#### Set up plot
png(here_fig("model-move-step-length-active.png"), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
xlim <- c(0, 1200)

#### Plot probability densities 
for (i in 1:3) {
  x <- seq(0, p$mobility[i] + 1, length.out = 1e5)
  # y <- dtrunc(x, "gamma", a = 0, b = p$mobility[i], shape = p$k2[i], scale = p$theta2[i])
  y <- dtrunc(x, "cauchy", a = 0, b = p$mobility[i], location = p$k2[i], scale = p$theta2[i])
  y <- y / max(y)
  if (i == 1) {
  pretty_plot(x, y,
              xlab = "", ylab = "",
              pretty_axis_args = list(pretty = list(n = 4), 
                                      axis = list(x = list(at = c(0, 250, 500, 750, 1000)), 
                                                  y = list(at = yat))
              ),
              xlim = xlim,
              type = "n")
  }
  lines(x, y, lwd = p$lwd[i], col = p$col[i])
  arrows(x0 = p$mobility[i], 
         x1 = p$mobility[i], 
         y0 = 0.075, 
         y1 = 0.025, 
         length = 0.05, 
         col = p$col[i])
}

if (FALSE) {
  # Compare new values to current model 
  x <- seq(0, p$mobility[i] + 1, length.out = 1e5)
  y <- dtrunc(x, "cauchy", a = 0, b = 1095, location = 5, scale = 100)
  y <- y / max(y)
  lines(x, y, col = "purple")
}

# Add legend
legend(500, 0.99,
       legend = p$label,
       lwd = p$lwd, 
       col = p$col, 
       box.lty = 3)

# Add axis for DLs-1
xat <- pretty_axis(lim = list(xlim / 120 / 1.75))[[1]]$axis$at
axis(1, at = xlim, lwd.ticks = 0, labels = FALSE, pos = -0.2)
axis(side = 1, at = xat * 120 * 1.75, labels = xat, pos = -0.2)
# mtext(side = 1, "Step length", line = 2)
# mtext(side = 2, "Density", line = 3.5)
dev.off()


###########################
#### Turning angles

png(here_fig("model-move-turning-angle.png"), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
x <- seq(-pi*1.1, pi*1.1, length.out = 1e5)
y <- dunif(x, -pi, pi)
y <- y / max(y)
plot(x, y,
     xlab = "", ylab = "",
     type = "l", 
     axes = FALSE)
axis(side = 1, at = c(-pi * 1.1, pi * 1.1), labels = FALSE, lwd.tick = 0, pos = 0)
axis(side = 1, 
     at = c(-pi, -pi/2, 0, pi/2, pi), 
     labels = c(expression(-pi), expression(-pi/2), expression(0), expression(pi/2), expression(pi)), 
     pos = 0)
axis(side = 2, prettyGraphics:::add_lagging_point_zero(yat), las = TRUE, pos = -pi*1.1)

# mtext(side = 1, "Turning angle (rad)", line = 2)
# mtext(side = 2, "Density", line = 3.5)
dev.off()


###########################
###########################
#### Acoustic observation model 

#### Define parameters
# pdetection <- data.table(receiver_alpha = c(4, 3, 5), 
#                          receiver_beta = c(-0.0094, -0.01, -0.009), 
#                          receiver_gamma = c(1750, 1500, 2000))
pinflate <- 1.25
pdeflate <- 0.75
pdetection <- data.table(receiver_alpha = c(4, 4 * pinflate, 4 * pdeflate), 
                         receiver_beta = c(-0.0094, -0.0094 * pdeflate, -0.0094 * pinflate), 
                         receiver_gamma = c(3000, 3000 * pinflate, 3000 * pdeflate))
p       <- copy(pdetection)
p$label <- c("Best-guess", "Restrictive", "Flexible")
p$lwd   <- lwds
p$col   <- cols

#### Set up plot
png(here_fig("model-obs-acoustic.png"), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
xlim <- c(0, 4000)

#### Plot probability densities 
for (i in 1:3) {
  x <- seq(0, p$receiver_gamma[i] + 1, length.out = 1e5)
  y <- ddetlogistic(x, p$receiver_alpha[i], p$receiver_beta[i], p$receiver_gamma[i])
  # y <- y / max(y)
  if (i == 1) {
    print(y[1])
    pretty_plot(x, y,
                xlab = "", ylab = "",
                pretty_axis_args = list(axis = list(x = list(at = c(0, 1000, 2000, 3000, 4000)), 
                                                    y = list(at = yat))
                ),
                xlim = xlim,
                ylim = c(0, 1),
                type = "n")
    axis(side = 1, at = xlim, lwd.ticks = 0, labels = FALSE, pos = 0)
    polygon(c(400, 450, 450, 400), c(0, 0, 1, 1), col = scales::alpha("grey", 0.75), border = NA)
    lines(xlim, c(0.5, 0.5), lty = 3, col = "dimgrey")
  }
  lines(x, y, lwd = p$lwd[i], col = p$col[i])
  arrows(x0 = p$receiver_gamma[i], 
         x1 = p$receiver_gamma[i], 
         y0 = 0.05, 
         y1 = 0, 
         length = 0.05, 
         col = p$col[i])
  # rug(p$receiver_gamma[i], pos = -0.1, lwd = 2, col = p$col[i])
}
dev.off()


###########################
###########################
#### Archival observation model

#### Define parameters
# Make plot for seabed depth (mu) = 50 m and 350 m
pdepth <- data.frame(mu = 200, 
                     depth_sigma = c(100, 100 * 0.5, 100 * 1.5), 
                     depth_deep_eps = c(350, 350, 350))
p       <- copy(pdepth)
p$label <- c("Best-guess", "Restrictive", "Flexible")
p$lwd   <- lwds
p$col   <- cols

#### Set up plot
png(here_fig(paste0("model-obs-archival-", p$mu[1], "m.png")), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
ylim <- c(-350, 0) # c(p$mu[1] * -1 - max(p$depth_deep_eps) - 20, 0)

#### Plot probability densities 
for (i in 1:3) {
  
  ##### Plot likelihood profile
  x <- seq(0, 350, length.out = 1e5)
  y <- dtrunc(x, "norm", a = 0, b = 350, mean = p$mu[i], sd = p$depth_sigma[i])
  # y <- dtrunc(x, "cauchy", a = 0, b = p$mu[i] + p$depth_deep_eps[i], location = p$mu[i], scale = 15)
  y <- y / max(y)
  x <- abs(x) * -1
  # y <- y / max(y)
  if (i == 1) {
    print(y[1])
    pretty_plot(y, x,
                xlab = "", ylab = "",
                pretty_axis_args = list(side = 3:2, 
                                        axis = list(x = list(at = yat), 
                                                    y = list(at = c(0, -100, -200, -300)))
                ),
                xlim = c(0, 1),
                ylim = ylim,
                type = "n")
    
    #### Visualise uncertainty components 
    # Define helper function 
    add_uncertainty <- function(mu, sds, cols, alphas = seq(0.8, 0.2, length.out = length(sds))) {
      # Iteratively add polygons 
      for (j in length(sds):1) {

        # Accumulate SDs and plot wider polygon (smaller polygons are added later)
        sigma <- sum(sds[1:j])
        x <- seq(-1000, 1000, length.out = 1e5)
        y <- dtrunc(x, "norm", a = 0, b = 350, mean = mu, sd = sigma)
        y <- y / max(y)
        x <- abs(x) * -1
        polygon(x = c(y, rev(y)), y = c(x, rep(0, length(x))), 
                col = scales::alpha(cols[j], alphas[j]), 
                border = NA)
      }
    }
    # Define uncertainties 
    zbathy <- 10 / 2        # bathy data are ± 10 m (95 % of data is within ± 2 SD of mean = ± 10/2 m of mean)
    ztide  <- 3 / 2         # tides are ± 3 m
    zsurge <- 2 / 2         # surges are ± 2 m
    ztag   <- 4.77 / 2      # tag accuracy is ± 4.77 m
    zagg   <- 75            # near-max SD in aggregated bathymetric grid cells (sufficient to cover majority of cells)
    # Collate uncertainties & colours
    sds  <- c(zbathy, ztide, zsurge, ztag, zagg) 
    cols <- rainbow(5, alpha = 0.25)
    sum(sds)
    # Add rect
    rect(xleft = 0, xright = 1, ybottom = p$mu[i] * -1, ytop = 0, col = scales::alpha("grey", 0.2), border = NA)
    add_uncertainty(mu = p$mu[i], sds = sds, cols = cols)
    # Seabed depth
    a <- 0.001; z <- p$mu[i] * -1
    lines(c(a, 1), c(z, z), lty = 3, col = "dimgrey")
  }
  lines(y, x, lwd = p$lwd[i], col = p$col[i])
  arrows(x0 = 0.05, 
         x1 = 0, 
         y0 = (p$mu[i] + p$depth_deep_eps[i]) * -1, 
         y1 = (p$mu[i] + p$depth_deep_eps[i]) * -1, 
         length = 0.05, 
         col = p$col[i])
}
dev.off()


###########################
###########################
#### Save parameter datasets

#### All parameters
pars <- list(pmovement = cbind(pmovement_rest[, .(k1, theta1)], pmovement_active),
             pdetection = pdetection, 
             pdepth = pdepth)
qs::qsave(pars, here_data("input", "pars.qs"))

#### Validity maps for mobility parameters (~6 mins)
pp <- par(mfrow = c(1, 3))
cl_lapply(pars$pmovement$mobility, function(mob) {
  vmap <- spatVmap(.map = map, .mobility = mob, .plot = TRUE)
  terra::writeRaster(vmap, 
                     here_data("spatial", glue("vmap-{mob}.tif")), 
                     overwrite = TRUE)
  
})
par(pp)


#### End of code. 
###########################
###########################
