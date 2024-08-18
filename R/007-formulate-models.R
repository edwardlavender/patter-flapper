###########################
###########################
#### formulate-models.R

#### Aims
# 1) Formulate movement and observation models

#### Prerequisites
# 1) Current speed analysis (process-data-spatial-fvcom.R)


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

#### Local pars
yat <- c(0, 0.25, 0.5, 0.75, 1)


###########################
###########################
#### Movement model

###########################
#### Resting

png(here_fig("model-move-step-length-resting.png"), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
x <- seq(0, 100, length.out = 1e5)
y <- dtrunc(x, "cauchy", a = 0, b = 1095, location = 0, scale = 10)
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
xat <- pretty_axis(lim = list(xlim / 120 / 1.75))[[1]]$axis$at
axis(1, at = xlim, lwd.ticks = 0, labels = FALSE, pos = -0.2)
axis(side = 1, at = xat * 120 * 1.75, labels = xat, pos = -0.2)
dev.off()


###########################
#### Active

#### Define parameters
p <- data.frame(k = c(5, 1, 10), 
                theta = c(80, 40, 150), 
                mobility = c(1095, 990, 1125),
                lwd = c(1.5, 1, 1), 
                col = c("black", "red", "blue"), 
                label = c("Best-guess", "Restrictive", "Flexible"))


#### Set up plot
png(here_fig("model-move-step-length-active.png"), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
xlim <- c(0, 1200)

#### Plot probability densities 
for (i in 1:3) {
  x <- seq(0, p$mobility[i] + 1, length.out = 1e5)
  y <- dtrunc(x, "gamma", a = 0, b = p$mobility[i], shape = p$k[i], scale = p$theta[i])
  y <- dtrunc(x, "cauchy", a = 0, b = p$mobility[i], location = p$k[i], scale = p$theta[i])
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
p <- data.frame(alpha = c(4, 3, 5), 
                beta = c(-0.0094, -0.01, -0.009), 
                gamma = c(1500, 1200, 1750),
                lwd = c(1.5, 1, 1), 
                col = c("black", "red", "blue"), 
                label = c("Best-guess", "Restrictive", "Flexible"))


#### Set up plot
png(here_fig("model-obs-acoustic.png"), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
xlim <- c(0, 1750)

#### Plot probability densities 
for (i in 1:3) {
  x <- seq(0, p$gamma[i] + 1, length.out = 1e5)
  y <- ddetlogistic(x, p$alpha[i], p$beta[i], p$gamma[i])
  # y <- y / max(y)
  if (i == 1) {
    print(y[1])
    pretty_plot(x, y,
                xlab = "", ylab = "",
                pretty_axis_args = list(axis = list(x = list(at = c(0, 500, 1000, 1500)), # c(0, 250, 500, 750, 1000, 1250, 1500, 1750)), 
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
}
dev.off()


###########################
###########################
#### Archival observation model

#### Define parameters
# Make plot for seabed depth (mu) = 50 m and 350 m
p <- data.frame(mu = 350, 
                sigma = c(20, 15, 30), 
                deep_depth_eps = c(20, 15, 30),
                lwd = c(1.5, 1, 1), 
                col = c("black", "red", "blue"), 
                label = c("Best-guess", "Restrictive", "Flexible"))


#### Set up plot
png(here_fig(paste0("model-obs-archival-", p$mu[1], "m.png")), 
    height = 4, width = 4, units = "in", res = 600)
set_par()
ylim <- c(p$mu[1] * -1 - max(p$deep_depth_eps) - 20, 0)

#### Plot probability densities 
for (i in 1:3) {
  x <- seq(0, 400, length.out = 1e5)
  y <- dtrunc(x, "norm", a = 0, b = p$mu[i] + p$deep_depth_eps[i], mean = p$mu[i], sd = p$sigma[i])
  # y <- dtrunc(x, "cauchy", a = 0, b = p$mu[i] + p$deep_depth_eps[i], location = p$mu[i], scale = 15)
  y <- y / max(y)
  x <- abs(x) * -1
  # y <- y / max(y)
  if (i == 1) {
    print(y[1])
    pretty_plot(y, x,
                xlab = "", ylab = "",
                pretty_axis_args = list(side = 3:2, 
                                        axis = list(x = list(at = yat), 
                                                    y = list(at = c()))
                ),
                xlim = c(0, 1),
                ylim = ylim,
                type = "n")
    z <- p$mu[i] * -1
    
    # Positive errors & negative errors 
    a <- 0.001
    cols <- rainbow(4, alpha = 0.25)
    fs  <- list(function(x, y) x + y, 
               function(x, y) x - y)
    lapply(1:2, 
           function(i) {
             f <- fs[[i]]
             zbathy <- f(z, 10)
             ztide  <- f(zbathy, 3)
             zsurge <- f(ztide, 2)
             ztag   <- f(ztide, 4.77)
             polygon(c(a, 1, 1, a), c(z, z, zbathy, zbathy), col = cols[1], border = NA)
             polygon(c(a, 1, 1, a), c(zbathy, zbathy, ztide, ztide), col = cols[2], border = NA)
             polygon(c(a, 1, 1, a), c(ztide, ztide, zsurge, zsurge), col = cols[3], border = NA)
             polygon(c(a, 1, 1, a), c(zsurge, zsurge, ztag, ztag), col = cols[4], border = NA)
             if (i == 1) {
               polygon(c(a, 1, 1, a), c(ztag, ztag, 0, 0), col = scales::alpha("lightgrey", 0.25), border = NA)
             }
           })
    
    # Seabed depth
    lines(c(a, 1), c(z, z), lty = 3, col = "dimgrey")
  }
  lines(y, x, lwd = p$lwd[i], col = p$col[i])
}
dev.off()


#### End of code. 
###########################
###########################
