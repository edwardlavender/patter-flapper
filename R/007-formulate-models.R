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
library(prettyGraphics)
library(fvcom.tbx)
library(truncdist)


###########################
###########################
#### Movement model

###########################
#### Resting


x <- seq(0, 100, length.out = 1e5)
y <- dtrunc(x, "cauchy", a = 0, b = 1095, location = 0, scale = 10)
y <- y / max(y)
pretty_plot(x, y,
            xlab = "", ylab = "",
            xlim = range(x),
            type = "l")


###########################
#### Active

# TO DO
# Save figs
# Add arrows with guides for other species

#### Define parameters
p <- data.frame(k = c(5, 1, 10), 
                theta = c(80, 40, 150), 
                mobility = c(1095, 990, 1125),
                lwd = c(2, 1, 1), 
                col = c("black", "red", "blue"))


#### Set up plot
# png(here_fig("model-move.png"),
#     height = 3, width = 6, units = "in", res = 600)
# pp <- set_par(mfrow = c(1, 2), mar = c(1.5, 1.5, 1.5, 1.5), oma = c(2, 3, 1, 1))
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
              xlim = xlim,
              type = "n")
  }
  lines(x, y, lwd = p$lwd[i], col = p$col[i])
}

# Add axis for DLs-1
xat <- pretty_axis(lim = list(xlim / 120 / 1.75))[[1]]$axis$at
axis(1, at = xlim, lwd.ticks = 0, labels = FALSE, pos = -0.1)
axis(side = 1, at = xat * 120 * 1.75, labels = xat, pos = -0.1)

# mtext(side = 1, "Step length", line = 2)
# mtext(side = 2, "Density", line = 3.5)
# dev.off()

stop()

###########################
#### Turning angles

x  <- seq(-pi*1.1, pi*1.1, length.out = 1e5)
y   <- dunif(x, -pi, pi)
pretty_plot(x, y,
            xlab = "", ylab = "",
            # pretty_axis_args = paa,
            type = "l")

lines(x, y)
# mtext(side = 1, "Turning angle (rad)", line = 2)
# mtext(side = 2, "Density", line = 3.5)
# dev.off()


###########################
###########################
#### Acoustic observation model 

#### TO DO (revise)

#### Observation models
png(here_fig("model-obs.png"),
    height = 3, width = 6, units = "in", res = 600)
pp <- set_par(mfrow = c(1, 2), mar = c(1.5, 1.5, 1.5, 1.5), oma = c(2, 3, 1, 1))
# Acoustic observations (0, 1)
a <- moorings$receiver_alpha[1]
b <- moorings$receiver_beta[1]
g <- moorings$receiver_gamma[1]
x <- seq(1, 1000, by = 1)
y <- dbinom(1, size = 1, prob = ifelse(x <= g, plogis(a * b * x), 0))
pretty_plot(x, y,
            pretty_axis_args = paa,
            xlab = "", ylab = "",
            type = "l")
add_dbn(x, y)
mtext(side = 1, "Distance (m)", line = 2)
mtext(side = 2, "Density", line = 3.5)
# Archival observations
x <- seq(-25, 25, by = 0.1)
y <- dunif(x, -20, 20)
axis_ls <- pretty_plot(x, y,
                       ylim = c(0, 0.04),
                       pretty_axis_args = list(pretty = list(n = 2)),
                       xlab = "", ylab = "",
                       type = "n")
# lines(c(0, 0), c(0, max(y)), lty = 3)
px <- par(xpd = NA)
arrows(0, 0, 0, 0.03, length = 0.1)
par(px)
add_dbn(x, y, border = "black", lty = 3)
mtext(side = 1, "Distance (m)", line = 2)
par(pp)
dev.off()


###########################
###########################
#### Archival observation model



#### End of code. 
###########################
###########################
