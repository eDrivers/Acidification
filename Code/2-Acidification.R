# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  LIBRARIES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(sp)
library(magrittr)
library(tidyverse)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                    DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('./data/rawData/bathy.RData')
load('./data/rawData/aragonite.RData')
load('./data/rawData/doxy.RData')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                         DIVIDE SOUTHERN AND NORTHERN GULF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
south <- (bathy$Y < 550000 & bathy$DEPH < 350 & bathy$X > 120000 & bathy$X < 650000)
south[bathy$X > 550000 & bathy$Y > 490000] <- FALSE
north <- !south
southGrid <- bathy[south, ]
northGrid <- bathy[north, ]

# Viz
png('./figures/delimStL.png', width = 1280, height = 920, res = 300, pointsize = 6)
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot(st_geometry(bathy), cex = .75)
plot(st_geometry(southGrid), col = '#fa9b7e', add = T, cex = .75)
plot(st_geometry(northGrid), col = '#d2fa7e', add = T, cex = .75)
dev.off()

# Select aragonite observations belonging to north vs south
aragSouth <- arag %>%
             st_intersects(southGrid, .) %>%
             unlist() %>%
             arag[., ]

aragNorth <- arag %>%
             st_intersects(northGrid, .) %>%
             unlist() %>%
             arag[., ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                DIVIDE NORTHERN GULF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Seperate channels from shallower areas
id <- which((aragNorth$Aragonite > .8 & aragNorth$DOXY < 4 |
      aragNorth$Aragonite > .69 & aragNorth$DOXY < 3 ))

# Linear models between omega aragonite and dissolved oxygen
mod1 <- lm(aragNorth$Aragonite[id] ~ aragNorth$DOXY[id])
mod2 <- lm(aragNorth$Aragonite[-id] ~ aragNorth$DOXY[-id])
r1 <- round(summary(mod1)$adj.r.squared, 2)
r2 <- round(summary(mod2)$adj.r.squared, 2)

# Visualize
col1 <- '#66c1c7'
col2 <- '#e7a867'

png('./figures/northArag-O2.png', width = 1280, height = 920, res = 300, pointsize = 6)
plot(aragNorth$Aragonite ~ aragNorth$DOXY, xlab = 'Dissolved oxygen', ylab = 'Omega Aragonite')
points(aragNorth$Aragonite[id] ~ aragNorth$DOXY[id], col = col1, pch = 20)
points(aragNorth$Aragonite[-id] ~ aragNorth$DOXY[-id], col = col2, pch = 20)
abline(mod1, col = col1)
abline(mod2, col = col2)
text(x = 2.5, y = 1.05, labels = 'Channels', col = col1, font = 2, cex = 1.25, adj = c(0,0))
text(x = 2.5, y = 1.025, labels = expression(R^2 ~ '= 0.96'), col = col1, font = 2, adj = c(0,0))
text(x = 5.5, y = .95, labels = 'Shallows', col = col2, font = 2, cex = 1.25, adj = c(0,0))
text(x = 5.5, y = .925, labels = expression(R^2 ~ '= 0.95'), col = col2, font = 2, adj = c(0,0))

# par(new = TRUE)
# par(fig = c(.49,1,.15,.65))
# plot(st_geometry(bathy), border = '#00000055')
# plot(st_geometry(aragNorth), add = T)
# plot(st_geometry(aragNorth[id,]), add = T, pch = 20, col = col1, cex = .75)
# plot(st_geometry(aragNorth[-id,]), add = T, pch = 20, col = col2, cex = .75)
dev.off()
# There indeed seem to be two different processes at play here.
# Let's see if we can use them to predict omega aragonite in the Northern Gulf using dissolved oxygen
# First, we need the dissolved oxygen data from the PMZA that we use to generate the hypoxia layer


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           INTERPOLATION DISSOLVED OXYGEN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Years
ys <- names(doxy)

# Get centroid of grid cells
coords <- st_centroid(bathy)

# Transform to sp object
coords <- as(coords, 'Spatial')
for(i in 1:length(doxy)) doxy[[i]] <- as(doxy[[i]], 'Spatial')

# Kriging Formula
varioForm <- as.formula("DOXY ~ DEPH")

# Interpolation
kriged <- vector('list', length(doxy))
for(i in 1:length(doxy)) {
  # Variogram
  varioFit <- automap::autofitVariogram(varioForm, doxy[[i]])
  plot(varioFit)

  # Krigging model
  kriged[[i]] <- gstat::krige(varioForm, doxy[[i]], coords, model = varioFit$var_model) %>%
                 st_as_sf(kriged)
}

# Single object
doxy <- bathy
for(i in 1:length(ys)) doxy <- cbind(doxy, kriged[[i]][, 1, drop = T])
colnames(doxy)[grep('kriged', colnames(doxy))] <- ys

# Annual mean of O2 saturation in deep waters
doxy$doxy <- rowMeans(doxy[, ys, drop = T])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                     PREDICT OMEGA ARAGONITE ~ DISSOLVED OXYGEN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Divide the dissolved oxygen dataset between deep channels and shallower areas
# Keep only Northern gulf
doxy <- doxy[north, ]

# Divide channels from shallower areas with depth > 350 meters
# This is debatable, but it is what provides the best estimates combined with
# knowledge expert so far
channel <- doxy$DEPH > 350
shallow <- !channel

# Visualize division
par(mfrow = c(1,1))
plot(st_geometry(doxy))
plot(st_geometry(doxy)[channel], col = '#ffffff', add = T)
plot(st_geometry(doxy)[shallow], col = '#000000', add = T)

# Predict aragonite ~ dissolved oxygen as a function of depth
doxy$ARAG <- 0
doxy$ARAG[channel] <- mod1$coefficients[1] + mod1$coefficients[2] * doxy$doxy[channel]
doxy$ARAG[shallow] <- mod2$coefficients[1] + mod2$coefficients[2] * doxy$doxy[shallow]

# Check in grid
envNorth <- bathy[north, ]
envNorth$ARAG <- doxy$ARAG
# envNorth$ARAG2 <- doxy$ARAG2

png('./figures/north.png', width = 1200, height = 950, res = 300, pointsize = 6)
plot(envNorth[, 'ARAG'], pch = 20, border = 'transparent')
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   SOUTHERN GULF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Southern Gulf
# Get centroid of grid cells
southGrid <- st_centroid(southGrid)

# Transform to sp object
southGrid <- as(southGrid, 'Spatial')
aragSouth <- as(aragSouth, 'Spatial')

# Kriging dataframe
# Formula
varioForm1 <- as.formula("Aragonite ~ 1")

# Variograms
varioFit1 <- automap::autofitVariogram(varioForm1, aragSouth)
plot(varioFit1)

# Krigging model
kriged <- gstat::krige(varioForm1, aragSouth, southGrid, model = varioFit1$var_model) %>%
          st_as_sf(kriged)

# Check in grid
envSouth <- bathy[south, ]
envSouth$ARAG <- kriged[, 1, drop = T]

png('./figures/south.png', width = 1200, height = 1000, res = 300, pointsize = 6)
plot(envSouth[, 'ARAG'], pch = 20, border = 'transparent')
dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                             COMBINE NORTH AND SOUTH
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aragPred <- rbind(envNorth[, 'ARAG'], envSouth[, 'ARAG'])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            PREDICTIONS VS OBSERVATIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Intersect aragonite dataset with the predictions
arag$pred <- arag %>%
             st_intersects(aragPred) %>%
             unlist() %>%
             aragPred$ARAG[.]

png('./figures/valid.png', width = 1500, height = 1000, res = 300, pointsize = 6)
par(mfrow = c(2,2))
# Compare observed vs predicted
modPred <- lm(arag$Aragonite ~ arag$pred)
summary(modPred)
arag$res <- resid(modPred)
par(mar = c(4,4,2,2))
plot(arag$Aragonite ~ arag$pred, pch = 20, col = '#33aff599',
     main = 'Observed vs predicted omega aragonite',
     ylab = 'Observed',
     xlab = 'Predicted')
abline(0, 1, col = '#00000055', lwd = 2)
text(x = .5, y = 1.6, adj = c(0,1), labels = expression(R^2 ~ '= 0.91'))

# Look at residuals
plot(arag$Aragonite, arag$res, pch = 20, col = '#33aff599',
     ylab="Residuals", xlab="Omega aragonite",
     main="Residuals")
abline(0, 0, col = '#00000055', lwd = 2)

# Spatial distribution of residuals
plot(st_geometry(arag), pch = 20, cex = abs(arag$res) * 10, col = '#33aff599')
text(x = mean(par('usr')[1:2]), par('usr')[4]-100000, 'Residual size', font = 2)
dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 ACIDIFICATION INDEX
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sig <- function(x, L = 3, Q = 1, k = 2) {
  dat <- -L / (.99 + (Q * exp(-k * x))) + L
  dat[dat < 0] <- 0
  dat[dat > 1] <- 1
  dat
}
x <- seq(0, 4, by = .01)
png('./figures/acidFunction.png', width = 1000, height = 800, res = 300, pointsize = 6)
par(mfrow = c(1,1), mar = c(4,4,2,2))
plot(x, sig(x), pch = 20, cex = .5, type = 'l', lwd = 1, xlab = 'Omega aragonite',
     ylab = 'Acidification index', frame = F, ylim = c(0,1), main = 'Acidification index')
dev.off()

# Acidification
aragPred$ACID <- sig(aragPred$ARAG)

# Remove cells with values = 0
acidification <- aragPred[aragPred$ACID > 0, 'ACID']

# Change name
acidification <- acidification %>%
                 rename(Acidification = ACID)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                  EXPORT DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Export object as .RData
save(acidification, file = './Data/Driver/Acidification.RData')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                 VISUALIZE DATA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png('./Figures/Acidification.png', width = 1280, height = 1000, res = 200, pointsize = 6)
plot(acidification[, 'Acidification'], border = 'transparent')
dev.off()
