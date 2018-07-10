rm(list = ls())
setwd("C:/Users/tuj53509/Dropbox/docs/Temple/Advanced Statistics for Urban Applications/Assignment4")
library(GISTools)
library(maptools)
library(lm.beta)
library(spdep)


kentucky = readShapeSpatial("ky_counties/ky_counties", proj4string=CRS("+proj=lcc"))
kentucky_census = read.csv("KY.bband.csv")

# Join census file to shapefile
kentucky@data <- data.frame(kentucky@data, kentucky_census[match(kentucky@data[,"NAME10"],kentucky_census[,"Name"]),])

#Linear Model
model01 <- lm(kentucky$adopt05 ~ kentucky$college + kentucky$hhpct_0401 + kentucky$popden2000)
#add standardized betas
model01.std <- lm.beta(model01)
summary(model01.std)
summary(model01$residuals)

pdf(file="f1.adopt05choropleth.pdf")
adopt05.shades <- auto.shading(kentucky$adopt05, cols=brewer.pal(5,"Blues"))
choropleth(kentucky, kentucky$adopt05, shading = adopt05.shades)
title("Figure 1: Kentucky Broadband Adoption")
choro.legend(3751596,4833529, adopt05.shades, title = "proportion of adoption")
map.scale(5021404, 3205804, miles2ft(200), "Miles", 5, 40)
dev.off()

#Lagged Means
#plot queen's case just for comparison's sake
# kentucky_neighbors_queen <- poly2nb(kentucky, queen = TRUE)
# plot(kentucky, main = "DIAG: Queen's Case Weighted Neighbor plot")
# plot(kentucky_neighbors_queen, coordinates(kentucky), add=T, col='blue')


kentucky_neighbors <- poly2nb(kentucky, queen=FALSE)
pdf(file = "f2.weightedneighborplot.pdf")
plot(kentucky, main = "Figure 2: Weighted Neighbor plot of Kentucky Broadband Adoption ")
plot(kentucky_neighbors, coordinates(kentucky), add=T, col='blue')
dev.off()
KY_Neigh_lw <- nb2listw(kentucky_neighbors)
bband <- kentucky$adopt05
bband_lagged_mean <- lag.listw(KY_Neigh_lw, bband)
#create a choropleth map from the lagged means result
par(mar = c(2,1,2,1))
pdf(file="f3.lagged_means_choropleth.pdf")
laggedmeanshades = auto.shading(bband_lagged_mean, cols=brewer.pal(5,"Blues"))
choropleth(kentucky, bband_lagged_mean, shading =laggedmeanshades)
title("Figure 3: Lagged Means plot of Kentucky Broadband Adoption")
choro.legend(3751596,4833529, laggedmeanshades, title = "Lagged Means")
dev.off()

#create a Lagged Means Plot (Moran's Plot)
pdf(file = "f4.moransplot.pdf")
moran.plot(kentucky$adopt05, KY_Neigh_lw)
title("Figure 4: Moran's Plot")
#this closes the file handle
dev.off()

#Run a Moran's I test
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(KY_Neigh_lw)
#approximate test statistic using normal distribution
moran.test(kentucky$adopt05, KY_Neigh_lw)
#calculate the test statistic using 10,000 random trials
moran.mc(kentucky$adopt05, KY_Neigh_lw, 10000)

#run the SAR model
model01.lag <- lagsarlm(kentucky$adopt05 ~ kentucky$college + kentucky$hhpct_0401 + kentucky$popden2000, data = kentucky, KY_Neigh_lw)
summary(model01.lag)
anova(model01.lag, model01)

