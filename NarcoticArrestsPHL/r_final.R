rm(list = ls())
setwd("C:/Users/tuj53509/Dropbox/docs/Temple/Advanced Statistics for Urban Applications/Final Project")

#load reqiured libraries
library(GISTools)
library(maptools)
library(spgwr)
library(spdep)
library(lm.beta)
library(pastecs)

#load data
philly <- readShapeSpatial("data/narc_join")
philly_acs <- read.csv("R11739357_SL140.csv")

#join CSV data to philly DF
philly@data <- data.frame(philly@data, philly_acs[match(philly@data[,"GEOID10"],philly_acs[,"Geo_FIPS"]),])

#SE_T083_001:    per capita income
#SE_T131_001 pop vet status 18 and over
#SE_T095_003 vacant houses

#SE_T033_006 unemployment pop over 16 - not used
#SE_T140_003 service job - not used
#SE_T145_002 no health insurance coverage - not significant
#SE_T098_001 median year structure built - unused

#normalize the count of narcotics incidents
philly$popinthousands <- philly$SE_T001_001 / 1000
philly$narcnorm2 <- philly$Count_ / philly$popinthousands
#unemployment
philly$unempnorm <- philly$SE_T033_006 / philly$SE_T033_001
#veteran
philly$vetnorm <- philly$SE_T131_002  / philly$SE_T131_001
#vacant houses
philly$vacantnorm <- philly$SE_T095_003 / philly$SE_T095_001
#service jobs
philly$servnorm <- philly$SE_T140_003 / philly$SE_T140_001
  #remove bad values from data created by 0 population tracts
na.inf <- function (x) {
  x[is.infinite(x)] <- 0
  return(x)
}
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

philly$narcnorm2 <- na.zero(philly$narcnorm2)
philly$unempnorm <- na.zero(philly$unempnorm)
philly$vetnorm <- na.zero(philly$vetnorm)
philly$vacantnorm <- na.zero(philly$vacantnorm)
philly$servnorm <- na.zero(philly$servnorm)
philly$incomefix <- na.zero(philly$SE_T083_001)
philly$narcnorm2 <- na.inf(philly$narcnorm2)
philly$unempnorm <- na.inf(philly$unempnorm)
philly$vetnorm <- na.inf(philly$vetnorm)
philly$vacantnorm <- na.inf(philly$vacantnorm)
philly$servnorm <- na.inf(philly$servnorm)


vals <- cbind(philly$narcnorm2, philly$vetnorm, philly$incomefix, philly$servnorm)
stat.desc(vals)
summary(vals)

pdf(file = "f1.narcoticsarrestschoropleth.pdf")
narcnorm.shades <- auto.shading(philly$narcnorm2, cols=brewer.pal(5,"Greens"))
choropleth(philly, philly$narcnorm2, shading = narcnorm.shades)
title("Figure 1: Narcotics Arrests, 2012-2016")
choro.legend(2725857,238144, narcnorm.shades, title = "arrests per thousand")
dev.off()

#create our linear model
# lm01 <- lm(philly$Count_ ~ philly$SE_T033_001)
# lm02 <- lm(philly$Count_ ~ philly$SE_T033_001 + philly$SE_T083_001)
# lm03 <- lm(philly$narcnorm2 ~ philly$SE_T033_001 + philly$SE_T083_001 + philly$SE_T131_001 + philly$SE_T095_003 + philly$SE_T140_003)
#lm04 <- lm(philly$narcnorm2 ~  philly$vetnorm + philly$incomefix + philly$vacantnorm )
lm04 <- lm(philly$narcnorm2 ~  philly$vetnorm + philly$incomefix + philly$servnorm )
summary(lm04)
#add standardized betas
lm04.std <- lm.beta(lm04)
summary(lm04.std)
#dianostic plots for lm04
pdf(file = "f2.lm04diag.pdf")
par(mar=c(2,2,2,2),mfrow=c(2,2))
plot(lm04)
par(mfrow=c(1,1))
dev.off()

#lm04 residuals
lm04residshades = auto.shading(lm04$residuals, cols=brewer.pal(5,"Greens"))
pdf(file = "f3.lm04residchoropleth.pdf")
choropleth(philly, lm04$residuals, shading = lm04residshades)
title("Figure 3: Residuals from linear model lm04")
choro.legend(2729559,251723.8, lm04residshades, title = "residuals")
dev.off()

#create neighbor list and plot
philly_neighbors <- poly2nb(philly)
pdf(file = "f4.neighbormap.pdf")
plot(philly, main = "Figure 4: Neighbor plot of Philadelphia Census Tracts")
plot(philly_neighbors, coordinates(philly), add=T, col='blue')
dev.off()

#create lagged means
philly_neigh_lw <- nb2listw(philly_neighbors)
narclag <- lag.listw(philly_neigh_lw, philly$narcnorm2)
#create lagged means plot
laggedmeanshades = auto.shading(narclag, cols=brewer.pal(5,"Greens"))
pdf(file = "f5.laggedmeanchoropleth.pdf")
choropleth(philly, narclag, shading =laggedmeanshades)
title("Figure 5: Lagged Means plot of Philadelphia Narcotics Arrests")
choro.legend(2729559,251723.8, laggedmeanshades, title = "Lagged Means")
dev.off()

#Moran's I
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(philly_neigh_lw)
#approximate test statistic using normal distribution
moran.test(philly$narcnorm2, philly_neigh_lw)
#calculate the test statistic using 10,000 random trials
moran.mc(philly$narcnorm2, philly_neigh_lw, 10000)
#Moran's Plot
pdf(file = "f6.moransplotnorm.pdf")
moran.plot(philly$narcnorm2, philly_neigh_lw)
title("Figure 6: Moran's Plot: Normalized Arrests")
dev.off()
pdf(file = "f7.moransplotraw.pdf")
moran.plot(philly$Count_, philly_neigh_lw)
title("Figure 7: Moran's Plot: Raw Arrest Count")
dev.off()

#GWR
philly.bw <- gwr.sel(philly$narcnorm2 ~ philly$vetnorm + philly$incomefix + philly$servnorm, data = philly, gweight = gwr.Gauss)
philly.gwr <- gwr(philly$narcnorm2 ~ philly$vetnorm + philly$incomefix + philly$servnorm, data = philly, bandwidth = philly.bw, gweight = gwr.Gauss, hatmatrix = T)
gwr.df <- slot(philly.gwr$SDF, "data")

summary(philly.gwr)
names(philly.gwr)
names(philly.gwr$SDF)
print(philly.gwr)
anova(philly.gwr)

#calculate the t-score and p value
gwr.df$vet_tval <- gwr.df$philly.vetnorm / gwr.df$philly.vetnorm_se
gwr.df$vet_pval <-pt(gwr.df$vet_tval, 380, lower.tail = FALSE)
gwr.df$income_tval <- gwr.df$philly.incomefix / gwr.df$philly.incomefix_se
gwr.df$income_pval <- pt(gwr.df$income_tval, 380, lower.tail = FALSE)
gwr.df$serv_tval <- gwr.df$philly.servnorm / gwr.df$philly.servnorm_se
gwr.df$serv_pval <- pt(gwr.df$serv_tval, 380, lower.tail = FALSE)

pdf(file = "f8.vetparamchoropleth.pdf")
vetparamshades = auto.shading(gwr.df$philly.vetnorm, cols=brewer.pal(5,"Greens"))
choropleth(philly, gwr.df$philly.vetnorm, shading=vetparamshades)
title("Figure 8: GWR Parameter Estimate for Veteran Status")
choro.legend(2729559,251723.8, vetparamshades, title = "vet status")
dev.off()

pdf(file = "f9.vetpvalchoro.pdf")
vetpvalshades <- shading(c(0.005, 0.01, 0.05), cols = rev(brewer.pal(4, "Greens")))
choropleth(philly, gwr.df$vet_pval, shading = vetpvalshades)
title("Figure 9: P-Value for Veteran Status")
choro.legend(2729559,251723.8, vetpvalshades, title = "p-value")
dev.off()

pdf(file = "f10.incomeparamchoropleth.pdf")
par(mar=c(1,1,2,1))
incomeparamshades = auto.shading(gwr.df$philly.incomefix, cols=brewer.pal(5,"Greens"))
choropleth(philly, gwr.df$philly.incomefix, shading=incomeparamshades)
title("Figure 10: GWR Parameter Estimate for Per Capita Income")
choro.legend(2711398,241846.7, incomeparamshades, title = "per cap income")
dev.off()

pdf(file = "f11.incomepvalchoro.pdf")
incomepvalshades <- shading(c(0.005, 0.01, 0.05), cols = rev(brewer.pal(4, "Greens")))
choropleth(philly, gwr.df$income_pval, shading = incomepvalshades)
title("Figure 11: P-Value for Per Capita Income")
choro.legend(2729559,251723.8, incomepvalshades, title = "p-value")
dev.off()

pdf(file = "f12.servparamchoropleth.pdf")
servparamshades = auto.shading(gwr.df$philly.servnorm, cols=brewer.pal(5,"Greens"))
choropleth(philly, gwr.df$philly.servnorm, shading=servparamshades)
title("Figure 12: GWR Parameter Estimate for Service Jobs")
choro.legend(2729559,251723.8, servparamshades, title = "service jobs")
dev.off()

pdf(file = "f13.servpvalchoro.pdf")
servpvalshades <- shading(c(0.005, 0.01, 0.05), cols = rev(brewer.pal(4, "Greens")))
choropleth(philly, gwr.df$serv_pval, shading=servpvalshades)
title("Figure 13: P-Value for Service Jobs")
choro.legend(2729559,251723.8, servpvalshades, title = "p-value")
dev.off()

#R^2
pdf(file = "f14.localr2.pdf")
localr2shades = auto.shading(gwr.df$localR2, cols=brewer.pal(5,"Greens"))
choropleth(philly, gwr.df$localR2, shading=localr2shades)
title("Figure 14: GWR Local R² values")
choro.legend(2729559,251723.8, localr2shades, title = "Local R²")
dev.off()

#Spatial Autoregressive Models
lm04.lag <- lagsarlm(philly$narcnorm2 ~ philly$vetnorm + philly$incomefix + philly$servnorm, data = philly, philly_neigh_lw)
summary(lm04.lag)
anova(lm04.lag, lm04)
