## Proof of whisky stored in a charred oak barrel against time in years.
##
## Schoeneman, R.; Dyer, R.; Earl, E. (1971), “Analytical profile of straight 
## bourbon whiskies”, Journal of the Association of Official Analytical 
## Chemists, vol. 54, pp. 1247-1261.
whiskey <- data.frame(
  age = c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8),
  proof = c(104.6, 104.1, 104.4, 105.0, 106.0, 
            106.8, 107.7, 108.7, 110.6, 112.1)
)

## Scatterplot of whiskey data
plot(proof ~ age, data = whiskey, pch = 1)
#with(whiskey, lines(lowess(age, proof), col = "dodgerblue2", lwd = 2))
whiskey.css <- with(whiskey, smooth.spline(age, proof)) # df = 6.281831
lines(whiskey.css, col = "dodgerblue2", lwd = 2)

## Prediction interval based on smooth 


## Load required packages
library(investr)

## Scatterplot with quadratic fit
whiskey.lm <- lm(proof ~ age + I(age^2), data = whiskey)
plotFit(whiskey.lm, pch = 19, col = "dodgerblue2", lwd.fit = 2, cex = 1.3, 
        interval = "prediction", lty.pred = 1)

## Calibration with proof = 110
abline(h = 110)
