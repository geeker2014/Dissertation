## Whiskey data from JSM 2000
whiskey <- data.frame(
  age = c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8),
  proof = c(104.6, 104.1, 104.4, 105.0, 106.0, 
            106.8, 107.7, 108.7, 110.6, 112.1)
)

## Load packages
library(investr)

## Fit models
mod1 <- lm(proof ~ age + I(age^2), data = whiskey)
mod2 <- with(whiskey, lmmSmooth(age, proof, degree = 2))
mod3 <- with(whiskey, lmmSmooth2(age, proof, degree = 2))

## Calibration intervals
invest(mod1, y0 = 108)                # 5.2329 (4.6776, 5.7352)
invest(mod2, y0 = 108)                # 5.2329 (4.6776, 5.7352)
invest(mod3, y0 = 108)                # 5.2329 (4.6776, 5.7352)

invest(mod2, y0 = 108, corrected = F) # 5.2329 (4.7223, 5.7016)
invest(mod3, y0 = 108, corrected = F) # 5.2329 (4.7223, 5.7016)