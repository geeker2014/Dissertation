## Load packages and data ------------------------------------------------------
library(lattice)
library(nlme)
dogs <- read.csv("/home/w108bmg/Desktop/dogs", header = T)
# head(dogs)

## Plot data -------------------------------------------------------------------
xyplot(y ~ x, data = dogs, groups = dog, pch = 19, type = c("p", "r"))
xyplot(y ~ x, data = dogs, groups = condition, pch = 19, type = c("p", "r"))

xyplot(y ~ x | as.factor(dog), data = dogs, pch = 19, type = c("p", "r"))
xyplot(y ~ x | as.factor(condition), data = dogs, pch = 19, type = c("p", "r"))

## Interaction plot
with(dogs, interaction.plot(dog, condition, y))

## Dog looks like a random effect
dogs.lmList <- lmList(y ~ I(x-mean(x)) | dog, data = dogs)
plot(intervals(dogs.lmList))

## Condition looks like a fixed effect
dogs.lmList <- lmList(y ~ I(x-mean(x)) | condition, data = dogs)
plot(intervals(dogs.lmList))

## Fit models ------------------------------------------------------------------
mod1 <- lm(y ~ x + as.factor(condition), data = dogs)
mod2 <- lme(y ~ x + as.factor(condition), data = dogs, random = ~ 1|dog)
AIC(mod1, mod2)

library(lme4)
mod.lmer


