# For reference, see
#http://www.ats.ucla.edu/stat/r/dae/logit.htm
#setwd("C:/Prostate_LogisticRegression/Analysis")
setwd("C:/Users/Owner/Documents/GitHub/MSDS_6372/Prostate_LogisticRegression/Analysis")

#install.packages("aod")
library(aod)

require(reshape2)

#install.packages("ggplot2")
library(ggplot2)
library(Rcpp)

#install.packages("car")
library(car)

#prostatedata <- read.csv("C:/Users/Johnny/Documents/6372/Project3/pros.csv")
prostatedata <- read.csv("pros.csv")

## view the first few rows of the data
head(prostatedata)
summary(prostatedata)
nrow(prostatedata) ##number of rows is 380

## NAs are present in VOL (1) and RACE (3), so we need to remove them
prostatedata.cleaned <- na.omit(prostatedata)
summary(prostatedata.cleaned) ## no NAs are left
nrow(prostatedata.cleaned) ##number of rows is 376

## histograms
# Observe attribute distributions
prostate.subset <- subset(prostatedata.cleaned, select = c(AGE, PSA, VOL, GLEASON))
prostate.melt <- melt(prostate.subset[sapply(prostate.subset, is.numeric)])
ggplot(data = prostate.melt, mapping = aes(x = value)) + 
  geom_histogram(bins = 10) + facet_wrap(~variable, scales = 'free_x')

## PSA is right skewed, so log transform it (VOL is but contains many 0s)
prostatedata.cleaned$logPSA <- log(prostatedata.cleaned$PSA)

# Observe attribute distributions after transform
prostate.subset <- subset(prostatedata.cleaned, select = c(AGE, logPSA, VOL, GLEASON))
prostate.melt <- melt(prostate.subset[sapply(prostate.subset, is.numeric)])
ggplot(data = prostate.melt, mapping = aes(x = value)) + 
  geom_histogram(bins = 10) + facet_wrap(~variable, scales = 'free_x')

## Mean, Median, and Standard Deviation
sapply(prostatedata.cleaned, mean)
sapply(prostatedata.cleaned, median)
sapply(prostatedata.cleaned, sd)

#correlation
# Observe correlations between variables
write.csv(cor(prostatedata.cleaned[sapply(prostatedata.cleaned, is.numeric)]), file = "Prostate_Correlations.csv")

# Set RACE, DPROS, and DCAPS to factors
prostatedata.cleaned$RACE <- factor(prostatedata.cleaned$RACE)
prostatedata.cleaned$DPROS <- factor(prostatedata.cleaned$DPROS)
prostatedata.cleaned$DCAPS <- factor(prostatedata.cleaned$DCAPS)

attach(prostatedata.cleaned)

# First Logit with all variables included
prostatelogit <- glm(CAPSULE ~ AGE + RACE + DPROS + DCAPS + logPSA + VOL + GLEASON, data = prostatedata.cleaned, family = "binomial")

summary(prostatelogit)

prostatelogit.interaction <- glm(CAPSULE ~  AGE*RACE + AGE + RACE + DPROS*DCAPS + DPROS + DCAPS + logPSA + VOL + GLEASON, data = prostatedata.cleaned, family = "binomial")
summary(prostatelogit.interaction)
## no interactions appear to exist between DCAPS and DPROS or AGE and RACE, so continue with no interactions in model

## Remove AGE, RACE, DCAPS, AND VOL as they are not statistically significant
prostatelogit2 <- glm(CAPSULE ~ DPROS + logPSA + GLEASON, data = prostatedata.cleaned, family = "binomial")
summary(prostatelogit2)

detach(prostatedata.cleaned)

##goodness of fit
with(prostatelogit, null.deviance - deviance)
with(prostatelogit, df.null - df.residual)
with(prostatelogit, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
logLik(prostatelogit)
AIC(prostatelogit)
BIC(prostatelogit)

with(prostatelogit.interaction, null.deviance - deviance)
with(prostatelogit.interaction, df.null - df.residual)
with(prostatelogit.interaction, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
logLik(prostatelogit.interaction)
AIC(prostatelogit.interaction)
BIC(prostatelogit.interaction)

with(prostatelogit2, null.deviance - deviance)
with(prostatelogit2, df.null - df.residual)
with(prostatelogit2, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
logLik(prostatelogit2)
AIC(prostatelogit2)
BIC(prostatelogit2)

## log odds and odds ratios with 95% CI
logOdds <- cbind(Log.Odds = coef(prostatelogit2), confint(prostatelogit2))
oddsRatio <- exp(cbind(OR = coef(prostatelogit2), confint(prostatelogit2)))
cbind(logOdds, oddsRatio)

## predictions
predprostatedata1 <- with(prostatedata.cleaned,
                 data.frame(DPROS = factor(1:4), logPSA = mean(logPSA), GLEASON = mean(GLEASON)))
predprostatedata1

predprostatedata1$DPROSP <- predict(prostatelogit2, newdata = predprostatedata1, type = "response")
predprostatedata1

predprostatedata2 <- with(prostatedata.cleaned,
                 data.frame(logPSA = rep(seq(from = min(prostatedata.cleaned$logPSA), to = max(prostatedata.cleaned$logPSA), length.out = 10), 4),
                            GLEASON = mean(GLEASON), DPROS = factor(rep(1:4, each = 10))))
predprostatedata3 <- cbind(predprostatedata2, predict(prostatelogit2, newdata = predprostatedata2, type="link", se=TRUE))
predprostatedata3 <- within(predprostatedata3, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})
predprostatedata3

prostatedata.cleaned$PredictedProb <- prostatedata.cleaned$CAPSULE
ggplot(predprostatedata3, aes(x = logPSA, y = PredictedProb)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = DPROS), alpha = .2) +
  geom_line(aes(colour = DPROS), size=1) +
  geom_jitter(data = prostatedata.cleaned, height = 0.05, width = 0) +
  annotate("text", x = -1, y = 0.85, label = "Note: Jittered dots represent\nactual response", hjust = 0) +
  labs(title = "Predicted Probabilities and 95% CIs for Prostate Cancer")
