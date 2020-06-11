#library("lme4") #for analysis of posterior draws of a Bayesian model
library(tidyverse) # for data manipulation and plots
library(haven) #for reading sav data
library(sjstats) #for calculating intra-class correlation (ICC)
library(ROCR) #for calculating area under the curve (AUC) statistics
library(brms) #for Bayesian (multilevel) generalised linear modelling
library(modelr) #for data manipulation
library(tidybayes) #for analysis of posterior draws of a Bayesian model
library(pROC)

d <- read.csv("unlogged.csv")

simplesm <- lm(NLR ~ 1 + time_weeks + log(Staphylococcus+2e-6), data=d)
summary(simplesm)

m1 <- lme(log(NLR)~log(Staphylococcus+2e-6),random=~1|cT,data=d)
anova(m1)



tmp <- d%>%mutate(CTT = c(T),logstaph = log(Staphylococcus+2e-6),logburkholderiaceae = log(Burkholderiaceae+2e-6))

fit <- lmer(NLR ~ logstaph + (1|CTT), data=tmp)
summary(fit)

fit <- lmer(NLR ~ logstaph  + (1|CTT), data=tmp)
summary(fit)

mbrms <- brm(NLR~1+logstaph+time_weeks+(1|subject),data=tmp)
summary(mbrms)

mcmc_plot(mbrms,type = "areas", prob = 0.9,)

library(sjPlot)
library(sjlabelled)
library(sjmisc)
plot_model(mbrms, sort.est = TRUE,type = "est")
