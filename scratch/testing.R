require(tidyverse)
require(cmdstanr)
set_cmdstan_path("C:/cmdstan")
options(brms.backend = "cmdstanr", mc.cores = 8)
library(bayesnec)

dat <- read.csv("vignettes/example_ogl.csv") |> 
  dplyr::mutate(y=Perc_change/100,
                groupvar=factor(paste(Replicate, Treatment, sep="_")),
                x=as.numeric(Treatment)) |> 
  na.omit()

dat2 <- read.csv("vignettes/example_ogl.csv") |> 
  dplyr::mutate(y=jitter(Perc_change/100),
                groupvar=factor(paste("rep2", Replicate, Treatment, sep="_")),
                x=as.numeric(Treatment)) |> 
  na.omit()


head(dat)

newdat <- rbind(dat, dat2)

table(newdat$groupvar)
# five individuals in each group

ggplot(dat, aes(x=Treatment, y=y, colour=groupvar)) +
  geom_jitter()

# let's see which model is good
fit_all <- bnec(y~crf(x, model = "decline"), data = dat)
range(dat$y)

fit <- bnec(y~crf(x, model = "ecxlin") + ogl(groupvar), data = dat)
fit.2 <- bnec(y~crf(x, model = "ecxlin"), data = dat)

fit.3 <- brm(y~x + (1|groupvar), data = dat)


fit.4 <- bnec(y~crf(x, model = "ecxlin") + (top | groupvar), data = dat)

priors <- pull_prior(fit.2)

y ~ ogl + top - exp(slope) * x 
top ~ 1
slope ~ 1
ogl ~ 1 + (1 | groupvar)

pred3 <- posterior_epred(fit.3)


fit.original <- pull_out(fit_all, model = "ecxlin")

autoplot(fit.original)

require(lme4)

tt <- lmer(y~x + (1|groupvar), data=dat)
tt2 <- lm(y~x, data = dat)
summary(tt)
summary(tt2)

cowplot::plot_grid(autoplot(fit.2), autoplot(fit.4), labels = c("no RE", "R Intercept"))
summary(fit.2$fit)$spec_pars
summary(fit.2$fit)$fixed


summary(fit.4$fit)$spec_pars
summary(fit.4$fit)$fixed
summary(fit.4$fit)$random

summary(fit$fit)$spec_pars
summary(fit$fit)$fixed
summary(fit$fit)$random
summary(fit)

f1 <- pull_brmsfit(fit.2)
f2 <- pull_brmsfit(fit.4)

loo_controls <- list(fitting = list(), weights = list(method = "pseudobma"))
loo_model_weights(f1, f2, method = "pseudobma")


