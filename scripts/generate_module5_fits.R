# ---------------------------------------------------------------------------
# Fits the module 5 examples that use datasets shipping with `bayesnec`.
#
# Module 5 teaches matching a statistical family to a response. Until
# 2026-09-22 four of its worked examples were constructed from `nec_data`:
# the Poisson and negative binomial counts were `abs(round(y * 100 + rnorm(0,
# 15)))`, the gamma was `exp(y)`, and the Gaussian was a logit of the
# proportion data. All four fitted their own families badly, and the reason was
# the construction rather than anything a participant would meet. RF raised it
# after delivering the course.
#
# The replacements are published datasets in the package:
#
#   beta               coral_pam, the 2018 climate scenario. PAM yield, a
#                      genuine physical proportion, with readings at exactly 0
#                      at the top of the series.
#   Poisson, negbin    alga, Cladocopium proliferum against contaminant A.
#                      Cell counts whose control standard deviation is 1454
#                      against the 136 a Poisson implies.
#   gamma              lum31, one zinc plate. Luminescence, positive
#                      throughout, whose coefficient of variation runs 0.10 to
#                      1.05 across the series.
#   Gaussian           alga again, the specific growth rate, which is genuinely
#                      negative at the two highest doses.
#
# The dispersion example follows the gamma on the same plate, so the module
# fails a check and repairs it on one dataset rather than three.
#
# Run with:
#   Rscript scripts/generate_module5_fits.R [--force]
# ---------------------------------------------------------------------------

suppressMessages({
  library(bayesnec)
  library(dplyr)
})

options(brms.backend = "cmdstanr", mc.cores = 4)
FORCE <- "--force" %in% commandArgs(trailingOnly = TRUE)
out_dir <- "vignettes/fits"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# The formula environment is reset before saving, for the reason recorded in
# scripts/generate_taught_fits.R: a saved formula keeps a reference to the
# environment it was written in, and here that environment holds every other
# fit made by this script.
strip_formula_env <- function(x) {
  clean <- function(f) { environment(f) <- new.env(parent = globalenv()); f }
  if (inherits(x, "bayesmanecfit")) {
    x$mod_fits <- lapply(x$mod_fits, function(m) {
      m$bayesnecformula <- clean(m$bayesnecformula); m
    })
  } else if (!is.null(x$bayesnecformula)) {
    x$bayesnecformula <- clean(x$bayesnecformula)
  }
  x
}

save_fit <- function(obj, name) {
  assign(name, strip_formula_env(obj))
  path <- file.path(out_dir, paste0("m5_", name, ".RData"))
  save(list = name, file = path)
  message(sprintf("  wrote %-22s %5.1f MB", basename(path),
                  file.size(path) / 1024^2))
}

needed <- function(name) {
  FORCE || !file.exists(file.path(out_dir, paste0("m5_", name, ".RData")))
}

data(coral_pam); data(alga); data(lum31)

# --- the data each example uses, prepared exactly as the module shows --------

# A log predictor needs the control somewhere other than log(0). One decade
# below the lowest tested level is the convention the module uses throughout.
pam_data <- subset(coral_pam, climate == "2018")
pam_data$log.x <- log(pam_data$diuron + min(pam_data$diuron[pam_data$diuron > 0]) / 10)

alga_data <- subset(alga, species == "c_proliferum" & contaminant == "A")
alga_data$log.x <- log(alga_data$dose + min(alga_data$dose[alga_data$dose > 0]) / 10)

# One plate rather than the whole series, because module 8 is where the
# plate-to-plate structure is modelled. This plate is chosen because every
# reading on it is positive, which a gamma requires; the copper plates include
# readings floored at zero.
lum_data <- subset(lum31, plate == "Zn 28Mar24 Rep1B")
lum_data$log.x <- log(lum_data$conc)

message("fitting: ", format(Sys.time()))
t0 <- Sys.time()

if (needed("exp_2")) {
  message("  beta: coral_pam yield")
  exp_2 <- bnec(yield ~ crf(log.x, model = c("ecxll3", "nec3param")),
                data = pam_data, seed = 333)
  save_fit(exp_2, "exp_2")
}

if (needed("exp_3")) {
  message("  poisson: alga cell density")
  exp_3 <- bnec(density ~ crf(log.x, model = c("ecxll3", "nec3param")),
                data = alga_data, family = "poisson", seed = 333)
  save_fit(exp_3, "exp_3")
}

if (needed("exp_3b")) {
  message("  negative binomial: alga cell density")
  exp_3b <- bnec(density ~ crf(log.x, model = c("ecxll3", "nec3param")),
                 data = alga_data, family = "negbinomial", seed = 333)
  save_fit(exp_3b, "exp_3b")
}

if (needed("exp_4")) {
  message("  gamma: lum31 one zinc plate")
  exp_4 <- bnec(rlu ~ crf(log.x, model = c("ecxll3", "nec3param")),
                data = lum_data, family = "Gamma", seed = 333)
  save_fit(exp_4, "exp_4")
}

if (needed("exp_5")) {
  message("  gaussian: alga specific growth rate")
  exp_5 <- bnec(sgr ~ crf(log.x, model = c("ecxll4", "nec4param")),
                data = alga_data, seed = 333)
  save_fit(exp_5, "exp_5")
}

# The dispersion sub-model, on the same plate the gamma section fits, so the
# module shows a check failing and being repaired on one dataset. The
# constant-dispersion baseline needs no separate fit: it is pulled out of
# exp_4 above, which makes the comparison like for like.
if (needed("exp_4_disp")) {
  message("  gamma with disp(\"power\"): lum31 one zinc plate")
  exp_4_disp <- bnec(rlu ~ crf(log.x, model = "nec3param") + disp("power"),
                     data = lum_data, family = "Gamma", seed = 333)
  save_fit(exp_4_disp, "exp_4_disp")
}

message("done in ", format(round(difftime(Sys.time(), t0), 1)))
