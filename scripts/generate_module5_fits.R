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
#   Poisson, negbin    alga, Cladocopium proliferum against contaminant B.
#                      Cell counts whose control standard deviation is 1417
#                      against the 169 a Poisson implies.
#   gamma              lum31, one zinc plate read at 15 minutes.
#                      Luminescence, positive throughout, whose coefficient
#                      of variation runs 0.02 to 0.21 across the series.
#   Gaussian           alga again, the specific growth rate, which is genuinely
#                      negative at the two highest doses.
#
# The equation fitted to each is the one scripts/screen_module5_models.R found
# best supported among those passing check_sampling(), with its evidence in
# notes/module5-model-screen/ (2026-09-23). The beta keeps nec3param beside
# ecxwb1 so the module can show a threshold equation losing the weight on a
# series with no flat region.
#
# The dispersion example follows the gamma on the same plate, so the module
# fails a check and repairs it on one dataset rather than three. The beta is
# repaired the same way, in its own section.
#
# The binomial pair and the hurdle example are here too, so that this script
# owns every fit module 5 loads. They were made by generate_taught_fits.R until
# 2026-09-22, which stopped being true when module 5 was removed from its
# MODULES list, leaving four objects with no generator at all.
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

# The data preparation is shared with scripts/screen_module5_models.R.
source("scripts/module5_data.R")

message("fitting: ", format(Sys.time()))
t0 <- Sys.time()

if (needed("exp_1nec")) {
  message("  binomial: threshold equation")
  exp_1nec <- bnec(suc | trials(tot) ~ crf(log.x, model = "nec3param"),
                   data = binom_data, seed = 333)
  save_fit(exp_1nec, "exp_1nec")
}

if (needed("exp_1ecx")) {
  message("  binomial: smooth equation")
  exp_1ecx <- bnec(suc | trials(tot) ~ crf(log.x, model = "ecxll3"),
                   data = binom_data, seed = 333)
  save_fit(exp_1ecx, "exp_1ecx")
}

if (needed("exp_1b")) {
  message("  beta-binomial: the same counts")
  exp_1b <- bnec(suc | trials(tot) ~ crf(log.x, model = c("ecxll3", "nec3param")),
                 data = binom_data, family = "beta_binomial", seed = 333)
  save_fit(exp_1b, "exp_1b")
}

if (needed("exp_2")) {
  message("  beta: coral_pam yield")
  exp_2 <- bnec(yield ~ crf(log.x, model = c("ecxwb1", "nec3param")),
                data = pam_data, seed = 333)
  save_fit(exp_2, "exp_2")
}

if (needed("exp_2_disp")) {
  message("  beta with disp(\"power\"): coral_pam yield")
  exp_2_disp <- bnec(yield ~ crf(log.x, model = "ecxwb1") + disp("power"),
                     data = pam_data, seed = 333)
  save_fit(exp_2_disp, "exp_2_disp")
}

if (needed("exp_3")) {
  message("  poisson: alga cell density")
  exp_3 <- bnec(density ~ crf(log.x, model = "ecxwb1"),
                data = alga_data, family = "poisson", seed = 333)
  save_fit(exp_3, "exp_3")
}

if (needed("exp_3b")) {
  message("  negative binomial: alga cell density")
  exp_3b <- bnec(density ~ crf(log.x, model = "ecxwb1"),
                 data = alga_data, family = "negbinomial", seed = 333)
  save_fit(exp_3b, "exp_3b")
}

if (needed("exp_4")) {
  message("  gamma: lum31 one zinc plate")
  exp_4 <- bnec(rlu ~ crf(log.x, model = "ecxll3"),
                data = lum_data, family = "Gamma", seed = 333)
  save_fit(exp_4, "exp_4")
}

if (needed("exp_5")) {
  message("  gaussian: alga specific growth rate")
  # adapt_delta is raised because ecxwb1 returns divergent transitions on these
  # data at the default of 0.8 and none at 0.99 (the screen, gaussian_b_ad).
  exp_5 <- bnec(sgr ~ crf(log.x, model = "ecxwb1"),
                data = alga_data, seed = 333,
                control = list(adapt_delta = 0.99))
  save_fit(exp_5, "exp_5")
}

# The dispersion sub-model, on the same plate the gamma section fits, so the
# module shows a check failing and being repaired on one dataset. exp_4 is the
# constant-dispersion baseline: the same equation on the same data.
if (needed("exp_4_disp")) {
  message("  gamma with disp(\"power\"): lum31 one zinc plate")
  exp_4_disp <- bnec(rlu ~ crf(log.x, model = "ecxll3") + disp("power"),
                     data = lum_data, family = "Gamma", seed = 333)
  save_fit(exp_4_disp, "exp_4_disp")
}

# adapt_delta is raised because the survival curve is determined by the top two
# doses and samples awkwardly at the default.
if (needed("exp_6")) {
  message("  hurdle gamma: nassarius contaminant A")
  exp_6 <- bnec_hurdle(growth ~ crf(log_dose, model = c("nec3param", "ecxll3")),
                       data = snail, seed = 333,
                       control = list(adapt_delta = 0.99))
  save_fit(exp_6, "exp_6")
}

message("done in ", format(round(difftime(Sys.time(), t0), 1)))
