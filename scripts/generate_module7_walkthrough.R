# ---------------------------------------------------------------------------
# Fits the two ametryn model sets that module 7 works through step by step, and
# saves them for the module to load.
#
# Module 7 follows the workflow of `bayesnec`'s example9 vignette: prepare the
# data and choose the family, fit the candidate set, screen the sampler, check
# the fit, report. Two fits are needed rather than one, because step 3 is a
# demonstration and not an assertion:
#
#   m7_am_reduced   iter = 4000 over 2 chains, below bnec()'s defaults
#   m7_am_fit       bnec()'s defaults, with adapt_delta and max_treedepth set
#
# These data are well behaved at the defaults, so a fit made there produces no
# screen failure to examine. Running short guarantees that some equations miss a
# threshold, which is what lets the module show the three kinds of failure and
# what each of them calls for. The reduced fit is a teaching device and every
# estimate module 7 reports comes from the fit made at the defaults.
#
# The objects are saved rather than fitted at render time because the decline
# set is fourteen equations and takes about a quarter of an hour. The module
# shows the bnec() call and the save(), marked `eval: false`, and runs the
# load(), which is the workflow modules 2 to 6 follow.
#
# Run with:
#   Rscript scripts/generate_module7_walkthrough.R [--force]
# ---------------------------------------------------------------------------

suppressMessages({
  library(bayesnec)
  library(brms)
})

options(brms.backend = "cmdstanr", mc.cores = 4)

FORCE <- "--force" %in% commandArgs(trailingOnly = TRUE)
out_dir <- "vignettes/fits"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SEED <- 17

# The equation pull_best() returns from the screened set; see the dispersion
# section below for why it is named here rather than derived.
BEST_MODEL <- "ecx4param"

# The formula environment is reset before saving, for the reason recorded in
# scripts/generate_taught_fits.R: a saved formula keeps a reference to the
# environment it was written in, and here that environment holds the other fit.
strip_formula_env <- function(x) {
  clean <- function(f) { environment(f) <- new.env(parent = globalenv()); f }
  if (inherits(x, "bayesmanecfit")) {
    x$mod_fits <- lapply(x$mod_fits, function(m) {
      m$bayesnecformula <- clean(m$bayesnecformula)
      m
    })
  } else if (!is.null(x$bayesnecformula)) {
    x$bayesnecformula <- clean(x$bayesnecformula)
  }
  x
}

save_fit <- function(obj, name) {
  obj <- strip_formula_env(obj)
  assign(name, obj)
  path <- file.path(out_dir, paste0(name, ".RData"))
  save(list = name, file = path)
  message(sprintf("  wrote %s (%.1f MB)", basename(path),
                  file.size(path) / 1024^2))
}

message("fitting: ", format(Sys.time()))
t0 <- Sys.time()

data(herbicide)
ametryn <- herbicide[herbicide$herbicide == "ametryn", ]

targets <- file.path(out_dir, c("m7_am_reduced.RData", "m7_am_fit.RData"))
if (!FORCE && all(file.exists(targets))) {
  message("  skipping, both fits already saved")
} else {
  # Below the defaults on purpose. bnec() takes chains = 4 and iter = 1e4, and
  # sets warmup to floor(iter / 5) * 4, so it retains 8000 draws; this call
  # retains 1600.
  message("  ametryn: reduced settings")
  m7_am_reduced <- bnec(fvfm ~ crf(log(concentration), model = "decline"),
                        data = ametryn, family = Beta(link = "identity"),
                        iter = 4000, chains = 2, seed = SEED)
  save_fit(m7_am_reduced, "m7_am_reduced")

  # bnec()'s defaults for iter and chains, with both control arguments set.
  # Neither is needed on these data -- no equation records a divergent
  # transition at either settings -- and the module fits with them anyway, so
  # that the worked example shows the two arguments a marginal fit needs rather
  # than only naming them. Step 3 of the module says which failure each one
  # answers.
  message("  ametryn: bnec() defaults, adapt_delta = 0.99, max_treedepth = 12")
  m7_am_fit <- bnec(fvfm ~ crf(log(concentration), model = "decline"),
                    data = ametryn, family = Beta(link = "identity"),
                    seed = SEED,
                    control = list(adapt_delta = 0.99, max_treedepth = 12))
  save_fit(m7_am_fit, "m7_am_fit")
}

# --- the dispersion sub-model ------------------------------------------------
# One equation refitted with the dispersion free to vary along the predictor, so
# module 7 can show what repairing the control does rather than only naming the
# remedy. The equation has to be the one pull_best() returns from the screened
# set, so that this is a like-for-like comparison against it and the baseline
# needs no separate fit.
#
# That is ecx4param as of 2026-09-21, at a weight of 0.405 against 0.390 for
# ecxll4, which it was until the default fit gained adapt_delta = 0.99. The two
# are close enough that the selection can change with the settings or the seed,
# so read pull_best() rather than assuming, and change BEST_MODEL below and the
# name in the module together if it moves again.
#
# disp(~ log(concentration)) is used rather than disp("power") because it fits
# these data better. Measured 2026-09-17, the control sd_ratio is 0.924 against
# 0.717 and one group of eight falls outside the posterior predictive interval
# against two.
disp_target <- file.path(out_dir, "m7_am_disp.RData")
if (FORCE || !file.exists(disp_target)) {
  message("  ametryn: ", BEST_MODEL, " with a dispersion sub-model")
  m7_am_disp <- bnec(fvfm ~ crf(log(concentration), model = BEST_MODEL) +
                       disp(~ log(concentration)),
                     data = ametryn, family = Beta(link = "identity"),
                     seed = SEED)
  save_fit(m7_am_disp, "m7_am_disp")
} else {
  message("  skipping, dispersion fit already saved")
}

# --- provenance -------------------------------------------------------------
writeLines(c(
  "Generated by scripts/generate_module7_walkthrough.R",
  paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("elapsed:", format(round(difftime(Sys.time(), t0), 1))),
  paste("seed:", SEED),
  paste("series: ametryn,", nrow(ametryn), "observations,",
        length(unique(ametryn$concentration)), "concentrations"),
  paste("R:", getRversion()),
  paste("bayesnec:", as.character(packageVersion("bayesnec"))),
  paste("brms:", as.character(packageVersion("brms"))),
  "reduced fit: iter = 4000, chains = 2 (1600 retained draws)",
  "default fit: bnec() defaults (8000 retained draws),",
  "              control = list(adapt_delta = 0.99, max_treedepth = 12)",
  paste0("dispersion fit: ", BEST_MODEL,
         " + disp(~ log(concentration)), bnec() defaults")
), "vignettes/data/module7_provenance.txt")

message("done: ", format(Sys.time()))
