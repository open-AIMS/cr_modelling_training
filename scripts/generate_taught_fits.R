# ---------------------------------------------------------------------------
# Fits the models the taught modules use, and saves them, so that a participant
# working through the material does not wait for the sampler.
#
# This is the only thing that fits them. The modules show each `bnec()` call and
# the `save()` that keeps it, marked `eval: false`, and execute a `load()`
# instead, which is the workflow an analysis follows. A render therefore samples
# nothing, and the objects have to come from somewhere: here.
#
# Why it matters. Measured on 2026-09-16 on a four-core WSL2 machine, rendering
# the development profile with the modules fitting took 35 minutes; with them
# loading it took 4 minutes 51 seconds and sampled nothing at all. A participant
# running module 5 alongside the presenter would otherwise wait 13 minutes, and
# module 6 nearly 16.
#
# Run with:
#   Rscript scripts/generate_taught_fits.R
#
# Expect roughly half an hour. Run it once, bundle the result with
# scripts/bundle_fits.R, and distribute that.
# ---------------------------------------------------------------------------

suppressMessages({
  library(bayesnec)
  library(brms)
  library(dplyr)
})

options(brms.backend = "cmdstanr", mc.cores = 4)

FORCE <- "--force" %in% commandArgs(trailingOnly = TRUE)

vig_dir <- "vignettes"
out_dir <- file.path(vig_dir, "fits")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

MODULES <- c(
  "2Fitting-a-CR-model-using-bayesnec",
  "4Model_averaging_and_multimodel_inference",
  "5Response_data_and_statistical_distributions",
  "6Priors_and_Bayesian_inference"
)

# Which chunks to run.
#
# An executed chunk is run, because later chunks depend on the data it reads and
# the objects it builds. A chunk marked `eval: false` is skipped, except where it
# saves into `fits/`: those are the fit chunks, and they are the whole point of
# this script. Skipping every `eval: false` chunk, which an earlier version did,
# left it unable to fit anything at all.
#
# knitr::purl() cannot make this distinction. Chunk option hooks are applied
# when a document is knitted and not when it is purled, so purling extracts the
# `eval: false` demonstrations as well -- module 4's parallel plan among them --
# and sourcing the result would run them.
runnable_chunks <- function(qmd) {
  lines <- readLines(qmd, warn = FALSE)
  starts <- grep("^```\\{r[ ,}]", lines)
  ends <- grep("^```\\s*$", lines)
  out <- character()
  seen <- character()
  for (s in starts) {
    e <- ends[ends > s][1]
    if (is.na(e)) next
    body <- lines[(s + 1):(e - 1)]
    opts <- grep("^#\\|", body, value = TRUE)
    code <- body[!grepl("^#\\|", body)]
    eval_false <- any(grepl("^#\\|\\s*eval:\\s*false", opts))
    targets <- sub('.*file\\s*=\\s*"([^"]+)".*', "\\1",
                   grep('save\\s*\\(.*file\\s*=\\s*"fits/', code, value = TRUE))
    if (eval_false && length(targets) == 0) next
    # A module can show the same fit call twice, once to fit and once to
    # illustrate the fit-save-load pattern, which would fit the object twice.
    # A chunk whose every target has already been written this run is skipped.
    if (length(targets) > 0 && all(targets %in% seen)) next
    seen <- c(seen, targets)
    out <- c(out, code, "")
  }
  out
}

# Resets the environment attached to each stored `bayesnecformula`.
#
# The module's own `save()` call writes the object exactly as an analysis would,
# and a formula keeps a reference to the environment it was written in. Every
# chunk of a module is evaluated here in one shared environment, so that
# reference reaches every object the module has built so far, and `save()`
# follows it. Measured on 2026-09-17: `m5_exp_5` serialised to 357 MB, of which
# 354.58 MB was the formula environment holding the seven fits made before it,
# against 2.8 MB for the same object with the environment reset. The whole set
# came to 534.5 MB rather than 24 MB.
#
# The environment is reset after the fact rather than avoided, because avoiding
# it means evaluating each fit in isolation and the later chunks of a module
# depend on objects the earlier ones made. The new parent is globalenv() rather
# than baseenv() so that `crf()` and the family functions still resolve if a
# method re-evaluates the formula. Estimates are unchanged by this: `nec()`,
# `ecx()` and `nsec()` on the stripped `m5_exp_5` returned the same values to
# five decimal places, and `check_fit()` and `plot()` both ran.
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

# Rewrites a saved file with its formula environments reset.
clean_saved <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  before <- file.size(path)
  e <- new.env()
  nms <- load(path, envir = e)
  for (nm in nms) assign(nm, strip_formula_env(get(nm, envir = e)), envir = e)
  save(list = nms, file = path, envir = e)
  message(sprintf("    %-26s %7.1f -> %6.1f MB", basename(path),
                  before / 1024^2, file.size(path) / 1024^2))
  invisible(TRUE)
}

message("fitting: ", format(Sys.time()))
t0 <- Sys.time()

for (mod in MODULES) {
  qmd <- file.path(vig_dir, paste0(mod, ".qmd"))
  stopifnot(file.exists(qmd))
  message("\n--- ", mod, " ---")
  m0 <- Sys.time()

  code <- runnable_chunks(qmd)
  targets <- sub('.*file\\s*=\\s*"([^"]+)".*', "\\1",
                 grep('save\\s*\\(.*file\\s*=\\s*"fits/', code, value = TRUE))
  message("  ", length(targets), " fit(s) to make")

  # A module whose fits are all on disk is skipped, so that correcting one
  # module costs that module rather than the whole set. --force redoes
  # everything.
  if (!FORCE && length(targets) > 0 &&
      all(file.exists(file.path(vig_dir, targets)))) {
    message("  skipping, every fit already saved")
    next
  }

  env <- new.env(parent = globalenv())
  owd <- setwd(vig_dir)          # chunks resolve data by relative path
  grDevices::pdf(NULL)
  ok <- tryCatch({
    eval(parse(text = code), envir = env)
    TRUE
  }, error = function(e) {
    message("  FAILED: ", conditionMessage(e))
    FALSE
  })
  grDevices::dev.off()
  setwd(owd)

  if (ok) for (tg in targets) clean_saved(file.path(vig_dir, tg))

  message("  ", if (ok) "done in " else "stopped after ",
          format(round(difftime(Sys.time(), m0), 1)))
}

saved <- list.files(out_dir, pattern = "\\.RData$", full.names = TRUE)
writeLines(c(
  "Generated by scripts/generate_taught_fits.R",
  paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("elapsed:", format(round(difftime(Sys.time(), t0), 1))),
  paste("R:", getRversion()),
  paste("bayesnec:", as.character(packageVersion("bayesnec"))),
  paste("brms:", as.character(packageVersion("brms"))),
  paste("cmdstanr:", as.character(packageVersion("cmdstanr"))),
  paste("objects:", length(saved)),
  paste("MB:", round(sum(file.size(saved)) / 1024^2, 1))
), file.path(out_dir, "PROVENANCE.txt"))

message("\ndone: ", format(Sys.time()))
message(length(saved), " objects in ", out_dir,
        " (", round(sum(file.size(saved)) / 1024^2, 1), " MB)")
