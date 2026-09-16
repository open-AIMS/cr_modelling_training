# ---------------------------------------------------------------------------
# Fits the seven herbicide concentration-response curves used by module 7 and
# writes the compact artefacts that module reads.
#
# Seven herbicides by the fourteen equations of the `decline` set is 98 model
# fits and takes hours. It is therefore run deliberately, not at render time,
# and only its outputs are committed. The fitted objects themselves are far too
# large for a public repository: the equivalent from 2023 was 451 MB.
#
# Run with:
#   Rscript scripts/generate_herbicide_fits.R
# ---------------------------------------------------------------------------

suppressMessages({
  library(bayesnec)
  library(dplyr)
  library(purrr)
  library(tidyr)
})

options(mc.cores = 4)
out_dir <- "vignettes/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SEED <- 17

message("fitting: ", format(Sys.time()))
t0 <- Sys.time()

data(herbicide)

# The predictor is logged before fitting, as in the published analysis. Every
# estimate returned is therefore on the log scale and is back-transformed with
# xform = exp where it is reported.
manecfit <- herbicide |>
  mutate(concentration = log(concentration)) |>
  split(f = ~ herbicide) |>
  map(function(x) {
    bnec(fvfm ~ crf(concentration, model = "decline"), data = x, seed = SEED)
  })

message("fitted in ", format(round(difftime(Sys.time(), t0), 1)))

# --- screen out any model that did not sample adequately ---------------------
# screen_models() applies all three sampler diagnostics -- Rhat, effective
# sample size and divergent transitions -- and drops what failed, reporting
# which models went and why. It replaces a local helper written over rhat()
# alone, which had two defects: it screened on one diagnostic of the three, and
# it shadowed bayesnec::failed_models(), which reports something else entirely
# (the equations that produced no draws at all).
#
# The full per-model table is written out, not only the failures, because a
# methods section has to state the thresholds applied as well as the result.
sampling <- imap_dfr(manecfit, function(fit, nm) {
  check_sampling(fit) |>
    mutate(herbicide = nm, .before = 1)
})
write.csv(sampling, file.path(out_dir, "herbicide_sampling.csv"),
          row.names = FALSE)

# Equations that produced no draws at all are a different failure from equations
# that sampled badly, and are reported separately.
never_fitted <- imap_dfr(manecfit, function(fit, nm) {
  nms <- names(failed_models(fit))
  if (length(nms) == 0) return(NULL)
  data.frame(herbicide = nm, model = nms)
})
write.csv(never_fitted, file.path(out_dir, "herbicide_never_fitted.csv"),
          row.names = FALSE)

# quiet = FALSE so the exclusions are named in the log this script writes.
cleaned_fits <- map(manecfit, function(fit) screen_models(fit))

# --- model weights ----------------------------------------------------------
imap_dfr(cleaned_fits, function(fit, nm) {
  fit$mod_stats |>
    select(model, waic, wi) |>
    mutate(herbicide = nm, .before = 1)
}) |>
  write.csv(file.path(out_dir, "herbicide_mod_stats.csv"), row.names = FALSE)

# --- fitted curves, for the panel figure ------------------------------------
# The predictor is returned on the fitted (log) scale and is exponentiated here,
# so the module plots concentrations directly.
imap_dfr(cleaned_fits, function(fit, nm) {
  fit$w_pred_vals$data |>
    mutate(herbicide = nm, concentration = exp(x), .before = 1) |>
    select(herbicide, concentration, Estimate, Q2.5, Q97.5)
}) |>
  write.csv(file.path(out_dir, "herbicide_pred_vals.csv"), row.names = FALSE)

# --- the observations, so the figure can show the data ----------------------
herbicide |>
  select(herbicide, concentration, fvfm) |>
  write.csv(file.path(out_dir, "herbicide_observations.csv"), row.names = FALSE)

# --- no-effect posteriors and their pairwise comparison ---------------------
# "n(s)ec" is the weighted mixture of NEC values from the threshold equations
# and NSEC values from the smooth ones, which is what these mixed sets hold and
# what the module reports. "nec" would ask for the threshold parameter alone.
post_comp <- compare_posterior(cleaned_fits, comparison = "n(s)ec")

post_comp$posterior_data |>
  mutate(concentration = exp(value)) |>
  select(model, value, concentration) |>
  write.csv(file.path(out_dir, "herbicide_ne_posterior.csv"), row.names = FALSE)

write.csv(post_comp$prob_diff,
          file.path(out_dir, "herbicide_prob_diff.csv"), row.names = FALSE)

# --- curve parameters -------------------------------------------------------
# summary() reports weights and the no-effect estimate but no parameter
# estimates. curve_params() returns them per equation, with the model weight
# beside each row. xform reaches nec and ec50, which are on the predictor axis;
# top, bot and beta are not and are left alone.
imap_dfr(cleaned_fits, function(fit, nm) {
  curve_params(fit, xform = exp) |>
    mutate(herbicide = nm, .before = 1)
}) |>
  write.csv(file.path(out_dir, "herbicide_curve_params.csv"), row.names = FALSE)

# --- summarised estimates ---------------------------------------------------
imap_dfr(cleaned_fits, function(fit, nm) {
  q <- quantile(fit$w_ne_posterior, probs = c(0.5, 0.025, 0.975))
  data.frame(herbicide = nm,
             log_estimate = q[[1]], log_lower = q[[2]], log_upper = q[[3]],
             estimate = exp(q[[1]]), lower = exp(q[[2]]), upper = exp(q[[3]]))
}) |>
  write.csv(file.path(out_dir, "herbicide_estimates.csv"), row.names = FALSE)

# --- ametryn, for the single-herbicide walkthrough --------------------------
# Module 7 walks through one herbicide before comparing all seven. The 2023
# version did not run that fit live either, so its console output is captured
# here and displayed verbatim rather than re-fitted at render time.
ametryn_fit <- cleaned_fits[["ametryn"]]

writeLines(capture.output(summary(ametryn_fit)),
           file.path(out_dir, "ametryn_summary.txt"))

# Per-model fitted curves, for the plot showing every equation in the set.
# pull_out() returns a bayesnecfit for a single model, which carries pred_vals.
map_dfr(ametryn_fit$success_models, function(nm) {
  pull_out(ametryn_fit, model = nm)$pred_vals$data |>
    mutate(model = nm, concentration = exp(x), .before = 1) |>
    select(model, concentration, Estimate, Q2.5, Q97.5)
}) |>
  write.csv(file.path(out_dir, "ametryn_allmodels.csv"), row.names = FALSE)

# --- provenance -------------------------------------------------------------
writeLines(c(
  "Generated by scripts/generate_herbicide_fits.R",
  paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("elapsed:", format(round(difftime(Sys.time(), t0), 1))),
  paste("seed:", SEED),
  paste("R:", getRversion()),
  paste("bayesnec:", as.character(packageVersion("bayesnec"))),
  paste("brms:", as.character(packageVersion("brms"))),
  paste("cmdstanr:", as.character(packageVersion("cmdstanr"))),
  paste("model set: decline,", length(models()$decline), "equations"),
  paste("sampler screen: rhat <= 1.01, ess >= 400, divergences <= 10"),
  paste("herbicides:", paste(names(cleaned_fits), collapse = ", "))
), file.path(out_dir, "herbicide_provenance.txt"))

message("done: ", format(Sys.time()))
message("wrote to ", out_dir)
