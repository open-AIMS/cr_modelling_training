# ---------------------------------------------------------------------------
# Fits the grouped and factor-covariate models used by module 8 and writes the
# compact artefacts that module reads.
#
# None of this runs at render time. The `bnec_group()` call alone is the
# `decline` set fitted separately at each of three hardness levels, and the two
# `brms` interaction models are slow on top of that. A workshop participant
# waiting for a chunk has stopped learning, so the fitting happens here and the
# module reads the results.
#
# Run with:
#   Rscript scripts/generate_grouping_fits.R
#
# Expect this to take hours. Run it deliberately, not as part of a render.
# ---------------------------------------------------------------------------

suppressMessages({
  library(bayesnec)
  library(brms)
  library(dplyr)
  library(purrr)
})

options(mc.cores = 4)
in_dir  <- "vignettes"
out_dir <- "vignettes/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SEED <- 70

message("fitting: ", format(Sys.time()))
t0 <- Sys.time()

# Compares two fits of the same data by leave-one-out cross-validation, and
# returns the elpd difference and the model weights as one row per fit. Both
# quantities are reported because they answer different questions: the weight
# says how much of an averaged prediction each fit would contribute, and the
# elpd difference says by how much and with what standard error.
compare_two <- function(fit_a, fit_b, label_a, label_b) {
  a <- pull_brmsfit(fit_a)
  b <- pull_brmsfit(fit_b)
  wts <- brms::loo_model_weights(a, b, method = "pseudobma")
  cmp <- as.data.frame(brms::loo_compare(brms::loo(a), brms::loo(b)))
  cmp$fit <- rownames(cmp)
  cmp$fit <- c(label_a, label_b)[match(cmp$fit, c("a", "b"))]
  data.frame(
    fit = c(label_a, label_b),
    weight = round(as.numeric(wts), 3),
    elpd_diff = round(cmp$elpd_diff[match(c(label_a, label_b), cmp$fit)], 1),
    se_diff = round(cmp$se_diff[match(c(label_a, label_b), cmp$fit)], 1)
  )
}

# --- 1. deviation on the intercept: coral colour, tank offsets ---------------
colour <- read.csv(file.path(in_dir, "example_ogl.csv")) |>
  mutate(y = Perc_change / 100,
         groupvar = factor(paste(Replicate, Treatment, sep = "_")),
         x = as.numeric(Treatment)) |>
  na.omit()

message("  ogl: plain")
ogl_plain <- bnec(y ~ crf(x, model = "ecxlin"), data = colour, seed = SEED)
message("  ogl: grouped")
ogl_group <- bnec(y ~ crf(x, model = "ecxlin") + ogl(groupvar),
                  data = colour, seed = SEED)

compare_two(ogl_plain, ogl_group, "no group term", "ogl(groupvar)") |>
  write.csv(file.path(out_dir, "grouping_ogl_compare.csv"), row.names = FALSE)

# --- 2. deviation on every parameter: bioluminescent assay plates ------------
plates <- read.csv(file.path(in_dir, "example_pgl.csv"))

message("  pgl: plain")
pgl_plain <- bnec(y ~ crf(log_x, model = "nechorme"), data = plates, seed = SEED)
message("  pgl: grouped")
pgl_group <- bnec(y ~ crf(log_x, model = "nechorme") + pgl(plateRep),
                  data = plates, seed = SEED)

compare_two(pgl_plain, pgl_group, "no group term", "pgl(plateRep)") |>
  write.csv(file.path(out_dir, "grouping_pgl_compare.csv"), row.names = FALSE)

# --- 3. the zinc by hardness factor covariate -------------------------------
# The response is already a percentage of each experiment's own control, so the
# division module 2 warns about has happened before these data were supplied.
# It is divided again by the observed maximum and shifted off 0 and 1 because
# brms is called directly below and does not make the adjustment bnec() makes.
zinc <- read.csv(file.path(in_dir, "example_fi.csv")) |>
  na.omit() |>
  mutate(y = (Pcon / max(Pcon) + 0.01) * 0.99,
         x = log(Mean.Zn),
         hardness = factor(round(Hardness))) |>
  select(x, y, hardness)

message("  zinc: pooled model set")
zinc_all <- bnec(y ~ crf(x, model = "decline"), data = zinc, seed = SEED)
zinc_all <- screen_models(zinc_all)

zinc_all$mod_stats |>
  select(model, waic, wi) |>
  write.csv(file.path(out_dir, "grouping_zinc_weights.csv"), row.names = FALSE)

best_name <- zinc_all$mod_stats$model[which.max(zinc_all$mod_stats$wi)]
message("  zinc: best pooled equation is ", best_name)

# --- 4. the interaction model in brms ---------------------------------------
# The formula and the priors are read off the bayesnec fit rather than written
# out again, so that the interaction model is the same equation with the same
# priors and differs only in the right-hand side of each parameter line.
best_brms <- pull_brmsfit(zinc_all, model = best_name)
priors <- brms::prior_summary(best_brms)

# The equation is taken from the fit that the weights selected, rather than
# written out here. Writing it out is how the two came apart: the literal
# formula this script carried was `ecxwb2`, inherited from the 2023 version of
# the module, while the weights on these data select `ecxwb1`. Both have the
# same four parameters, so the priors matched by name and nothing errored; the
# interaction was simply fitted with a different equation from the one reported
# as best.
bf_plain <- stats::formula(best_brms)
pars <- names(bf_plain$pforms)
bf_inter <- do.call(brms::bf,
                    c(list(bf_plain$formula),
                      lapply(pars, function(p)
                        stats::as.formula(paste(p, "~ hardness"))),
                      list(nl = TRUE)))
message("  zinc: interaction on ", best_name, ", parameters ",
        paste(pars, collapse = ", "))

message("  zinc: brms, no interaction")
brm_plain <- brms::brm(bf_plain, data = zinc, family = Beta(link = "identity"),
                       prior = priors, iter = 5000, seed = 700,
                       save_pars = brms::save_pars(all = TRUE))
message("  zinc: brms, with interaction")
brm_inter <- brms::brm(bf_inter, data = zinc, family = Beta(link = "identity"),
                       prior = priors, iter = 5000, seed = 700,
                       save_pars = brms::save_pars(all = TRUE))

data.frame(
  fit = c("no interaction", "hardness interaction"),
  weight = round(as.numeric(
    brms::loo_model_weights(brm_plain, brm_inter, method = "pseudobma")), 3)
) |>
  write.csv(file.path(out_dir, "grouping_zinc_model_weights.csv"),
            row.names = FALSE)

# Fitted curves per level, on a grid spanning the observed predictor range.
grid <- expand.grid(
  x = seq(min(zinc$x), max(zinc$x), length.out = 100),
  hardness = levels(zinc$hardness))
preds <- cbind(grid, brms::posterior_epred(brm_inter, newdata = grid) |>
                 apply(2, quantile, probs = c(0.5, 0.025, 0.975)) |>
                 t() |>
                 as.data.frame() |>
                 setNames(c("Estimate", "Q2.5", "Q97.5")))
write.csv(preds, file.path(out_dir, "grouping_zinc_preds.csv"),
          row.names = FALSE)

# --- 5. each level fitted separately, and compared --------------------------
# bnec_group() fits the whole model set within each level, so the shape of the
# curve is free to differ between levels. That is the assumption the single
# interaction model above cannot make, and on these data it is the one that
# matters: the highest hardness level is not described well by the equation the
# pooled data selected.
message("  zinc: bnec_group across hardness")
zinc_group <- bnec_group(y ~ crf(x, model = "decline"), data = zinc,
                         group_var = "hardness", seed = SEED)

post_comp <- compare_posterior(zinc_group, comparison = "n(s)ec")
post_fitted <- compare_posterior(zinc_group, comparison = "fitted")

post_comp$posterior_data |>
  write.csv(file.path(out_dir, "grouping_zinc_ne_posterior.csv"),
            row.names = FALSE)
write.csv(post_comp$prob_diff,
          file.path(out_dir, "grouping_zinc_prob_diff.csv"), row.names = FALSE)
write.csv(post_fitted$diff_data,
          file.path(out_dir, "grouping_zinc_diff.csv"), row.names = FALSE)

# --- provenance -------------------------------------------------------------
writeLines(c(
  "Generated by scripts/generate_grouping_fits.R",
  paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("elapsed:", format(round(difftime(Sys.time(), t0), 1))),
  paste("seed:", SEED),
  paste("R:", getRversion()),
  paste("bayesnec:", as.character(packageVersion("bayesnec"))),
  paste("brms:", as.character(packageVersion("brms"))),
  paste("cmdstanr:", as.character(packageVersion("cmdstanr"))),
  paste("pooled zinc set: decline,", length(models()$decline), "equations"),
  paste("best pooled equation:", best_name),
  paste("hardness levels:", paste(levels(zinc$hardness), collapse = ", ")),
  "sampler screen: rhat <= 1.01, ess >= 400, divergences <= 10"
), file.path(out_dir, "grouping_provenance.txt"))

message("done: ", format(Sys.time()))
message("wrote to ", out_dir)
