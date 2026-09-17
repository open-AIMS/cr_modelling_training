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

# The backend is set here rather than left to whatever profile happens to be
# active. On 2026-09-17 this script was run twice within the hour, once with a
# profile setting cmdstanr and once without; under the rstan default the
# interaction fit failed outright with "some chains had errors", so the two runs
# were not comparable and one of them produced nothing.
options(brms.backend = "cmdstanr", mc.cores = 4)
in_dir  <- "vignettes"
out_dir <- "vignettes/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SEED <- 70

# Resumability. Each section below is skipped where the artefacts it writes are
# already present, so correcting one part of the module costs that part rather
# than the whole script. Pass --force to redo everything.
#
# This exists because it was needed. A fix to the interaction section had to be
# hand-written as a separate one-off script rather than re-run here, the
# alternative being 30 minutes of refitting to change two files.
FORCE <- "--force" %in% commandArgs(trailingOnly = TRUE)
need <- function(...) {
  targets <- file.path(out_dir, c(...))
  if (FORCE || !all(file.exists(targets))) return(TRUE)
  message("  skipping, already present: ", paste(c(...), collapse = ", "))
  FALSE
}

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

if (need("grouping_ogl_compare.csv")) {
  message("  ogl: plain")
  ogl_plain <- bnec(y ~ crf(x, model = "ecxlin"), data = colour, seed = SEED)
  message("  ogl: grouped")
  ogl_group <- bnec(y ~ crf(x, model = "ecxlin") + ogl(groupvar),
                    data = colour, seed = SEED)

  compare_two(ogl_plain, ogl_group, "no group term", "ogl(groupvar)") |>
    write.csv(file.path(out_dir, "grouping_ogl_compare.csv"), row.names = FALSE)
}

# --- 2. deviation on every parameter: bioluminescent assay plates ------------
plates <- read.csv(file.path(in_dir, "example_pgl.csv"))

if (need("grouping_pgl_compare.csv")) {
  message("  pgl: plain")
  pgl_plain <- bnec(y ~ crf(log_x, model = "nechorme"), data = plates,
                    seed = SEED)
  message("  pgl: grouped")
  pgl_group <- bnec(y ~ crf(log_x, model = "nechorme") + pgl(plateRep),
                    data = plates, seed = SEED)

  compare_two(pgl_plain, pgl_group, "no group term", "pgl(plateRep)") |>
    write.csv(file.path(out_dir, "grouping_pgl_compare.csv"), row.names = FALSE)
}

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

if (need("grouping_zinc_weights.csv")) {
  message("  zinc: pooled model set")
  zinc_all <- bnec(y ~ crf(x, model = "decline"), data = zinc, seed = SEED)
  zinc_all <- screen_models(zinc_all)

  zinc_all$mod_stats |>
    select(model, waic, wi) |>
    write.csv(file.path(out_dir, "grouping_zinc_weights.csv"), row.names = FALSE)
}

# Read back rather than kept in memory, so that the interaction section below
# runs whether or not the pooled set was refitted on this pass.
zinc_w <- read.csv(file.path(out_dir, "grouping_zinc_weights.csv"))
best_name <- zinc_w$model[which.max(zinc_w$wi)]
message("  zinc: best pooled equation is ", best_name)

# --- 4. the interaction model in brms ---------------------------------------
# The formula and the priors are read off the bayesnec fit rather than written
# out again, so that the interaction model is the same equation with the same
# priors and differs only in the right-hand side of each parameter line.
if (need("grouping_zinc_model_weights.csv", "grouping_zinc_preds.csv",
         "grouping_zinc_inter_coefs.csv", "grouping_zinc_inter_sampling.txt")) {
  # One equation refitted rather than the whole set pulled from memory, so this
  # section stands alone: it needs the selected equation's formula and priors and
  # nothing else, and refitting one model is about a minute.
  best_one <- bnec(y ~ crf(x, model = best_name), data = zinc, seed = SEED)
  best_brms <- pull_brmsfit(best_one)
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

  # Scales for the hardness contrasts, by what the parameter measures. `top`
  # and `bot` are on the response scale, which is a proportion; `ec50` is on the
  # log predictor scale, which spans about four units here; `beta` is a log
  # slope. Anything else falls back to 2, which is weak on every scale these
  # equations use.
  #
  # These are needed because bayesnec's own priors cannot be reused unchanged on
  # an interaction. Each is written with an empty `coef`, so it applies to the
  # whole population-level coefficient vector for that parameter. Under `~ 1`
  # that vector is the intercept alone; under `~ hardness` it is the intercept
  # plus the two contrasts, and each contrast would take a prior written for the
  # parameter's own value. For `ecxwb1` on these data that puts beta(5, 2) on
  # the `top` contrasts, which is zero outside (0, 1) and cannot represent a
  # level whose `top` is below the reference level's. Measured 2026-09-17 with
  # the contrasts freed: top_hardness406 is -0.179 [-0.236, -0.107], a value the
  # inherited prior excludes outright.
  contrast_sd <- c(top = 0.25, bot = 0.1, ec50 = 2, beta = 2)

  # The intercept keeps whatever bayesnec chose for that parameter; the
  # contrasts are centred on no difference between levels.
  inter_priors <- do.call(c, lapply(pars,
    function(pp) {
      row <- priors[priors$nlpar == pp & priors$class == "b" &
                      !nzchar(priors$coef), ]
      sd_pp <- if (pp %in% names(contrast_sd)) contrast_sd[[pp]] else 2
      c(brms::set_prior(row$prior[1], nlpar = pp, coef = "Intercept"),
        brms::set_prior(paste0("normal(0, ", sd_pp, ")"), nlpar = pp))
    }))

  message("  zinc: brms, no interaction")
  brm_plain <- brms::brm(bf_plain, data = zinc, family = Beta(link = "identity"),
                         prior = priors, iter = 5000, seed = 700,
                         save_pars = brms::save_pars(all = TRUE))
  # Starting values, taken from the fit without the interaction: each parameter
  # begins at its pooled estimate and every level difference begins at zero.
  #
  # Without these, three of the four chains die at initialisation and brms
  # returns the one that survived. The Beta family with an identity link needs
  # the mean inside (0, 1) at every observation, and a random start does not
  # satisfy that. bnec() supplies its own inits, which is why the model set
  # above has no such trouble; a formula handed straight to brm() does not
  # inherit them. Measured 2026-09-17: with random inits, 1 chain of 4 survived
  # and rhat was computed on that chain alone; with these inits, 4 of 4 finish
  # and rhat is 1.0012.
  fx_plain <- brms::fixef(brm_plain)
  n_lev <- nlevels(zinc$hardness)
  init_fun <- function() {
    stats::setNames(
      lapply(pars, function(pp)
        as.array(c(fx_plain[paste0(pp, "_Intercept"), "Estimate"],
                   rep(0, n_lev - 1)))),
      paste0("b_", pars))
  }

  # adapt_delta is raised because this model samples awkwardly even when it
  # starts well: `bot` sits at 0.025, close to the boundary the identity link
  # imposes. Measured 2026-09-17, that leaves 238 divergent transitions in
  # 10000 draws at adapt_delta = 0.99, which the estimates are insensitive to.
  # bnec_group() below fits each level separately and does not meet the problem.
  message("  zinc: brms, with interaction")
  brm_inter <- brms::brm(bf_inter, data = zinc, family = Beta(link = "identity"),
                         prior = inter_priors, iter = 5000, seed = 700,
                         init = init_fun, control = list(adapt_delta = 0.99),
                         save_pars = brms::save_pars(all = TRUE))

  # Chains are checked rather than assumed. brms returns whatever finished, so
  # a fit that lost three of four chains otherwise reads as a normal result.
  n_chains <- posterior::nchains(brms::as_draws_array(brm_inter))
  if (n_chains < 4) {
    stop("the interaction fit kept only ", n_chains,
         " of 4 chains; inits or adapt_delta need revisiting", call. = FALSE)
  }

  # Recorded so the module can report the sampler's behaviour rather than omit
  # it, since the course's own screen in module 4 sets divergences at 10.
  np <- brms::nuts_params(brm_inter)
  writeLines(c(
    paste("divergences:", sum(np$Parameter == "divergent__" & np$Value == 1)),
    paste("post-warmup draws:", sum(np$Parameter == "divergent__")),
    paste("max rhat:", round(max(brms::rhat(brm_inter), na.rm = TRUE), 4)),
    paste("chains:", n_chains),
    paste("adapt_delta: 0.99")
  ), file.path(out_dir, "grouping_zinc_inter_sampling.txt"))

  # The hardness contrasts, so the module can quote them rather than restate
  # them by hand.
  fx <- as.data.frame(summary(brm_inter)$fixed)
  fx$term <- rownames(fx)
  write.csv(fx[, c("term", "Estimate", "l-95% CI", "u-95% CI", "Rhat")],
            file.path(out_dir, "grouping_zinc_inter_coefs.csv"), row.names = FALSE)

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
}

# --- 5. each level fitted separately, and compared --------------------------
# bnec_group() fits the whole model set within each level, so the shape of the
# curve is free to differ between levels. That is the assumption the single
# interaction model above cannot make. On these data the interaction model
# describes all three levels adequately (RMSE 0.076, 0.055, 0.063 at hardness 5,
# 30 and 406), so this section demonstrates the route rather than rescuing a
# poor fit. The 2023 module claimed the highest level fitted badly; that was
# measured against ecxwb2, which is not the equation the weights select.
if (need("grouping_zinc_ne_posterior.csv", "grouping_zinc_prob_diff.csv",
         "grouping_zinc_diff.csv")) {
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
}

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
