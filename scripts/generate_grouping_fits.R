# ---------------------------------------------------------------------------
# Fits the grouped and factor-covariate models used by module 8 and writes the
# compact artefacts that module reads.
#
# None of this runs at render time. The `bnec_group()` calls alone are complete
# model sets fitted at each level of a factor, and the two `brms` interaction
# models are slow on top of that. A workshop participant waiting for a chunk has
# stopped learning, so the fitting happens here and the module reads the results.
#
# Run with:
#   Rscript scripts/generate_grouping_fits.R
#
# Expect this to take hours. Run it deliberately, not as part of a render.
#
# One fit is not made here. The copper-against-zinc `bnec_group()` call is the
# call `vignette("example8")` in `bayesnec` makes, and that vignette's fits are
# run as cluster array tasks by the `open-AIMS/grouping-structures` compendium
# and kept in a store. Set LUM31_TOX_FIT to the stored `.rds` and this script
# reads it instead of sampling; leave it unset and the call is fitted here. The
# store is keyed on the call and a digest of the data, so a stored fit answers a
# call only where both match.
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

# No future plan is set, and that is deliberate. `bnec()` fits a model set under
# whatever plan is active (bayesnec #184), and under `multisession` on this
# machine on 2026-09-17 the copper plate set deadlocked: the parent had 439
# bytes queued to one worker and that worker 415 bytes queued back, both
# polling, with no Stan process running and no CPU used by any of the six
# processes for fourteen minutes. Every fit below therefore runs in sequence,
# with its four chains in parallel through `mc.cores`.

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

# The `lum31` data set reaches `bayesnec` through pull request #228 and is not
# on `dev` until that merges. Checked here rather than at the call, because the
# first section takes half an hour and the failure would otherwise arrive after
# it. `data()` returns the name it loaded whether or not it found anything, so
# the test is on the object.
if (!"lum31" %in% data(package = "bayesnec")$results[, "Item"]) {
  stop("this bayesnec has no `lum31` data set. Install the branch of ",
       "open-AIMS/bayesnec#228, or dev once that has merged.", call. = FALSE)
}

message("fitting: ", format(Sys.time()))
t0 <- Sys.time()

# Compares fits of the same data by leave-one-out cross-validation, and returns
# the elpd difference and the model weights as one row per fit. Both quantities
# are reported because they answer different questions: the weight says how much
# of an averaged prediction each fit would contribute, and the elpd difference
# says by how much and with what standard error.
#
# `fits` is a named list, and the names are what the module's table reports, so
# they name the term rather than the object holding it.
compare_fits <- function(fits) {
  brmsfits <- lapply(fits, pull_brmsfit)
  wts <- do.call(brms::loo_model_weights,
                 c(unname(brmsfits), list(method = "pseudobma")))
  loos <- lapply(brmsfits, brms::loo)
  # loo_compare() is handed the named list rather than the objects one by one.
  # Supplied individually it labels its rows from the deparsed arguments, which
  # here are `unname(loos)[[1]]` and not the term being compared; from a named
  # list the names become the row names and survive the sort by elpd.
  cmp <- as.data.frame(loo::loo_compare(loos))
  cmp$fit <- rownames(cmp)
  # Pareto k is reported rather than assumed acceptable. A censored Gamma
  # likelihood with a group-level term over sixteen levels is the shape that
  # produces influential observations, and an elpd difference computed over
  # unreliable pointwise values is not a comparison.
  bad_k <- vapply(loos, function(l) sum(l$diagnostics$pareto_k > 0.7), numeric(1))
  data.frame(
    fit = names(fits),
    weight = round(as.numeric(wts), 3),
    elpd_diff = round(cmp$elpd_diff[match(names(fits), cmp$fit)], 1),
    se_diff = round(cmp$se_diff[match(names(fits), cmp$fit)], 1),
    pareto_k_over_0.7 = as.integer(bad_k)
  )
}

# The plotting frame behind `autoplot()`, stacked over several fits with a
# column naming which fit each row came from. The module draws its own figures
# from these, because it does not hold the fitted objects at render time.
#
# `xform` is not used on any call below. `ggbnec_data()` returns the predictor
# on the scale it was recorded on, having already inverted a transformation
# written inside `crf()`, so a figure drawn from this frame is already in the
# recorded units and the axis is transformed instead.
curve_frame <- function(fits, group = NULL) {
  keep <- c("x_e", "y_e", "y_ci", "x_r", "y_r", "nec_vals", "group", "panel")
  do.call(rbind, Map(function(f, nm) {
    d <- if (is.null(group)) ggbnec_data(f) else ggbnec_data(f, group = group)
    data.frame(fit = nm, d[, intersect(names(d), keep)])
  }, fits, names(fits)))
}

# The curve of a single-equation fit, taken from the fit's own prediction grid.
#
# `ggbnec_data()` has two sources and they differ by a factor of ten in
# resolution. A `bayesmanecfit`'s curve comes from `w_pred_vals`, which `bnec()`
# computed on its own grid at `resolution = 1000`; a `bayesnecfit`'s comes from
# `brms::conditional_effects()`, whose default resolution is 100. Both grids are
# evenly spaced on the *recorded* predictor, not on the scale the model was
# fitted on.
#
# For a fit written as `crf(log(conc), ...)` and drawn on a log axis that
# decides the figure. Over the copper range, 0.002 to 10.2 mg/L, a 1000-point
# grid puts 10 points below 0.1 and a 100-point grid puts 1, so the single fit's
# curve was drawn as one straight segment across the whole lower decade and
# appeared to decline steadily where `ecxll4` is flat. Measured 2026-09-17: the
# `bayesnecfit` frame held 100 grid points with one below 0.1 mg/L, and the
# `bayesmanecfit` frame for the same data held 1000 with ten.
#
# `pred_vals$data` is the 1000-point grid and is on the fit, so this uses it.
single_curve <- function(fits) {
  do.call(rbind, Map(function(f, nm) {
    d <- f$pred_vals$data
    data.frame(fit = nm, x = d$x, Estimate = d$Estimate,
               Q2.5 = d$Q2.5, Q97.5 = d$Q97.5)
  }, fits, names(fits)))
}

# The toxicity estimates a module reports, on the concentration scale. `xform`
# is needed because `ecx()` and `nsec()` return the scale the model was fitted
# on, which is logged in every call below.
#
# `ecx_vals` is a vector because one comparison needs an EC50: @luter2025 report
# EC50s for these data and nothing else, so an EC50 is what the published values
# can be read against. Each additional value is another pass of `ecx()` over
# every equation of the set, which is the cost bayesnec #306 measured, so the
# default asks for one.
estimate_frame <- function(fits, xform = exp, ecx_vals = 10) {
  do.call(rbind, Map(function(f, nm) do.call(rbind, c(
    list(data.frame(fit = nm, estimate = "NSEC", t(nsec(f, xform = xform)))),
    lapply(ecx_vals, function(v)
      data.frame(fit = nm, estimate = paste0("EC", v),
                 t(ecx(f, ecx_val = v, xform = xform))))
  )), fits, names(fits)))
}

# --- 1. deviation on the intercept: coral colour, tank offsets ---------------
colour <- read.csv(file.path(in_dir, "example_ogl.csv")) |>
  mutate(y = Perc_change / 100,
         groupvar = factor(paste(Replicate, Treatment, sep = "_")),
         x = as.numeric(Treatment)) |>
  na.omit()

if (need("grouping_ogl_compare.csv", "grouping_ogl_curves.csv")) {
  message("  ogl: plain")
  ogl_plain <- bnec(y ~ crf(x, model = "ecxlin"), data = colour, seed = SEED)
  message("  ogl: grouped")
  ogl_group <- bnec(y ~ crf(x, model = "ecxlin") + ogl(groupvar),
                    data = colour, seed = SEED)

  compare_fits(list(`no group term` = ogl_plain, `ogl(groupvar)` = ogl_group)) |>
    write.csv(file.path(out_dir, "grouping_ogl_compare.csv"), row.names = FALSE)

  # The grouped frame carries the grouping, so the figure can mark the mean of
  # each container. The ungrouped fit has no grouping to carry, and its `group`
  # column is filled in afterwards so that the two frames stack.
  rbind(
    transform(curve_frame(list(`no group term` = ogl_plain)),
              group = NA_character_),
    curve_frame(list(`ogl(groupvar)` = ogl_group), group = "groupvar")) |>
    write.csv(file.path(out_dir, "grouping_ogl_curves.csv"), row.names = FALSE)
}

# --- 2. copper against zinc: a factor covariate on the assay ----------------
# Two reference toxicants, run repeatedly to establish that the assay is
# reproducible. They are two contaminants of interest rather than a sample from
# a population of metals, so they are levels and not a grouping, and
# `bnec_group()` fits each separately.
#
# `model = "all"` rather than `decline`, because zinc rises above its control at
# low concentrations and the `decline` set holds only equations whose curves fall
# from the control. Excluding the hormesis equations would be a decision about
# the shape of the zinc curve taken before fitting.
data(lum31)
zn_cu_15 <- droplevels(subset(lum31, minutes == 15))
gfam <- Gamma(link = "identity")

if (need("grouping_tox_weights.csv", "grouping_tox_estimates.csv",
         "grouping_tox_prob_diff.csv", "grouping_tox_posterior.csv",
         "grouping_tox_curves.csv")) {
  store <- Sys.getenv("LUM31_TOX_FIT")
  if (nzchar(store)) {
    message("  toxicant: reading the stored fit from ", store)
    fits_tox <- readRDS(store)
    if (!inherits(fits_tox, "bayesnecgroupfit")) {
      stop("LUM31_TOX_FIT is not a bayesnecgroupfit", call. = FALSE)
    }
  } else {
    # The sampler settings and the seed are the vignette's rather than this
    # script's SEED, so that a fit made here and a fit read from the store are
    # the same call. The store is keyed on the call, so they would otherwise be
    # two different analyses reported under one set of numbers.
    message("  toxicant: bnec_group across copper and zinc")
    fits_tox <- bnec_group(
      rlu_cens | cens(censoring) ~ crf(log(conc), "all") + ogl(plate) +
        disp("power"),
      data = zn_cu_15, group_var = "toxicant", family = gfam,
      chains = 4, iter = 5000, warmup = 4000, seed = 228,
      control = list(adapt_delta = 0.99))
  }
  tox_fits <- setNames(fits_tox$fits, fits_tox$levels)

  # Model weights at each level, which is what says whether the two metals are
  # described by the same equation. `crossed_group_weights()` answers that over
  # the whole table at once; the per-level weights are what the module shows,
  # because the disagreement between levels is visible in them directly.
  do.call(rbind, Map(function(f, lev) {
    w <- summary(f)$mod_weights[, "wi"]
    data.frame(toxicant = lev, model = names(w), wi = round(as.numeric(w), 3))
  }, fits_tox$fits, fits_tox$levels)) |>
    write.csv(file.path(out_dir, "grouping_tox_weights.csv"), row.names = FALSE)

  estimate_frame(tox_fits, ecx_vals = c(10, 50)) |>
    rename(toxicant = fit) |>
    write.csv(file.path(out_dir, "grouping_tox_estimates.csv"), row.names = FALSE)

  cmp_tox <- compare_posterior(fits_tox, comparison = "nec")
  write.csv(cmp_tox$prob_diff, file.path(out_dir, "grouping_tox_prob_diff.csv"),
            row.names = FALSE)
  # Returned on the fitted scale, log(mg/L), and back-transformed here so that
  # the module's density figure reads in the units the estimates are reported in.
  tox_post <- cmp_tox$posterior_data
  tox_post$value <- exp(tox_post$value)
  write.csv(tox_post, file.path(out_dir, "grouping_tox_posterior.csv"),
            row.names = FALSE)

  write.csv(curve_frame(tox_fits),
            file.path(out_dir, "grouping_tox_curves.csv"), row.names = FALSE)
}

# --- 3. the plate term on the copper arm ------------------------------------
# The recorded bioluminescence, rather than the normalised response this module
# used until 2026-09-17. That response was each plate divided by its own four
# control wells and then by the plate maximum, which is the practice Ritz,
# Gerhard and Streibig (2026) measure the bias and coverage of, and the file
# holding it covered four plates. `lum31` holds all sixteen copper plates as
# recorded.
#
# Copper rather than zinc, and 15 minutes rather than 30. Copper's mean response
# falls below the control at every concentration on every one of the sixteen
# plates, so no hormesis equation is needed; the 15-minute read has 106 of its
# 704 readings at the instrument floor against 197 at 30 minutes, so the
# censoring bound decides less of the answer.
cu15 <- droplevels(subset(lum31, toxicant == "Cu" & minutes == 15))

# The equation is taken from the copper level of the fit above rather than from
# a model set fitted again here. That fit is the same 704 readings with the same
# plate term and the same dispersion sub-model, so a second set would be the
# same sampling twice: eleven equations of about twenty minutes each on this
# machine. The module says where the weights come from.
plate_best <- {
  w <- read.csv(file.path(out_dir, "grouping_tox_weights.csv"))
  w <- w[w$toxicant == "Cu", ]
  w$model[which.max(w$wi)]
}
message("  plate: best copper equation under ogl(plate) is ", plate_best)

# The three fits are kept rather than discarded. They are the slowest thing this
# script does, about half an hour each, and every artefact below is a summary of
# them: a change to how the curve is drawn or which estimate is reported should
# not cost a refit. `ignore/` is git-ignored, so this is a local cache and not a
# second copy of the analysis.
fit_cache <- "ignore/grouping_fits"
dir.create(fit_cache, showWarnings = FALSE, recursive = TRUE)
plate_rds <- file.path(fit_cache, "plate_fits.rds")

if (need("grouping_plate_compare.csv", "grouping_plate_curves.csv",
         "grouping_plate_estimates.csv", "grouping_plate_sd.csv",
         "grouping_plate_obs.csv")) {
  if (file.exists(plate_rds)) {
    message("  plate: reading the three fits from ", plate_rds)
    plate_fits <- readRDS(plate_rds)
  } else {
    # One equation, fitted three ways. What is being compared is the group-level
    # term and not the equation, so the equation is held fixed. Comparing two
    # model averages by leave-one-out would confound the two.
    message("  plate: ", plate_best, ", no group term")
    plate_plain <- bnec(rlu_cens | cens(censoring) ~ crf(log(conc), plate_best) +
                          disp("power"),
                        data = cu15, family = gfam, seed = SEED)
    message("  plate: ", plate_best, ", ogl(plate)")
    plate_ogl <- bnec(rlu_cens | cens(censoring) ~ crf(log(conc), plate_best) +
                        ogl(plate) + disp("power"),
                      data = cu15, family = gfam, seed = SEED)
    message("  plate: ", plate_best, ", pgl(plate)")
    plate_pgl <- bnec(rlu_cens | cens(censoring) ~ crf(log(conc), plate_best) +
                        pgl(plate) + disp("power"),
                      data = cu15, family = gfam, seed = SEED)

    plate_fits <- list(`no group term` = plate_plain,
                       `ogl(plate)` = plate_ogl,
                       `pgl(plate)` = plate_pgl)
    saveRDS(plate_fits, plate_rds)
  }

  compare_fits(plate_fits) |>
    write.csv(file.path(out_dir, "grouping_plate_compare.csv"), row.names = FALSE)

  estimate_frame(plate_fits) |>
    write.csv(file.path(out_dir, "grouping_plate_estimates.csv"), row.names = FALSE)

  single_curve(plate_fits) |>
    write.csv(file.path(out_dir, "grouping_plate_curves.csv"), row.names = FALSE)

  # The observations, written once. They are the same in all three fits, so
  # repeating them per fit would triple the file to say the same thing.
  cu15[, c("conc", "rlu_cens", "plate")] |>
    write.csv(file.path(out_dir, "grouping_plate_obs.csv"), row.names = FALSE)

  # The group-level standard deviations. `ogl()` estimates one; `pgl()`
  # estimates one for every parameter of the equation, which is what the module
  # compares them on. `summary()` does not report them, so they are read from
  # the underlying brmsfit.
  grouped <- plate_fits[c("ogl(plate)", "pgl(plate)")]
  do.call(rbind, Map(function(f, nm) {
    v <- brms::VarCorr(pull_brmsfit(f))$plate$sd
    data.frame(fit = nm, term = rownames(v), round(as.data.frame(v), 3),
               row.names = NULL)
  }, grouped, names(grouped))) |>
    write.csv(file.path(out_dir, "grouping_plate_sd.csv"), row.names = FALSE)
}

# --- 4. the zinc by hardness factor covariate -------------------------------
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

# --- 5. the interaction model in brms ---------------------------------------
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
  # imposes. It still leaves divergent transitions, which the estimates are
  # insensitive to; the count is written to grouping_zinc_inter_sampling.txt
  # and the module prints that file rather than quoting a number from here.
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

# --- 6. each level fitted separately, and compared --------------------------
# bnec_group() fits the whole model set within each level, so the shape of the
# curve is free to differ between levels. That is the assumption the single
# interaction model above cannot make, and the per-level weights written below
# are what says whether it holds on these data: an interaction is a reasonable
# description where every level puts its weight on the same equation.
#
# A per-level residual root mean square stood here until 2026-09-17, quoted in
# the module as 0.068, 0.052 and 0.057 and in this comment as 0.076, 0.055 and
# 0.063. Neither was produced by anything that ran, the two disagreed, and the
# weights answer the question the number was standing in for. The 2023 module
# claimed the highest level fitted badly; that was measured against ecxwb2,
# which is not the equation the weights select.
if (need("grouping_zinc_ne_posterior.csv", "grouping_zinc_prob_diff.csv",
         "grouping_zinc_diff.csv", "grouping_zinc_group_weights.csv",
         "grouping_zinc_group_estimates.csv",
         "grouping_zinc_group_curves.csv")) {
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

  # Which equation describes each level, which is the question the single
  # interaction model cannot ask: it holds the functional form fixed and lets
  # only the parameter values differ.
  do.call(rbind, Map(function(f, lev) {
    w <- summary(f)$mod_weights[, "wi"]
    data.frame(hardness = lev, model = names(w), wi = round(as.numeric(w), 3))
  }, zinc_group$fits, zinc_group$levels)) |>
    write.csv(file.path(out_dir, "grouping_zinc_group_weights.csv"),
              row.names = FALSE)

  # The predictor is logged in the data rather than inside crf(), so the
  # estimates are on the log scale and exp() returns them to concentration.
  estimate_frame(setNames(zinc_group$fits, zinc_group$levels)) |>
    rename(hardness = fit) |>
    write.csv(file.path(out_dir, "grouping_zinc_group_estimates.csv"),
              row.names = FALSE)

  # The fitted curves, so the module can put this route's figure beside the
  # interaction model's. The predictor was logged in the data rather than inside
  # crf(), so these are on the log scale and the module's axis matches the
  # interaction figure's without back-transformation.
  write.csv(curve_frame(setNames(zinc_group$fits, zinc_group$levels)),
            file.path(out_dir, "grouping_zinc_group_curves.csv"),
            row.names = FALSE)
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
  paste("best copper equation under ogl(plate):", plate_best),
  paste("copper plates:", nlevels(cu15$plate), "at 15 minutes,",
        nrow(cu15), "readings,", sum(cu15$censoring == "left"), "left-censored"),
  paste("copper against zinc:",
        if (nzchar(Sys.getenv("LUM31_TOX_FIT")))
          paste("read from", Sys.getenv("LUM31_TOX_FIT")) else "fitted here"),
  paste("pooled zinc set: decline,", length(models()$decline), "equations"),
  paste("best pooled equation:", best_name),
  paste("hardness levels:", paste(levels(zinc$hardness), collapse = ", ")),
  "sampler screen: rhat <= 1.01, ess >= 400, divergences <= 10"
), file.path(out_dir, "grouping_provenance.txt"))

message("done: ", format(Sys.time()))
message("wrote to ", out_dir)
