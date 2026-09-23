# ---------------------------------------------------------------------------
# Fits every admissible equation to each module 5 dataset, to choose which
# equations the module's worked examples use.
#
# Module 5 teaches matching a family to a response, and each example was fitted
# with two equations (an ecx and a nec) chosen without looking at the data. On
# 2026-09-23 RF found that several of the resulting check_fit() plots showed a
# poor fit, and asked whether that was the family or the equation. This script
# answers it: for each example it fits model = "all" under the family the
# module uses, runs check_sampling() and check_fit() on every equation, and
# records the weight each receives. The module is not changed by this script;
# it records the evidence the choice of equations is made from.
#
# Output, per example:
#   notes/module5-model-screen/<example>.csv   one row per equation
#   ignore/module5_screen/<example>.rds        the fitted set (local only)
#
# Run one example, or all of them in sequence:
#   Rscript scripts/screen_module5_models.R gamma
#   Rscript scripts/screen_module5_models.R
# An example whose CSV exists is skipped unless --force is given.
# ---------------------------------------------------------------------------

suppressMessages({
  library(bayesnec)
  library(dplyr)
})

set.seed(333)
options(brms.backend = "cmdstanr", mc.cores = 4)

args <- commandArgs(trailingOnly = TRUE)
FORCE <- "--force" %in% args
which_examples <- setdiff(args, "--force")

source("scripts/module5_data.R")

csv_dir <- "notes/module5-model-screen"
rds_dir <- "ignore/module5_screen"
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# Each example as the module fits it, with the equation set widened to "all".
# bnec() drops the equations a family cannot represent and records why, so
# "all" is the full admissible set for that family and response.
examples <- list(
  binomial = function() bnec(suc | trials(tot) ~ crf(log.x, model = "all"),
                             data = binom_data, seed = 333),
  betabinomial = function() bnec(suc | trials(tot) ~ crf(log.x, model = "all"),
                                 data = binom_data, family = "beta_binomial",
                                 seed = 333),
  beta = function() bnec(yield ~ crf(log.x, model = "all"),
                         data = pam_data, seed = 333),
  poisson = function() bnec(density ~ crf(log.x, model = "all"),
                            data = alga_data_a, family = "poisson", seed = 333),
  negbinomial = function() bnec(density ~ crf(log.x, model = "all"),
                                data = alga_data_a, family = "negbinomial",
                                seed = 333),
  gamma = function() bnec(rlu ~ crf(log.x, model = "all"),
                          data = lum_data, family = "Gamma", seed = 333),
  gaussian = function() bnec(sgr ~ crf(log.x, model = "all"),
                             data = alga_data_a, seed = 333),
  # Alternatives screened after the examples above, on 2026-09-23. Every
  # equation failed the scale check on the alga contaminant A counts and growth
  # rates (alga_data_a above), and on the coral_pam yield, so these test
  # whether the dataset or the dispersion is the cause rather than the
  # equation. alga_data is contaminant B, which module 5 now fits.
  poisson_b = function() bnec(density ~ crf(log.x, model = "all"),
                              data = alga_data, family = "poisson", seed = 333),
  negbinomial_b = function() bnec(density ~ crf(log.x, model = "all"),
                                  data = alga_data, family = "negbinomial",
                                  seed = 333),
  gaussian_b = function() bnec(sgr ~ crf(log.x, model = "all"),
                               data = alga_data, seed = 333),
  # The three equations that describe the contaminant B growth rate best all
  # failed check_sampling() on divergent transitions at the default
  # adapt_delta; refitted here with a smaller step size.
  gaussian_b_ad = function() bnec(sgr ~ crf(log.x, model = c("ecxwb1", "ecx4param",
                                                              "ecxll4")),
                                  data = alga_data, seed = 333,
                                  control = list(adapt_delta = 0.99)),
  # The best-weighted beta equation with the precision modelled on the fitted
  # mean. Two equations are fitted because summarise_set() expects a set.
  beta_disp = function() bnec(yield ~ crf(log.x, model = c("ecxwb1", "ecxll3")) +
                                disp("power"),
                              data = pam_data, seed = 333),
  # adapt_delta as in the module's hurdle fit, where the survival block is
  # determined by the top two doses and samples awkwardly at the default.
  hurdle = function() bnec_hurdle(growth ~ crf(log_dose, model = "all"),
                                  data = snail, seed = 333,
                                  control = list(adapt_delta = 0.99))
)

if (length(which_examples)) {
  unknown <- setdiff(which_examples, names(examples))
  if (length(unknown)) stop("unknown example: ", paste(unknown, collapse = ", "))
  examples <- examples[which_examples]
}

# check_fit() flags a group when the observed statistic lies outside the
# central 95% of its simulated distribution, which is what the bars in its plot
# show. The same cut-off is used here so the counts match the plots.
outside <- function(p) sum(p < 0.025 | p > 0.975, na.rm = TRUE)

summarise_set <- function(fit, example) {
  # A hurdle fit holds two sets; the growth block is the one module 5 reads
  # against its family, so that is the one screened.
  set <- if (inherits(fit, "bayesnechurdlefit")) fit$growth else fit
  if (!inherits(set, "bayesmanecfit")) {
    stop("expected a model set for ", example, "; only one equation fitted")
  }

  cf <- check_fit(set, plot = FALSE)
  fit_summary <- cf |>
    group_by(model) |>
    summarise(n_groups = n(),
              location_outside = outside(ppp_mean),
              scale_outside = outside(ppp_sd),
              max_abs_log_sd_ratio = max(abs(log(sd_ratio)), na.rm = TRUE),
              .groups = "drop")

  cs <- check_sampling(set)
  stats <- data.frame(model = set$mod_stats$model,
                      waic = set$mod_stats$waic,
                      wi_all = set$mod_stats$wi)

  # Weights recomputed over the equations that pass check_sampling(), which
  # is what a set is reported on after screen_models().
  screened <- tryCatch(screen_models(set), error = function(e) NULL)
  wi_screened <- if (is.null(screened) || !inherits(screened, "bayesmanecfit")) {
    data.frame(model = character(), wi_screened = numeric())
  } else {
    data.frame(model = screened$mod_stats$model,
               wi_screened = screened$mod_stats$wi)
  }

  out <- stats |>
    left_join(cs, by = "model") |>
    left_join(fit_summary, by = "model") |>
    left_join(wi_screened, by = "model") |>
    mutate(example = example, .before = 1) |>
    arrange(desc(wi_all))

  rec <- bnec_record(set)
  attr(out, "excluded") <- rec$excluded
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

for (ex in names(examples)) {
  csv <- file.path(csv_dir, paste0(ex, ".csv"))
  if (!FORCE && file.exists(csv)) {
    message("skip ", ex, " (", csv, " exists)")
    next
  }
  message("fitting ", ex, ": ", format(Sys.time()))
  t0 <- Sys.time()
  fit <- examples[[ex]]()
  saveRDS(fit, file.path(rds_dir, paste0(ex, ".rds")))
  res <- summarise_set(fit, ex)
  res$minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  write.csv(res, csv, row.names = FALSE)
  excl <- attr(res, "excluded")
  if (NROW(excl)) {
    write.csv(excl, file.path(csv_dir, paste0(ex, "_excluded.csv")),
              row.names = FALSE)
  }
  message("  done ", ex, " in ", res$minutes[1], " min")
  print(res[, c("model", "wi_all", "wi_screened", "failed",
                "location_outside", "scale_outside", "n_groups")])
}
