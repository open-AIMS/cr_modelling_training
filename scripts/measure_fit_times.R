# ---------------------------------------------------------------------------
# Records how long each saved fit sampled for, into scripts/fit_times.csv.
#
# scripts/generate_live_scripts.R uses this to decide which fit calls a
# participant can run in the room. Without it every fit call is commented out,
# which is the safe default but a poor one: most of these fits sample in
# seconds and the participant never sees bnec() work.
#
# The sampling time is read from the saved object rather than timed here, so
# this costs seconds rather than the half hour a refit would take.
# rstan::get_elapsed_time gives warmup and sampling per chain. The chains run
# at the same time, so the fit takes as long as its slowest chain, and a model
# set takes the sum over its models because bnec() fits them in turn.
#
# What is NOT in the object is Stan compilation, which is most of the wall time
# of a small fit: module 2's fit samples in 2.8 seconds and the module tells the
# participant it took a little over two minutes. The generator adds an estimate
# for compilation, measured separately; see COMPILE_SECONDS there.
#
# The output is tracked, while vignettes/fits/ is git-ignored, so that
# generating the live scripts on a clone without the bundle gives the same
# result as generating them here. Commit it with the rebuilt archive.
#
# Run with:
#   Rscript scripts/measure_fit_times.R
#
# Re-run whenever the fits are rebuilt.
# ---------------------------------------------------------------------------

suppressMessages(library(bayesnec))

fit_dir <- file.path("vignettes", "fits")
stopifnot(dir.exists(fit_dir))

# The wall time of one model: its slowest chain, warm-up included.
model_seconds <- function(nf) {
  et <- tryCatch(rstan::get_elapsed_time(nf$fit$fit), error = function(e) NULL)
  if (is.null(et)) return(NA_real_)
  max(rowSums(et))
}

rows <- list()
for (f in sort(list.files(fit_dir, pattern = "\\.RData$"))) {
  e <- new.env()
  for (nm in load(file.path(fit_dir, f), envir = e)) {
    x <- get(nm, envir = e)
    # A hurdle fit is two model sets, a response block and a survival block,
    # and both are sampled. Counting only one of them halves the estimate the
    # live-script generator judges the call against.
    seconds_of <- function(z) {
      if (inherits(z, "bayesmanecfit")) {
        vapply(z$mod_fits, model_seconds, numeric(1))
      } else if (inherits(z, "bayesnecfit")) {
        model_seconds(z)
      } else {
        numeric(0)
      }
    }
    per <- if (inherits(x, c("bayesnechurdlefit", "bayesnecjointfit"))) {
      c(seconds_of(x$growth), seconds_of(x$survival))
    } else {
      seconds_of(x)
    }
    if (!length(per)) next
    rows[[length(rows) + 1L]] <- data.frame(
      file = paste0("vignettes/fits/", f),
      object = nm,
      n_models = length(per),
      sample_seconds = round(sum(per, na.rm = TRUE), 1)
    )
  }
}

out <- do.call(rbind, rows)
out <- out[order(-out$sample_seconds), ]
write.csv(out, file.path("scripts", "fit_times.csv"), row.names = FALSE)

message(nrow(out), " fits measured, ", sum(out$n_models), " models")
message("slowest: ", out$object[1], " at ", out$sample_seconds[1], "s")
message("written to scripts/fit_times.csv")
