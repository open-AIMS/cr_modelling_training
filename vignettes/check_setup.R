# ---------------------------------------------------------------------------
# Pre-workshop setup check
# Concentration-response modelling course, SETAC Australasia 2026
#
# Run this in a fresh R session AFTER following the setup instructions.
# It installs nothing. It reports what you have and whether it works.
#
# Send the whole output back before the workshop. There is no time to
# diagnose installations on the day.
# ---------------------------------------------------------------------------

# bayesnec is installed from its development branch, so there is no fixed
# version to check against: dev moves, and pinning a number here would report a
# mismatch for everyone the day after it was written. The version is recorded
# instead, so that a result can be attributed to a state of the package.

# Stages run in increasing order of expense, so a failure is localised to the
# cheapest stage that shows it. Compilation is stage 4 and sampling stage 5;
# everything before them is inspection only.
results <- new.env(parent = emptyenv())
record <- function(stage, ok, detail = "") {
  assign(stage, list(ok = ok, detail = detail), envir = results)
  cat(sprintf("[%s] %-28s %s\n",
              if (isTRUE(ok)) "PASS" else if (is.na(ok)) "WARN" else "FAIL",
              stage, detail))
}

rule <- function(x) cat("\n", x, "\n", strrep("-", nchar(x)), "\n", sep = "")

# --- Stage 1: R and platform --------------------------------------------------
rule("Stage 1: R and platform")

cat("R version : ", R.version.string, "\n", sep = "")
cat("Platform  : ", R.version$platform, "\n", sep = "")
cat("OS        : ", utils::sessionInfo()$running, "\n", sep = "")
cat("HOME      : ", Sys.getenv("HOME"), "\n", sep = "")
cat("Library   : ", .libPaths()[1], "\n", sep = "")

r_ok <- getRversion() >= "4.4.0"
record("r_version", r_ok,
       if (r_ok) as.character(getRversion())
       else paste0(getRversion(), " is too old; install R 4.4.0 or later"))

# HOME is checked because the CmdStan build fails on paths that are synced by
# OneDrive, contain a space, or contain non-ASCII characters, and the default
# install location is ~/.cmdstan. The failure appears as a compiler error, so
# it is worth naming the cause here rather than leaving it to be diagnosed
# from build output.
home <- Sys.getenv("HOME")
home_problems <- c(
  if (grepl("OneDrive", home, ignore.case = TRUE)) "is inside OneDrive",
  if (grepl(" ", home)) "contains a space",
  if (grepl("[^\x01-\x7f]", home)) "contains non-ASCII characters"
)
record("home_path",
       if (length(home_problems)) NA else TRUE,
       if (length(home_problems))
         paste0("HOME ", paste(home_problems, collapse = " and "),
                "; install CmdStan to a plain path, e.g. ",
                'install_cmdstan(dir = "C:/cmdstan")')
       else "plain path")

# --- Stage 2: packages --------------------------------------------------------
rule("Stage 2: packages")

# rstan must be installed because brms imports it, but it does not have to be
# the sampling backend. The course uses cmdstanr, which is the more reliable
# route on Windows.
needed <- c("bayesnec", "brms", "cmdstanr", "rstan", "posterior", "ggplot2", "dplyr")
versions <- vapply(needed, function(p) {
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
}, character(1))

for (p in needed) {
  cat(sprintf("  %-10s %s\n", p, if (is.na(versions[[p]])) "NOT INSTALLED" else versions[[p]]))
}

missing <- needed[is.na(versions)]
record("packages", length(missing) == 0,
       if (length(missing)) paste("missing:", paste(missing, collapse = ", "))
       else "all present")

bn <- versions[["bayesnec"]]
record("bayesnec_version", !is.na(bn),
       if (is.na(bn))
         "not installed; see the setup instructions for the dev install"
       else paste0(bn, " (from the dev branch)"))

# --- Stage 3: C++ toolchain and CmdStan ---------------------------------------
rule("Stage 3: toolchain and CmdStan")

if (is.na(versions[["cmdstanr"]])) {
  record("toolchain", FALSE, "cmdstanr not installed, cannot check")
  record("cmdstan", FALSE, "cmdstanr not installed, cannot check")
} else {
  tc <- tryCatch({ cmdstanr::check_cmdstan_toolchain(quiet = TRUE); TRUE },
                 error = function(e) conditionMessage(e))
  record("toolchain", isTRUE(tc), if (isTRUE(tc)) "ok" else as.character(tc))

  cs <- tryCatch(as.character(cmdstanr::cmdstan_version()),
                 error = function(e) NA_character_)
  record("cmdstan", !is.na(cs),
         if (!is.na(cs)) paste0(cs, " at ", cmdstanr::cmdstan_path())
         else "not found; run cmdstanr::install_cmdstan(overwrite = TRUE, quiet = FALSE)")

  # A check that the CmdStan TBB directory is on PATH was added here and then
  # removed. A compiled Stan program does link against tbb.dll, but CmdStanRun
  # puts that directory on the PATH of the chain process itself, so the entry is
  # not required in the user's environment. Measured on Windows 2026-09-11:
  # sampling succeeded with the directory removed from PATH. The check reported
  # a warning on a correctly installed machine, which is worse than reporting
  # nothing.

  # cmdstanr decides which Rtools to look for from the R minor version, and the
  # mapping changed. 0.9.0 sends R 4.5.x to RTOOLS45_HOME; 0.8.0 sent every R
  # from 4.4 onward to RTOOLS44_HOME. A participant following the setup
  # instructions installs Rtools45, so an older cmdstanr searches for a toolchain
  # they were never told to install and reports it as missing.
  if (!is.na(versions[["cmdstanr"]]) && .Platform$OS.type == "windows") {
    cmdstanr_ok <- utils::packageVersion("cmdstanr") >= "0.9.0"
    record("cmdstanr_version", if (cmdstanr_ok) TRUE else NA,
           if (cmdstanr_ok) as.character(utils::packageVersion("cmdstanr"))
           else paste0(utils::packageVersion("cmdstanr"),
                       " expects Rtools44; install 0.9.0 or later for Rtools45"))
  }
}

# --- Stage 4: compile a Stan model --------------------------------------------
rule("Stage 4: compile (this takes a minute or two)")

# The bundled bernoulli example is used rather than a brms model because it
# isolates the toolchain: if this compiles, any later failure is in brms or
# bayesnec rather than in the C++ setup.
compiled <- NULL
if (isTRUE(get0("cmdstan", envir = results, ifnotfound = list(ok = FALSE))$ok)) {
  stan_file <- file.path(cmdstanr::cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")

  # force_recompile = TRUE because cmdstan_model() otherwise reuses any
  # executable newer than the .stan file, including one built months earlier by
  # a toolchain since removed. Without it this stage reports PASS without
  # invoking the compiler, which is the one thing it exists to test. The price
  # is a compile on every run rather than only the first.
  #
  # A successful build emits several hundred lines of warnings from CmdStan's
  # bundled TBB. They are captured rather than printed, because the participant
  # is asked to return the whole output and it has to stay readable. On a
  # failure they are the diagnosis, so the tail is printed then.
  build_log <- utils::capture.output(
    compiled <- tryCatch(
      cmdstanr::cmdstan_model(stan_file, force_recompile = TRUE, quiet = TRUE),
      error = function(e) conditionMessage(e)),
    type = "message")

  compile_ok <- inherits(compiled, "CmdStanModel")
  record("compile", compile_ok,
         if (compile_ok) "bernoulli example compiled" else as.character(compiled))

  if (!compile_ok && length(build_log)) {
    cat("\n  last 40 lines of compiler output:\n")
    cat(paste0("  | ", utils::tail(build_log, 40)), sep = "\n")
    cat("\n")
  }
} else {
  record("compile", FALSE, "skipped, CmdStan not available")
}

# --- Stage 5: sample ----------------------------------------------------------
rule("Stage 5: sample")

if (inherits(compiled, "CmdStanModel")) {
  fit <- tryCatch(
    compiled$sample(data = list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1)),
                    chains = 1, iter_warmup = 200, iter_sampling = 200,
                    refresh = 0, show_messages = FALSE),
    error = function(e) conditionMessage(e))
  # $sample() does not throw when every chain fails to start: it warns and
  # returns a CmdStanMCMC whose draws cannot be read. Testing with inherits()
  # alone therefore records PASS and then halts the script at the next line,
  # where the summary is read -- on exactly the machines this script exists to
  # diagnose. The result is read inside its own tryCatch so the failure is
  # recorded and the run still reaches the summary.
  theta <- if (inherits(fit, "CmdStanMCMC")) {
    tryCatch(fit$summary("theta")$mean, error = function(e) conditionMessage(e))
  } else NULL

  ok <- is.numeric(theta) && length(theta) == 1L && is.finite(theta)
  record("sample", ok,
         if (ok) sprintf("theta mean %.3f", theta)
         else if (is.character(theta)) theta
         else as.character(fit))

  # A CmdStanMCMC whose chains failed emits "No CmdStan config files found" from
  # its finaliser when the session ends, placing an error after the summary and
  # after the NOT READY line. Dropping the reference keeps it out of the output
  # the participant is asked to return.
  if (!ok) {
    rm(fit)
    invisible(gc())
  }
} else {
  record("sample", FALSE, "skipped, model did not compile")
}

# --- Stage 6: brms handoff ----------------------------------------------------
rule("Stage 6: brms with the cmdstanr backend")

# This is the step that matters for the course: bayesnec calls brms, and brms
# has to hand the model to CmdStan. A toolchain that compiles the bernoulli
# example can still fail here if the backend option is not set.
if (!is.na(versions[["brms"]]) && isTRUE(get0("sample", envir = results,
                                              ifnotfound = list(ok = FALSE))$ok)) {
  options(brms.backend = "cmdstanr")
  set.seed(42)
  d <- data.frame(x = rep(1:5, 4), y = rnorm(20))
  bfit <- tryCatch(
    brms::brm(y ~ x, data = d, chains = 1, iter = 400, refresh = 0, silent = 2),
    error = function(e) conditionMessage(e))
  ok <- inherits(bfit, "brmsfit")
  record("brms", ok, if (ok) "fitted via cmdstanr" else as.character(bfit))
} else {
  record("brms", FALSE, "skipped, earlier stage failed")
}

# --- Summary ------------------------------------------------------------------
rule("Summary")

all_stages <- ls(results)
status <- vapply(all_stages, function(s) get(s, envir = results)$ok, logical(1))
n_fail <- sum(!status & !is.na(status))
n_warn <- sum(is.na(status))

if (n_fail == 0 && n_warn == 0) {
  cat("READY. Every check passed.\n")
} else if (n_fail == 0) {
  cat("READY, with ", n_warn, " warning(s) above. Read them, but you can follow along.\n", sep = "")
} else {
  cat("NOT READY. ", n_fail, " check(s) failed:\n", sep = "")
  for (s in all_stages[!status & !is.na(status)]) {
    cat("  - ", s, ": ", get(s, envir = results)$detail, "\n", sep = "")
  }
  cat("\nSend this entire output back before the workshop.\n")
}

cat("\nGenerated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
