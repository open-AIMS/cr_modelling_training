# ---------------------------------------------------------------------------
# Writes the live scripts, one per taught module, into scripts/live/.
#
# A live script is the plain .R file holding one module's executable code. The
# presenter runs it in Positron beside the projected page, and every
# participant has the identical file open. `notes/delivery-format-human.md`
# decides the format and this script is step 2 of its order of work.
#
# The scripts are generated rather than written by hand, for the reason that
# note gives for rejecting a parallel slide deck: a module and a hand-kept copy
# of its code diverge the first time either is edited, and the revision edits
# every module. Regenerate after any change to a module. The header of each
# generated file says so, and says which module it came from.
#
# Run with:
#   Rscript scripts/generate_live_scripts.R
#
# Takes a second. It reads the .qmd sources and fits nothing.
# ---------------------------------------------------------------------------

source("scripts/qmd_chunks.R")

# The packages the modules attach. They are loaded here so that `exists()` can
# tell a package export from a reference to an object a script never creates.
suppressMessages({
  library(bayesnec)
  library(brms)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

vig_dir <- "vignettes"
out_dir <- file.path("scripts", "live")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SITE <- "https://open-aims.github.io/cr_modelling_training"

# How long a fit may take and still be run in the room.
#
# Set on 2026-09-19 (RF) at about 45 seconds, against a first value of 120. Two
# minutes is the right rule for a single fit and the wrong one for a module: at
# 120 the nine fits of module 5 ran, about 10 minutes of a block the agenda
# gives 30 minutes as a demonstration, and module 5 is about recognising which
# distribution suits an endpoint rather than about executing a call. Below 73
# the seven two-equation sets in module 5 are commented and its two
# single-equation fits still run, while modules 2, 4, 6 and 7 keep all of
# theirs.
#
# The value is 50 rather than 45 because the estimates are not evenly spread.
# The ten fits that run estimate between 36 and 45.1 seconds and the next is
# 73, so anything from 46 to 72 gives the same ten. A threshold of exactly 45
# would turn on the 0.1 second by which module 6's `fit-a` exceeds it, which is
# well inside the error of an estimate built on a compile time measured three
# times.
#
# Raise it to have participants fit more for themselves. The estimate each
# chunk is judged against is printed beside it in the script, so what a change
# would do is visible without regenerating.
MAX_LIVE_SECONDS <- 50

# Stan compilation, per distinct model, in seconds.
#
# This is the part of a fit that the saved object does not record and that
# dominates a small one. Measured on 2026-09-19 on the four-core WSL2 machine:
# a nec3param fit on nec_data took 38.5 s wall against 2.9 s sampling, and an
# ecxll3 fit took 33.5 s against 2.4 s. A second fit of a model already
# compiled in the same session cost 6.6 s, so this applies once per distinct
# equation and not once per call.
#
# A participant's machine is not this machine. Windows with Rtools is generally
# slower, so the figure is a floor and the estimates below are optimistic.
COMPILE_SECONDS <- 35

# What each saved fit sampled for, written by scripts/measure_fit_times.R.
# Absent, every fit call is commented out, which is the safe default.
times_file <- file.path("scripts", "fit_times.csv")
fit_times <- if (file.exists(times_file)) {
  read.csv(times_file, stringsAsFactors = FALSE)
} else {
  message("note: ", times_file, " is absent, so every fit call is commented out.")
  data.frame(file = character(), n_models = integer(), sample_seconds = numeric())
}

# The taught modules, in the order of `notes/workshop-agenda.md`, with the title
# used on the page. `drc-reference` is reference rather than taught and is not
# generated; add it here if that changes.
MODULES <- c(
  "1Software-stack"                              = "The software stack",
  "2Fitting-a-CR-model-using-bayesnec"           = "Fitting a single model",
  "3Toxicity_estimation_and_available_models"    = "Toxicity estimation and the model set",
  "4Model_averaging_and_multimodel_inference"    = "Model averaging and multi-model inference",
  "5Response_data_and_statistical_distributions" = "Response data and statistical distributions",
  "6Priors_and_Bayesian_inference"               = "Priors and Bayesian inference",
  "7Example_case_study"                          = "A worked case study",
  "8Factor_covariates_and_groupings"             = "Factor covariates and groupings"
)

# A section marker Positron and RStudio put in the document outline: a comment
# ending in four or more dashes. The headings are carried through from the
# module so that a question about a named section can be answered by jumping to
# it in the script rather than scrolling.
section <- function(text, level) {
  prefix <- paste0("# ", strrep("  ", level - 2L), text, " ")
  paste0(prefix, strrep("-", max(4L, 78L - nchar(prefix))))
}

# Why a chunk the page does not run is commented out here, if it is.
#
# `eval: false` on the page is not the test. The page withholds a chunk for
# whatever reason suits the render, and most of those chunks are fits, but some
# are a single `options()` call that a participant should run and that returns
# instantly. Commenting those out teaches the wrong thing and wastes the one
# mechanism available for saying "do not run this".
#
# Two reasons stand, and both are about what running the chunk would do in the
# room rather than about how it is marked:
#
#   "samples"    the chunk calls bnec(), bnec_group(), brm() or amend(), so it
#                takes minutes. This is the case the runs-here callout exists
#                for.
#   "undefined"  the chunk refers to objects the script never creates, so
#                running it raises an error. Module 2's `ecx(fit, ...)` is the
#                example: `fit` is a stand-in name, and the script has
#                `bnec_fit`.
#
# A chunk whose body is a bare formula is in the second class. It is written to
# show the shape of a term and its names are placeholders.
#
# Everything else is emitted live.
# bnec_hurdle and bnec_joint are named in full and before `bnec`, because
# `\\bbnec\\s*\\(` does not match `bnec_hurdle(` and module 5's hurdle fit was
# emitted live on 2026-09-21 as a result, at three and a half minutes.
SAMPLING <- "\\b(bnec_hurdle|bnec_joint|bnec_group|bnec|brm|amend)\\s*\\("

# Formulas are removed before the free variables are counted. A formula is not
# evaluated when the call around it is, so the names inside it are not
# references the script has to satisfy: `brms::bf(y ~ bot + ..., nl = TRUE)`
# runs against an empty environment. Counting them would comment out a call
# that takes no measurable time and cannot fail.
strip_formulas <- function(e) {
  if (is.call(e)) {
    # TRUE rather than NULL: assigning NULL into a call removes that element
    # and shrinks the call underneath the loop walking it. A constant is
    # ignored by findGlobals, which is what is wanted here.
    if (identical(e[[1]], as.name("~"))) return(TRUE)
    return(as.call(lapply(as.list(e), strip_formulas)))
  }
  e
}

# The names a block leaves behind, so that a later block referring to them is
# not mistaken for one referring to nothing.
assigned_names <- function(exprs) {
  out <- character()
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    op <- as.character(e[[1]])[1]
    if (op %in% c("<-", "=", "<<-") && length(e) >= 2) {
      lhs <- e[[2]]
      # `x$a <- v` and `x[i] <- v` both still define x.
      while (is.call(lhs)) lhs <- lhs[[2]]
      if (is.name(lhs)) out <<- c(out, as.character(lhs))
    }
    if (op == "for" && length(e) >= 2 && is.name(e[[2]])) {
      out <<- c(out, as.character(e[[2]]))
    }
    # data(nec_data) puts nec_data in the session.
    if (op == "data") {
      for (a in as.list(e)[-1]) if (is.name(a)) out <<- c(out, as.character(a))
    }
    for (i in seq_along(e)) if (!is.null(e[[i]])) walk(e[[i]])
    invisible(NULL)
  }
  for (e in exprs) walk(e)
  unique(out)
}

# The targets of replacement assignments, which are read before they are
# written: `fixed_prior$prior[i] <- v` needs `fixed_prior` to exist already.
# codetools counts one of these as a local assignment and does not report it,
# so it is collected separately.
replacement_targets <- function(exprs) {
  out <- character()
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    if (as.character(e[[1]])[1] %in% c("<-", "=", "<<-") &&
        length(e) >= 2 && is.call(e[[2]])) {
      lhs <- e[[2]]
      while (is.call(lhs)) lhs <- lhs[[2]]
      if (is.name(lhs)) out <<- c(out, as.character(lhs))
    }
    for (i in seq_along(e)) walk(e[[i]])
    invisible(NULL)
  }
  for (e in exprs) walk(e)
  unique(out)
}

# The free variables a block needs and does not create itself.
#
# codetools::findGlobals is used rather than all.vars, because it already
# excludes a function's own arguments, its local assignments and its loop
# variables. all.vars reports every one of those as a reference.
free_vars <- function(exprs) {
  exprs <- lapply(exprs, strip_formulas)
  if (!length(exprs)) return(character())
  f <- eval(call("function", NULL, as.call(c(list(as.name("{")), exprs))))
  tryCatch(codetools::findGlobals(f, merge = FALSE)$variables,
           error = function(e) character())
}

# The number of models an amend() call fits.
#
# amend() adds equations to a set that is already fitted and refits only those,
# so the saved result's model count overstates what the call costs. The models
# added are named in `add =` and can be counted from the source.
amend_added <- function(exprs) {
  n <- 0L
  walk <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    if (identical(e[[1]], as.name("amend"))) {
      a <- tryCatch(eval(e[["add"]]), error = function(err) NULL)
      if (!is.null(a)) n <<- n + length(a)
    }
    for (i in seq_along(e)) walk(e[[i]])
    invisible(NULL)
  }
  for (e in exprs) walk(e)
  n
}

# Seconds a fit chunk would take, or NA where it cannot be estimated.
#
# The sampling time comes from the object the chunk saves. A chunk that fits
# without saving, and module 8's fits, which are made on the cluster and never
# saved here, have no measurement and stay commented.
fit_seconds <- function(ch, exprs) {
  if (!length(ch$saves)) return(NA_real_)
  rows <- fit_times[fit_times$file %in% ch$saves, , drop = FALSE]
  if (nrow(rows) != length(unique(ch$saves))) return(NA_real_)
  n_compile <- if (!is.null(exprs) && amend_added(exprs) > 0) {
    amend_added(exprs)
  } else {
    sum(rows$n_models)
  }
  n_compile * COMPILE_SECONDS + sum(rows$sample_seconds)
}

# A readable duration for the note beside a fit.
duration <- function(secs) {
  if (secs < 90) sprintf("%.0f s", secs) else sprintf("%.0f min", secs / 60)
}

# Where a name the script has not created can be got from a saved object.
#
# Module 3 demonstrates ecnsec() on `bnec_fit` and never loads it, because the
# page reads on from module 2 where it was fitted. Loading it lets the chunk
# run, which serves a participant better than a commented block naming an
# object they have no way to make. The file chosen is the one saved by the
# latest module at or before this one, so module 3 gets module 2's fit rather
# than module 6's, which is a different fit under the same name.
#
# character(0) where any of the names cannot be resolved: a partial repair
# would leave the chunk failing on the rest.
resolve_missing <- function(missing, mod_num, restores) {
  files <- character()
  for (nm in missing) {
    cand <- names(restores)[vapply(restores, function(v) nm %in% v, logical(1))]
    n <- suppressWarnings(as.integer(sub("^vignettes/fits/m(\\d+).*$", "\\1", cand)))
    keep <- !is.na(n) & n <= mod_num
    if (!any(keep)) return(character())
    files <- c(files, cand[keep][which.max(n[keep])])
  }
  unique(files)
}

# Comments out a block, one `#` per line.
#
# `if (FALSE) { ... }` was written first and rejected. It keeps the syntax
# highlighting, but in Positron and RStudio a cursor on a line inside the block
# still sends that line to the console, so a participant working down the file
# line by line starts the fit anyway. That is the failure the runs-here callout
# exists to prevent, and a form of it that is safe to touch is worth more than
# one that looks better projected. The projector shows the rendered page in any
# case, which is where the call is read from.
inert <- function(code) {
  code <- code[seq_len(max(which(nzchar(trimws(code))), 0L))]
  ifelse(nzchar(trimws(code)), paste0("# ", code), "#")
}

# knitr's own configuration is dropped. It does nothing outside a render and
# tells a participant to set something they have no reason to set.
strip_render_only <- function(code) {
  code[!grepl("^\\s*knitr::opts_chunk\\$set", code)]
}

# What each saved file restores, read from the save() call that writes it.
#
# Built across every module rather than per module, because a load() can
# restore an object a different module saved: module 3 loads module 2's
# `bnec_fit`. The save() chunks are ones the page shows and does not run, so
# this is taken from the source and works without the bundle on disk.
restores <- list()
for (stem in names(MODULES)) {
  for (ch in qmd_chunks(file.path(vig_dir, paste0(stem, ".qmd")))) {
    for (ln in grep('save\\s*\\(.*file\\s*=\\s*"vignettes/fits/', ch$code, value = TRUE)) {
      f <- sub('.*file\\s*=\\s*"([^"]+)".*', "\\1", ln)
      nm <- trimws(strsplit(sub('^\\s*save\\s*\\(', "", ln), ",")[[1]][1])
      if (grepl("^[A-Za-z._][A-Za-z0-9._]*$", nm)) {
        restores[[f]] <- unique(c(restores[[f]], nm))
      }
    }
  }
}

for (i in seq_along(MODULES)) {
  stem <- names(MODULES)[i]
  title <- unname(MODULES[i])
  qmd <- file.path(vig_dir, paste0(stem, ".qmd"))
  stopifnot(file.exists(qmd))

  chunks <- qmd_chunks(qmd)
  body <- character()
  n_run <- 0L
  n_inert <- 0L
  n_skip <- 0L
  loads_fits <- FALSE
  reasons <- character()
  n_fit_live <- 0L
  secs_live <- 0
  fitted_live <- character()
  mod_num <- as.integer(substr(stem, 1, 1))

  # The names the script has created by the point each chunk is reached. Only a
  # chunk that is emitted live adds to it: a commented chunk creates nothing.
  defined <- character()

  for (ch in chunks) {
    for (h in ch$headings) {
      level <- nchar(sub("^(#+).*$", "\\1", h))
      body <- c(body, "", section(sub("^#+\\s*", "", h), level), "")
    }

    code <- strip_render_only(ch$code)
    code <- code[seq_len(max(which(nzchar(trimws(code))), 0L))]
    if (!length(code)) next

    # A chunk that is nothing but an image include exists to place a figure on
    # the page. It has no meaning in a console and no later chunk depends on it.
    if (ch$image_only) {
      n_skip <- n_skip + 1L
      next
    }

    label <- if (nzchar(ch$label)) ch$label else sprintf("line %d", ch$line)

    exprs <- tryCatch(as.list(parse(text = code)), error = function(e) NULL)

    # The classification. A chunk the page runs is emitted live whatever it
    # contains, because the page has already run it. A chunk the page withholds
    # is emitted live unless it samples or cannot run.
    reason <- ""
    prelude <- character()
    est <- NA_real_
    if (ch$eval_false) {
      if (is.null(exprs)) {
        reason <- "fragment"
      } else if (any(grepl(SAMPLING, code))) {
        est <- fit_seconds(ch, exprs)
        if (is.na(est) || est > MAX_LIVE_SECONDS) reason <- "samples"
      } else if (all(vapply(exprs, function(e)
                     is.call(e) && identical(e[[1]], as.name("~")), logical(1)))) {
        reason <- "syntax"
      } else {
        missing <- setdiff(union(free_vars(exprs), replacement_targets(exprs)),
                           defined)
        missing <- missing[!vapply(missing, exists, logical(1),
                                   envir = globalenv())]
        if (length(missing)) {
          fixes <- resolve_missing(missing, mod_num, restores)
          if (length(fixes)) {
            prelude <- c(
              sprintf("#> %s comes from an earlier module. The page reads on from there; this line is added so the block runs.",
                      paste(missing, collapse = ", ")),
              sprintf('load("%s")', fixes))
            loads_fits <- TRUE
            defined <- unique(c(defined, unlist(restores[fixes])))
          } else {
            reason <- paste0("undefined:", paste(missing, collapse = ","))
          }
        }
      }
    }

    if (nzchar(reason)) {
      note <- if (reason == "samples" && !is.na(est)) {
        sprintf("#> %s -- not run here: fitting this takes about %s.", label, duration(est))
      } else if (reason == "samples") {
        sprintf("#> %s -- not run here: this samples, and is not among the fits that were timed.", label)
      } else if (reason == "syntax") {
        sprintf("#> %s -- the shape of the term. The names in it are stand-ins.", label)
      } else if (reason == "fragment") {
        sprintf("#> %s -- shown on the page as an illustration and not run.", label)
      } else {
        sprintf("#> %s -- not run here: it refers to %s, which this script does not create.",
                label, gsub(",", ", ", sub("^undefined:", "", reason)))
      }
      tail_note <- if (reason == "samples" && length(ch$saves)) {
        "#> Remove the leading # to fit it yourself; the load() below restores the saved object."
      } else if (reason == "samples") {
        "#> Remove the leading # to fit it yourself, and expect to wait."
      } else {
        "#> Remove the leading # to run it against your own objects."
      }
      body <- c(body, note, tail_note, inert(code), "")
      n_inert <- n_inert + 1L
      reasons <- c(reasons, sub(":.*$", "", reason))
    } else {
      # A chunk hidden on the page is hidden for one of two reasons, and they
      # call for different things from a participant. `include: false` is the
      # setup block, which has to run before anything else and has nothing to
      # look at. `echo: false` draws a figure whose code the page withholds, so
      # the output is on the screen and the code is only here.
      note <- if (ch$setup) {
        sprintf("#> %s -- the libraries and options the rest of the module runs under.", label)
      } else if (!is.na(est)) {
        sprintf("#> %s -- fits the model here, in about %s, most of it compiling.",
                label, duration(est))
      } else if (ch$hidden) {
        sprintf("#> %s -- draws a figure on the page; the page shows the figure and not this code.", label)
      } else {
        sprintf("#> %s", label)
      }
      # The save() inside a fit chunk is commented out where the fit itself
      # runs. It writes into vignettes/fits/, which holds the distributed
      # objects: running the script would replace them with fits made under
      # whatever versions the participant has, and the checksum in
      # vignettes/fits.sha256 that a download is verified against would no
      # longer match. The fit still happens and the object is still in the
      # session; only the overwrite is prevented.
      if (!is.na(est)) {
        is_save <- grepl('save\\s*\\(.*file\\s*=\\s*"vignettes/fits/', code)
        if (any(is_save)) {
          code[is_save] <- paste0("# ", code[is_save],
                                  "   # not run: it would overwrite the distributed fit")
        }
      }

      # A load() of something just fitted above restores the distributed copy.
      # Said plainly, because otherwise it reads as the participant's own fit
      # being thrown away for no reason.
      reload <- character()
      for (f in fitted_live) {
        if (any(grepl(paste0('load\\s*\\(\\s*"', f, '"'), code))) {
          reload <- "#> This restores the distributed object, so what follows matches the page."
        }
      }
      body <- c(body, note, reload, prelude, code, "")
      n_run <- n_run + 1L
      if (!is.na(est)) {
        n_fit_live <- n_fit_live + 1L
        secs_live <- secs_live + est
        fitted_live <- c(fitted_live, ch$saves)
      }
      defined <- unique(c(defined, if (!is.null(exprs)) assigned_names(exprs)))
      for (ln in grep('load\\s*\\(\\s*"vignettes/fits/', code, value = TRUE)) {
        loads_fits <- TRUE
        f <- sub('.*load\\s*\\(\\s*"([^"]+)".*', "\\1", ln)
        defined <- unique(c(defined, restores[[f]]))
      }
    }
  }

  # Every path in a module is written from the project root, as
  # `vignettes/example_binomial.csv`, and a render uses the same working
  # directory (`execute-dir: project` in _quarto.yml). Nothing is rewritten
  # here and no working directory is set: a line in this script is the same
  # text as the line on the page, and both resolve for a participant who has
  # opened cr_modelling_training.Rproj. Until 2026-09-21 the scripts called
  # setwd("vignettes") instead, which left the page and the console needing
  # different paths for the same file.
  preamble <- c(
    'if (!dir.exists("vignettes")) {',
    '  stop("Open cr_modelling_training.Rproj first: every path below is ",',
    '       "written from the project root, as vignettes/...")',
    '}',
    ''
  )

  # No backend is set here. Every module that fits sets
  # options(brms.backend = "cmdstanr") in its own setup chunk, so the line
  # arrives through the module like any other, and a line in this script stays
  # the same text as the line on the page. The generator set it in the preamble
  # for a few hours on 2026-09-19, before the modules were corrected.

  if (loads_fits) {
    preamble <- c(preamble,
      '# The saved model objects. Sampling is what takes the minutes, so the fits',
      '# are made once and loaded here. If this stops the script, fetch them:',
      '#   source("vignettes/fetch_fits.R")',
      'if (!length(list.files("vignettes/fits", pattern = "\\\\.RData$"))) {',
      '  stop("No fitted objects in vignettes/fits/. Run ",',
      '       "source(\\"vignettes/fetch_fits.R\\") first.")',
      '}',
      ''
    )
  }

  header <- c(
    strrep("#", 78),
    paste0("# Module ", substr(stem, 1, 1), " -- ", title),
    "#",
    paste0("# Generated from ", qmd),
    paste0("# by scripts/generate_live_scripts.R on ", format(Sys.Date()), "."),
    "#",
    "# Do not edit this file. Edit the module and run the generator again, or the",
    "# script and the page will disagree.",
    "#",
    paste0("# The page: ", SITE, "/vignettes/", stem, ".html"),
    "#",
    "# Lines beginning `#>` name the chunk each block came from, so a block here",
    "# can be matched to the place on the page it is discussed.",
    "#",
    "# Every block runs unless the line above it says otherwise. A fit is",
    paste0("# commented out where it would take more than ", MAX_LIVE_SECONDS,
           " seconds, and so is a"),
    "# block that refers to objects this script does not create; the line above",
    "# each one says which, and gives the estimate for a fit.",
    "#",
    "# Where a fit does run, its save() is commented out so that the fits you",
    "# downloaded are not overwritten, and the load() after it restores the",
    "# distributed copy so the output below matches the page.",
    strrep("#", 78),
    ""
  )

  path <- file.path(out_dir, paste0(stem, ".R"))
  writeLines(c(header, preamble, body), path)

  # A script that does not parse is a script that fails in the room.
  ok <- tryCatch({ parse(path); TRUE },
                 error = function(e) { message("  PARSE ERROR: ", conditionMessage(e)); FALSE })

  tally <- table(factor(reasons, c("samples", "syntax", "fragment", "undefined")))
  message(sprintf("%-46s %3d run (%d fits, ~%s) | commented: %2d fits, %d syntax, %d undefined%s",
                  basename(path), n_run, n_fit_live,
                  if (n_fit_live) duration(secs_live) else "-",
                  tally[["samples"]],
                  tally[["syntax"]] + tally[["fragment"]], tally[["undefined"]],
                  if (ok) "" else "   [DOES NOT PARSE]"))
}

# The index. It is generated with the scripts so that the module list is
# written once, in MODULES above, rather than kept in step by hand.
readme <- c(
  "# Live scripts",
  "",
  "One `.R` file per taught module, holding the code from that module's page.",
  "Open the file beside the page and run it as the module is worked through.",
  "",
  "| Script | Module |",
  "|---|---|",
  sprintf("| `%s.R` | %s |", names(MODULES), unname(MODULES)),
  "",
  "## Running them",
  "",
  "Open `cr_modelling_training.Rproj` first. Every path in these scripts is",
  "written from the project root, as `vignettes/example_binomial.csv`, which",
  "is the same text the module page shows. Nothing sets the working",
  "directory, so a line copied from the page into the console resolves.",
  "",
  "The model fits are loaded rather than sampled. Fetch them once with",
  "`source(\"vignettes/fetch_fits.R\")`, which takes about half a minute; a",
  "script that loads a fit stops with an instruction if they are absent.",
  "",
  "Every block runs unless the line above it says otherwise. A fit is commented",
  paste0("out where it would take more than ", MAX_LIVE_SECONDS,
         " seconds, and so is a block that refers"),
  "to objects the script does not create. The line above each one says which,",
  "and gives the estimate for a fit. Remove the leading `#` to run a commented",
  "block yourself.",
  "",
  "The fits that do run are the short ones, and most of their time is Stan",
  "compiling the model rather than sampling it. Their `save()` is commented out",
  "so that the fits you downloaded are not overwritten, and the `load()` after",
  "each one restores the distributed copy so that what follows matches the page.",
  "",
  "## Regenerating",
  "",
  "These files are generated from the module sources and are not edited:",
  "",
  "```bash",
  "Rscript scripts/generate_live_scripts.R",
  "```",
  "",
  "Run it after any change to a module, or the script and the page will",
  "disagree while the room is looking at both."
)
writeLines(readme, file.path(out_dir, "README.md"))

message("\n", length(MODULES), " scripts and a README in ", out_dir)
