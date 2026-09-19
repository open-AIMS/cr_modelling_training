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

vig_dir <- "vignettes"
out_dir <- file.path("scripts", "live")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SITE <- "https://open-aims.github.io/cr_modelling_training"

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

# Comments out a block that the page shows and does not run.
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

    if (ch$eval_false) {
      note <- if (length(ch$saves)) {
        sprintf("#> %s -- the fit. Shown on the page and not run: it samples for minutes.", label)
      } else {
        sprintf("#> %s -- shown on the page and not run.", label)
      }
      body <- c(body, note,
                "#> Remove the leading # from the lines below to run it yourself.",
                inert(code), "")
      n_inert <- n_inert + 1L
    } else {
      # A chunk hidden on the page is hidden for one of two reasons, and they
      # call for different things from a participant. `include: false` is the
      # setup block, which has to run before anything else and has nothing to
      # look at. `echo: false` draws a figure whose code the page withholds, so
      # the output is on the screen and the code is only here.
      note <- if (ch$setup) {
        sprintf("#> %s -- the libraries and options the rest of the module runs under.", label)
      } else if (ch$hidden) {
        sprintf("#> %s -- draws a figure on the page; the page shows the figure and not this code.", label)
      } else {
        sprintf("#> %s", label)
      }
      body <- c(body, note, code, "")
      n_run <- n_run + 1L
      if (any(grepl('load\\s*\\(\\s*"fits/', code))) loads_fits <- TRUE
    }
  }

  # Every path in a module resolves from vignettes/, which is where the module
  # is rendered from and where fits/ and the CSV files sit. The working
  # directory is set rather than the paths rewritten, so that a line in this
  # script is the same text as the line on the page. A participant comparing
  # the two should find no difference to explain.
  preamble <- c(
    'if (basename(getwd()) != "vignettes") {',
    '  if (dir.exists("vignettes")) {',
    '    setwd("vignettes")',
    '  } else {',
    '    stop("Open the cr_modelling_training project first: every path below is ",',
    '         "relative to its vignettes/ folder.")',
    '  }',
    '}',
    ''
  )

  if (loads_fits) {
    preamble <- c(preamble,
      '# The saved model objects. Sampling is what takes the minutes, so the fits',
      '# are made once and loaded here. If this stops the script, fetch them:',
      '#   source("fetch_fits.R")',
      'if (!length(list.files("fits", pattern = "\\\\.RData$"))) {',
      '  stop("No fitted objects in vignettes/fits/. Run source(\\"fetch_fits.R\\") first.")',
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
    "# A block that is commented out is one the page shows and does not run. Those",
    "# are the model fits and the calls that are there to be read rather than",
    "# executed. Below each fit is the load() that restores the saved object.",
    strrep("#", 78),
    ""
  )

  path <- file.path(out_dir, paste0(stem, ".R"))
  writeLines(c(header, preamble, body), path)

  # A script that does not parse is a script that fails in the room.
  ok <- tryCatch({ parse(path); TRUE },
                 error = function(e) { message("  PARSE ERROR: ", conditionMessage(e)); FALSE })

  message(sprintf("%-46s %3d run, %2d shown only, %2d image chunks dropped%s",
                  basename(path), n_run, n_inert, n_skip,
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
  "Open `cr_modelling_training.Rproj` first. Each script sets the working",
  "directory to `vignettes/`, because that is where the data files and the",
  "saved model objects sit, and every path in the module is written relative",
  "to it.",
  "",
  "The model fits are loaded rather than sampled. Fetch them once with",
  "`source(\"vignettes/fetch_fits.R\")`, which takes about half a minute; a",
  "script that loads a fit stops with an instruction if they are absent.",
  "",
  "A block that is commented out is one the page shows and does not run. Those",
  "are the `bnec()` calls, which sample for minutes, and the calls written to",
  "be read rather than executed. Remove the leading `#` to run one yourself.",
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
