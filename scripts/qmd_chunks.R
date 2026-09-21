# ---------------------------------------------------------------------------
# Chunk extraction from a Quarto module.
#
# Two scripts read a module's code. `generate_taught_fits.R` runs the chunks
# that fit and save the models; `generate_live_scripts.R` writes the .R file
# the presenter runs in the room. Both must see the same chunks in the same
# order, so the parsing lives here rather than being written twice. A live
# script that has drifted from the page it came from is worse than no live
# script, because the room is looking at both at once.
#
# knitr::purl() is not used, and this is the reason the parsing is done by
# hand. Chunk option hooks are applied when a document is knitted and not when
# it is purled, so purling extracts the `eval: false` demonstrations with
# nothing to mark them as code that is not meant to run -- module 4's parallel
# plan among them.
# ---------------------------------------------------------------------------

# Returns one record per `{r}` chunk, in document order, with the chunk options
# already interpreted. The callers filter; nothing is dropped here.
#
# Fields:
#   line        the line the opening fence sits on
#   label       the chunk label, or "" where the chunk has none
#   headings    the section headings between the previous chunk and this one
#   code        the chunk body with the `#|` option lines removed
#   eval_false  the chunk is displayed on the page and not run
#   hidden      the chunk runs on the page but its code is not shown
#   setup       `include: false`: neither the code nor its output appears
#   image_only  the body is nothing but knitr::include_graphics()
#   saves       the vignettes/fits/ paths the chunk writes, if any
qmd_chunks <- function(qmd) {
  lines <- readLines(qmd, warn = FALSE)
  starts <- grep("^```\\{r[ ,}]", lines)
  ends <- grep("^```\\s*$", lines)
  # The title is a single #, so only ## and below are section headings.
  head_at <- grep("^#{2,4} \\S", lines)

  out <- list()
  prev <- 0L
  for (s in starts) {
    e <- ends[ends > s][1]
    if (is.na(e)) next
    body <- if (e - s > 1) lines[(s + 1):(e - 1)] else character()
    opts <- grep("^#\\|", body, value = TRUE)
    code <- body[!grepl("^#\\|", body)]
    label <- grep("^#\\|\\s*label:", opts, value = TRUE)

    here <- head_at[head_at > prev & head_at < s]
    body_code <- code[nzchar(trimws(code))]

    out[[length(out) + 1L]] <- list(
      line = s,
      label = if (length(label)) trimws(sub("^#\\|\\s*label:\\s*", "", label[1])) else "",
      headings = lines[here],
      code = code,
      eval_false = any(grepl("^#\\|\\s*eval:\\s*false", opts)),
      hidden = any(grepl("^#\\|\\s*(echo|include):\\s*false", opts)),
      setup = any(grepl("^#\\|\\s*include:\\s*false", opts)),
      image_only = length(body_code) > 0 &&
        all(grepl("include_graphics", body_code)),
      saves = sub('.*file\\s*=\\s*"([^"]+)".*', "\\1",
                  grep('save\\s*\\(.*file\\s*=\\s*"vignettes/fits/', code, value = TRUE))
    )
    prev <- e
  }
  out
}

# The code `generate_taught_fits.R` evaluates, as one character vector.
#
# An executed chunk is included, because later chunks depend on the data it
# reads and the objects it builds. A chunk marked `eval: false` is skipped,
# except where it saves into `vignettes/fits/`: those are the fit chunks, and
# they are the whole point of that script. Skipping every `eval: false` chunk,
# which an earlier version did, left it unable to fit anything at all.
runnable_chunks <- function(qmd) {
  chunks <- qmd_chunks(qmd)
  out <- character()
  seen <- character()
  for (ch in chunks) {
    if (ch$eval_false && length(ch$saves) == 0) next
    # A module can show the same fit call twice, once to fit and once to
    # illustrate the fit-save-load pattern, which would fit the object twice.
    # A chunk whose every target has already been written this run is skipped.
    if (length(ch$saves) > 0 && all(ch$saves %in% seen)) next
    seen <- c(seen, ch$saves)
    out <- c(out, ch$code, "")
  }
  out
}
