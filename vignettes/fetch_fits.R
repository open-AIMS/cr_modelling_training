# ---------------------------------------------------------------------------
# Downloads the fitted models the workshop modules use, and extracts them into
# `vignettes/fits/`.
#
# Run this once, before the workshop, from the cr_modelling_training folder:
#
#   source("vignettes/fetch_fits.R")
#
# Fitting the models yourself works and is slow: modules 5 and 6 take about
# thirteen and sixteen minutes each on four cores. The archive is about 23 MB. With these objects in place
# the same modules run in seconds, because every fit call finds a saved result.
#
# Nothing here is required. A participant without the archive fits the models
# instead, and the material behaves identically.
# ---------------------------------------------------------------------------

# The archive is deposited on Zenodo, which gives it a permanent identifier and
# needs no account to download. Set ZENODO_RECORD to the numeric record id once
# the deposit exists; CR_FITS_URL overrides it for testing against any host.
ZENODO_RECORD <- NA_character_
ARCHIVE <- "cr_modelling_fits.zip"

fits_url <- function() {
  override <- Sys.getenv("CR_FITS_URL", unset = "")
  if (nzchar(override)) return(override)
  if (is.na(ZENODO_RECORD)) {
    stop("The download location is not set yet.\n",
         "  Edit ZENODO_RECORD at the top of vignettes/fetch_fits.R, or set\n",
         "  CR_FITS_URL to the archive's address. Until then, fit the models\n",
         "  as the modules do; everything works, it is only slower.",
         call. = FALSE)
  }
  sprintf("https://zenodo.org/records/%s/files/%s?download=1",
          ZENODO_RECORD, ARCHIVE)
}

fetch_fits <- function(dest = "vignettes/fits", quiet = FALSE) {
  if (!dir.exists("vignettes")) {
    stop("Run this from the cr_modelling_training folder, the one holding ",
         "'vignettes'.", call. = FALSE)
  }

  url <- fits_url()
  tmp <- tempfile(fileext = ".zip")
  if (!quiet) message("Downloading the fitted models, about 23 MB. This is ",
                      "done once.")

  ok <- tryCatch({
    utils::download.file(url, tmp, mode = "wb", quiet = quiet)
    TRUE
  }, error = function(e) {
    message("The download failed: ", conditionMessage(e))
    message("Fit the models instead by working through the modules as written. ",
            "Nothing is lost; it takes longer.")
    FALSE
  })
  if (!ok) return(invisible(FALSE))

  # The checksum is the point of downloading one at all: a truncated archive has
  # a plausible size and unpacks into objects that fail later and obscurely.
  sha_file <- file.path("vignettes", "fits.sha256")
  if (file.exists(sha_file)) {
    expected <- sub("\\s.*$", "", readLines(sha_file, warn = FALSE)[1])
    got <- if (requireNamespace("digest", quietly = TRUE)) {
      digest::digest(tmp, algo = "sha256", file = TRUE)
    } else NA_character_
    if (!is.na(got) && !identical(got, expected)) {
      unlink(tmp)
      stop("The downloaded archive does not match its checksum, so it is ",
           "incomplete or corrupt. Delete nothing and simply run this again.",
           call. = FALSE)
    }
  }

  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  utils::unzip(tmp, exdir = ".")
  unlink(tmp)

  n <- length(list.files(dest, pattern = "\\.RData$"))
  if (!quiet) message("Ready: ", n, " fitted models in ", normalizePath(dest))
  invisible(TRUE)
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  fetch_fits()
}
