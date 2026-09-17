# ---------------------------------------------------------------------------
# Downloads the fitted models the workshop modules use, and extracts them into
# `vignettes/fits/`.
#
# Run this once, before the workshop, from the cr_modelling_training folder:
#
#   source("vignettes/fetch_fits.R")
#
# Fitting the models yourself works and is slow: modules 5 and 6 take about
# thirteen and sixteen minutes each on four cores. The archive is about 39 MB.
# With these objects in place the same modules run in seconds, because every
# fit call finds a saved result.
#
# Nothing here is required. A participant without the archive fits the models
# instead, and the material behaves identically.
# ---------------------------------------------------------------------------

# Where the archive is hosted. A shared folder is used while the material is
# still being revised, because the file can be replaced in place and the link
# stays the same. A Zenodo deposit is the destination once the content is
# frozen, since it gives a permanent identifier, and the only change needed
# then is the address below.
#
# A OneDrive or SharePoint share link serves a web page rather than the file, so
# the address below ends in `&download=1`. Test any replacement by running this
# script rather than by opening the link in a browser: measured 2026-09-17, the
# plain share link returns 200 with 58 KB of HTML, which lands in a file named
# .zip and fails later at unzip rather than here.
#
# CR_FITS_URL overrides this for testing against any host.
FITS_URL <- paste0(
  "https://aimsgovau-my.sharepoint.com/:u:/g/personal/r_fisher_aims_gov_au/",
  "IQDag6rnW-2cQqI2C8Dozb9jARrz7DFsmv22c_Kvq5le-WQ?e=ucW2LI&download=1")
ARCHIVE <- "cr_modelling_fits.zip"

# The number of saved fits the archive holds. Checked after extraction, because
# a truncated download has a plausible size and unpacks into a subset that
# fails later and obscurely.
#
# Count this against vignettes/fits/ when the archive is rebuilt rather than
# trusting the number below. Modules are being revised in parallel and each
# revision that changes a fitted model adds or replaces objects here; the value
# was 18, then 20, and was 25 on 2026-09-17. An under-set value does not break
# the fetch, it only guards less.
N_FITS <- 25L

fits_url <- function() {
  override <- Sys.getenv("CR_FITS_URL", unset = "")
  if (nzchar(override)) return(override)
  if (is.na(FITS_URL)) {
    stop("The download location is not set yet.\n",
         "  Edit FITS_URL at the top of vignettes/fetch_fits.R, or set\n",
         "  CR_FITS_URL to the archive's address. Until then, fit the models\n",
         "  as the modules do; everything works, it is only slower.",
         call. = FALSE)
  }
  FITS_URL
}

fetch_fits <- function(dest = "vignettes/fits", quiet = FALSE) {
  if (!dir.exists("vignettes")) {
    stop("Run this from the cr_modelling_training folder, the one holding ",
         "'vignettes'.", call. = FALSE)
  }

  url <- fits_url()
  tmp <- tempfile(fileext = ".zip")
  if (!quiet) message("Downloading the fitted models, about 39 MB. This is ",
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

  # The checksum is advisory rather than fatal. The archive is replaced while
  # the material is being revised, so a clone taken before the last upload
  # holds an older checksum and would otherwise refuse a download that is
  # perfectly good. Whether the archive unpacks into the expected objects is
  # the test that catches a truncated file, and it is applied below.
  sha_file <- file.path("vignettes", "fits.sha256")
  if (file.exists(sha_file) && requireNamespace("digest", quietly = TRUE)) {
    expected <- sub("\\s.*$", "", readLines(sha_file, warn = FALSE)[1])
    got <- digest::digest(tmp, algo = "sha256", file = TRUE)
    if (!identical(got, expected) && !quiet) {
      message("The archive does not match the checksum recorded in this copy ",
              "of the repository.\n",
              "  That is expected where the fits have been updated since you ",
              "cloned it. Run `git pull` to take the current checksum if you ",
              "want the check to pass.")
    }
  }

  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  unpacked <- tryCatch(utils::unzip(tmp, exdir = "."),
                       error = function(e) character(0),
                       warning = function(w) character(0))
  unlink(tmp)

  n <- length(list.files(dest, pattern = "\\.RData$"))
  if (!length(unpacked) || n < N_FITS) {
    stop("The archive is incomplete: ", n, " of ", N_FITS, " fitted models ",
         "were extracted, so the download was truncated or corrupt.\n",
         "  Delete nothing and simply run this again.", call. = FALSE)
  }

  if (!quiet) message("Ready: ", n, " fitted models in ", normalizePath(dest))
  invisible(TRUE)
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  fetch_fits()
}
