# ---------------------------------------------------------------------------
# Downloads the fitted models the workshop modules use, and extracts them into
# `vignettes/fits/`.
#
# Run this once, before the workshop, from the cr_modelling_training folder:
#
#   source("vignettes/fetch_fits.R")
#
# Fitting the models yourself works and is slow: modules 5 and 6 take about
# thirteen and sixteen minutes each on four cores. The archive is about 42 MB.
# With these objects in place the same modules run in seconds, because every
# fit call finds a saved result.
#
# Nothing here is required. A participant without the archive fits the models
# instead, and the material behaves identically.
# ---------------------------------------------------------------------------

# Where the archive is hosted. A release asset on the course repository, under
# a tag that does not change, so the address below stays valid when the archive
# is rebuilt: `gh release upload fits dist/cr_modelling_fits.zip --clobber`
# replaces the file and the link still resolves.
#
# It is fetched anonymously. That is the property being bought here, and it is
# what a share link on institutional OneDrive or SharePoint does not give:
# measured 2026-09-17, the AIMS share link used before this one redirected to
# the file and then returned 403 to a client that was not signed in to the
# tenant, so no participant outside AIMS could have run this script.
#
# Test any replacement by running this script rather than by opening the link
# in a browser. A SharePoint link without `&download=1` returns 200 with 58 KB
# of HTML, which lands in a file named .zip and fails later at unzip.
#
# CR_FITS_URL overrides this for testing against any host.
FITS_URL <- paste0(
  "https://github.com/open-AIMS/cr_modelling_training/",
  "releases/download/fits/cr_modelling_fits.zip")
ARCHIVE <- "cr_modelling_fits.zip"

# The number of saved fits the archive holds. Checked after extraction, because
# a truncated download has a plausible size and unpacks into a subset that
# fails later and obscurely.
#
# Count this against vignettes/fits/ when the archive is rebuilt rather than
# trusting the number below. Modules are being revised in parallel and each
# revision that changes a fitted model adds or replaces objects here; the value
# was 18, then 20, 25 on 2026-09-17 and 26 on 2026-09-19. On 2026-09-21 four
# objects were dropped, because module 2's `ecxll3` fit moved to module 6 under
# a new name and module 6's prior sections for a model set were cut, and two
# were added, `m5_exp_6` for the hurdle example and `m6_ecxll3_fit`. On
# 2026-09-22 module 5's constructed examples were replaced with published
# datasets, which merged its two dispersion fits into one. An under-set value
# does not break the fetch, it only guards less.
N_FITS <- 23L

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
  if (!quiet) message("Downloading the fitted models, about 42 MB. This is ",
                      "done once.")

  # R's default download timeout is 60 seconds, which is not enough for an
  # archive this size on anything but a fast connection. Measured 2026-09-21 on
  # a domestic connection: the default stopped at 9.7 MB of 42.2 MB with
  # "Timeout of 60 seconds was reached", and the script then reported the
  # download as failed.
  #
  # The value is an hour rather than something tidier because the timeout is a
  # total transfer limit, not an idle limit, so it has to cover the whole
  # download at whatever rate the connection gives. Measured 2026-09-22 on the
  # same machine, the link ran at 15 KB/s for a sustained period, which puts a
  # 39 MB archive at about 45 minutes; a 30 minute limit set the day before
  # would have failed on a download that was progressing normally. Conference
  # wifi on the morning is the case this exists for.
  #
  # The timeout is raised for this call and restored afterwards, so a
  # participant's own session setting is left as they had it.
  old_timeout <- getOption("timeout")
  options(timeout = max(old_timeout, 3600))
  on.exit(options(timeout = old_timeout), add = TRUE)

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
