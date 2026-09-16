# ---------------------------------------------------------------------------
# Packages the fitted models in `fits/` into one archive for distribution, with
# a checksum and a manifest.
#
# The archive is what participants download before the workshop, and what goes
# on a USB stick for the day. It is not committed: the objects are large and the
# repository is public.
#
# Run with:
#   Rscript scripts/bundle_fits.R
#
# Produces, in `dist/`:
#   cr_modelling_fits.zip       the archive
#   cr_modelling_fits.sha256    its checksum
#
# and writes vignettes/fits.sha256, which is committed and is what
# fetch_fits.R verifies a download against.
#   cr_modelling_fits.txt       what is inside it and what built it
# ---------------------------------------------------------------------------

fits_dir <- "vignettes/fits"
dist_dir <- "dist"
stem <- "cr_modelling_fits"

if (!dir.exists(fits_dir) || length(list.files(fits_dir, pattern = "\\.RData$")) == 0) {
  stop("No saved fits in '", fits_dir, "'. Render the site, or run ",
       "scripts/generate_taught_fits.R, before bundling.", call. = FALSE)
}

dir.create(dist_dir, showWarnings = FALSE, recursive = TRUE)

rds <- list.files(fits_dir, pattern = "\\.RData$", full.names = TRUE)
sizes <- file.size(rds)
message("bundling ", length(rds), " objects, ",
        round(sum(sizes) / 1024^2, 1), " MB uncompressed")

zip_path <- file.path(dist_dir, paste0(stem, ".zip"))
if (file.exists(zip_path)) unlink(zip_path)

# Paths are stored relative to the repository root, so the archive extracts to
# `fits/` wherever the participant unpacks it. This runs from the repository
# root, so no directory change is needed; an earlier version changed into the
# parent of `fits/` and then wrote the archive to `../dist/`, which is outside
# the repository and does not exist.
utils::zip(zipfile = zip_path,
           files = file.path(fits_dir, basename(rds)),
           flags = "-q9X")

# A checksum rather than a size, because a truncated download is the failure
# this guards against and a truncated file has a plausible size.
sha <- if (requireNamespace("digest", quietly = TRUE)) {
  digest::digest(zip_path, algo = "sha256", file = TRUE)
} else {
  out <- system2("sha256sum", shQuote(zip_path), stdout = TRUE)
  sub("\\s.*$", "", out)
}
writeLines(paste(sha, basename(zip_path)),
           file.path(dist_dir, paste0(stem, ".sha256")))

# The same checksum is written into the repository, because that is where
# fetch_fits.R reads it from: a participant who cloned the repository has this
# file and can therefore verify a download against it. `dist/` is not committed,
# so a checksum kept only there could never be checked.
writeLines(paste(sha, basename(zip_path)), file.path("vignettes", "fits.sha256"))

writeLines(c(
  "Fitted models for the concentration-response modelling workshop.",
  "",
  paste("built:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("objects:", length(rds)),
  paste("archive MB:", round(file.size(zip_path) / 1024^2, 1)),
  paste("sha256:", sha),
  paste("R:", getRversion()),
  paste("bayesnec:", as.character(packageVersion("bayesnec"))),
  paste("brms:", as.character(packageVersion("brms"))),
  "",
  "Extract into the root of the cr_modelling_training folder, so that the",
  "objects end up in vignettes/fits/. vignettes/fetch_fits.R does this for you.",
  "",
  "Each file is named after a hash of the fit call that produced it. A module",
  "whose fit call has changed since this was built simply refits that one",
  "model; nothing silently returns a stale object."
), file.path(dist_dir, paste0(stem, ".txt")))

message("wrote ", zip_path, " (", round(file.size(zip_path) / 1024^2, 1), " MB)")
message("sha256 ", sha)
