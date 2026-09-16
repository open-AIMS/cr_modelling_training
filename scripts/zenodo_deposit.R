# ---------------------------------------------------------------------------
# Uploads the fitted-model archive to Zenodo as a draft deposition.
#
# It stops short of publishing. A published Zenodo record cannot be deleted,
# only a new version issued, so the last step is left to a person looking at the
# draft in a browser.
#
# Needs a personal access token with the `deposit:write` and `deposit:actions`
# scopes, from https://zenodo.org/account/settings/applications/tokens/new/:
#
#   export ZENODO_TOKEN=...
#   Rscript scripts/zenodo_deposit.R
#
# Add --sandbox to work against https://sandbox.zenodo.org, which is a separate
# service with separate tokens and is the right place to rehearse this.
#
# On success it prints the record id. Put that in ZENODO_RECORD at the top of
# vignettes/fetch_fits.R and the download route is live.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
sandbox <- "--sandbox" %in% args
base <- if (sandbox) "https://sandbox.zenodo.org" else "https://zenodo.org"

token <- Sys.getenv("ZENODO_TOKEN", unset = "")
if (!nzchar(token)) {
  stop("Set ZENODO_TOKEN. See the header of this file.", call. = FALSE)
}
for (pkg in c("httr", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is needed and is not installed.", call. = FALSE)
  }
}

zip_path <- "dist/cr_modelling_fits.zip"
txt_path <- "dist/cr_modelling_fits.txt"
if (!file.exists(zip_path)) {
  stop("No archive at ", zip_path, ". Run scripts/bundle_fits.R first.",
       call. = FALSE)
}

api <- function(method, path, ...) {
  r <- httr::VERB(method, paste0(base, path),
                  httr::add_headers(Authorization = paste("Bearer", token)),
                  ...)
  if (httr::status_code(r) >= 300) {
    stop("Zenodo returned ", httr::status_code(r), ": ",
         httr::content(r, "text", encoding = "UTF-8"), call. = FALSE)
  }
  httr::content(r, "parsed", encoding = "UTF-8")
}

message("creating a draft deposition on ", base)
dep <- api("POST", "/api/deposit/depositions",
           body = "{}", httr::content_type_json())
id <- dep$id
bucket <- dep$links$bucket

message("uploading ", basename(zip_path), " (",
        round(file.size(zip_path) / 1024^2, 1), " MB)")
r <- httr::PUT(paste0(bucket, "/", basename(zip_path)),
               httr::add_headers(Authorization = paste("Bearer", token)),
               body = httr::upload_file(zip_path))
if (httr::status_code(r) >= 300) {
  stop("Upload failed: ", httr::content(r, "text", encoding = "UTF-8"),
       call. = FALSE)
}

# The description is read by someone who found the record without the course, so
# it says what the objects are, what produced them, and that they are optional.
sha <- readLines("dist/cr_modelling_fits.sha256", warn = FALSE)[1]
description <- paste0(
  "<p>Fitted Bayesian concentration-response models for the training course ",
  "<em>Concentration-response modelling: estimating no-effect toxicity values ",
  "from concentration-response data in R</em>, taught at SETAC Australasia ",
  "2026.</p>",
  "<p>The course modules fit these models with the R package ",
  "<code>bayesnec</code>. Fitting them takes roughly half an hour in total, so ",
  "they are distributed here and the modules load them instead of sampling. ",
  "They are a convenience and nothing depends on them: a participant without ",
  "this archive fits the models as the modules are written.</p>",
  "<p>Extract into the root of the course folder so that the objects land in ",
  "<code>fits/</code>, or run <code>source(\"vignettes/fetch_fits.R\")</code> ",
  "from that folder. Source and modules: ",
  "<a href=\"https://github.com/open-AIMS/cr_modelling_training\">",
  "open-AIMS/cr_modelling_training</a>.</p>",
  "<p>sha256: <code>", sub("\\s.*$", "", sha), "</code></p>")

meta <- list(metadata = list(
  title = "Fitted models for the concentration-response modelling course",
  upload_type = "dataset",
  description = description,
  creators = list(list(name = "Fisher, Rebecca",
                       affiliation = "Australian Institute of Marine Science")),
  license = "cc-zero",
  keywords = list("ecotoxicology", "concentration-response", "bayesnec",
                  "no-effect concentration", "Bayesian"),
  related_identifiers = list(list(
    identifier = "https://github.com/open-AIMS/cr_modelling_training",
    relation = "isSupplementTo", scheme = "url"))
))

message("setting metadata")
api("PUT", paste0("/api/deposit/depositions/", id),
    body = jsonlite::toJSON(meta, auto_unbox = TRUE),
    httr::content_type_json())

cat("\n")
cat("Draft created, not published.\n")
cat("  review : ", base, "/deposit/", id, "\n", sep = "")
cat("  record : ", id, "\n", sep = "")
cat("\nPublish it in the browser after checking the metadata, then put\n")
cat("  ZENODO_RECORD <- \"", id, "\"\n", sep = "")
cat("at the top of vignettes/fetch_fits.R.\n")
cat("\nPublishing cannot be undone, which is why this script does not do it.\n")
