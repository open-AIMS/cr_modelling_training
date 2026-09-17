# Specification for the delivery format

Companion to `delivery-format-human.md`, which states the decisions. This file
holds the file layout, the generator, the per-module inventory, the edge cases
and the evidence. Every decision is stated in full here and in one sentence
there.

Written 2026-09-16. Every measurement below was taken on that date on a
four-core WSL2 machine unless stated otherwise.

---

# The participant bundle

## Layout

Participants receive one zip, not a git clone. A clone requires the fits bundle
to be fetched separately, requires git, and puts the module `.qmd` sources in
front of people who will open them instead of the live scripts.

```
cr-workshop/
  cr-workshop.Rproj
  fit_cache.R                 copy of vignettes/fit_cache.R
  live/
    02-fitting-a-model.R
    03-toxicity-estimation.R
    04-model-averaging.R
    05-response-distributions.R
    06-priors.R
    07-case-study.R
    08-groupings.R
  fits/
    <hash>.rds                written by cached()
    PROVENANCE.txt
  example_ogl.csv             the four tracked CSVs, verbatim
  example_pgl.csv
  example_fi.csv
  example_binomial.csv
  example_proportion.csv
  site/                       unzipped copy of docs/
    index.html
```

## The path invariant

The bundle root mirrors `vignettes/`, so every relative path in a module chunk
resolves unchanged from the bundle root. The module chunks are run with
`vignettes/` as the working directory under the render, and the data files sit
beside them.

This is what allows a live script to be verbatim module code. Any change that
breaks it (relocating the CSVs into a `data/` subdirectory, for instance) forces
the generator to rewrite paths, and a generator that rewrites code is a
generator that can produce code the module does not contain.

`fits/` is the one addition, and it is reached through the environment variable
rather than a relative path. See *Cache directory resolution*.

## Distributing a git clone instead

Rejected because `fits/` is git-ignored (correctly; the repository is public and
the objects are large), so a clone is never sufficient on its own. A participant
would need the clone and the bundle and would have to place one inside the
other. That is two steps to get wrong in a room with no troubleshooting slot.

---

# Live script generation

## Location of the extractor

`executable_chunks()` currently sits in `scripts/generate_taught_fits.R` at
lines 53 to 68. Move it, unchanged, to `scripts/chunk_extract.R` and have both
callers source it.

Do not reimplement it and do not substitute `knitr::purl()`. The existing
comment records why `purl()` cannot do this job: chunk option hooks are applied
when a document is knitted and not when it is purled, so `purl()` extracts an
`eval: false` block and the extracted script then runs it. Module 4's
`parallel-nested` chunk and module 6's `constant-prior` chunk are both
demonstrations that must not execute.

## Signature

```r
# scripts/chunk_extract.R
executable_chunks <- function(qmd)          # character vector of code lines
```

Unchanged from the current definition. It returns the body of every ```` ```{r} ````
block whose options do not include `eval: false`, with the `#|` option lines
stripped and a blank line between blocks.

## The generator

```r
# scripts/generate_live_scripts.R
#
# Writes one .R file per taught module into scripts/live/, from the module's
# own chunks. Run after any module edit; the scripts are generated artefacts
# and are never edited by hand.

source("scripts/chunk_extract.R")

MODULES <- c(
  "2Fitting-a-CR-model-using-bayesnec"            = "02-fitting-a-model",
  "3Toxicity_estimation_and_available_models"     = "03-toxicity-estimation",
  "4Model_averaging_and_multimodel_inference"     = "04-model-averaging",
  "5Response_data_and_statistical_distributions"  = "05-response-distributions",
  "6Priors_and_Bayesian_inference"                = "06-priors",
  "7Example_case_study"                           = "07-case-study",
  "8Factor_covariates_and_groupings"              = "08-groupings"
)
```

Each generated file opens with a fixed preamble and then the extracted chunks
verbatim:

```r
# ---------------------------------------------------------------------------
# Module 5, response data and statistical distributions.
#
# Generated from vignettes/5Response_data_and_statistical_distributions.qmd
# by scripts/generate_live_scripts.R on <date>. Do not edit; edit the module.
#
# Run from the bundle root, with fits/ beside this directory.
# The page for this module is at site/vignettes/<stem>.html
# ---------------------------------------------------------------------------

Sys.setenv(CR_FIT_CACHE = "fits")
source("fit_cache.R")
```

Module 1 has three chunks and no fits, and is not given a live script. It is
taught from the page.

## Section markers

Insert `# ---- <chunk label> ----` before each extracted block, taking the label
from the chunk's `#| label:` option where it has one. Positron and RStudio both
fold on that comment form, so the presenter can jump to a named section rather
than scrolling. A chunk without a label gets no marker rather than a generated
one; an invented name that does not appear on the page is worse than none.

Every fit chunk in modules 2, 4, 5 and 6 already has a label, listed in
*Fit call inventory*.

## Invariant

A generated script must run start to finish from the bundle root against a
populated `fits/` without sampling. Check this in CI-free fashion by running
each one after generation and asserting that no `Fitted and saved` message
appears. `cached()` emits `Loading a saved fit` on a hit and `Fitted and saved`
on a miss, so the two are distinguishable from the message stream.

---

# Wrapping the fit calls

This section is superseded as of 2026-09-17 and must not be implemented. The
`cached()` helper was written, used, and withdrawn. Each taught module now shows the
`bnec()` call and the `save()` beside it with `eval: false`, followed by the
`load()` that the render executes. The reason is in `intermediate-revision-human.md`,
section *Constraint on new code*: a course-specific caching abstraction hides
the fit-save-load sequence that participants are here to learn. `vignettes/fit_cache.R`
does not exist, and neither does the hash-keyed cache directory the file layout
above describes. The section is kept because it records what was tried and why
it was dropped, and because the fit call inventory below is still accurate about
which calls execute.

## The change to each call

Wrap the right-hand side of each executing fit call in `cached()`. The call
reads the same and the module gains one function name:

```r
# before
exp_2 <- bnec(resp ~ crf(sqrt.x, model = c("ecxll3", "nec3param")),
              data = beta_data, family = Beta(link = "identity"))

# after
exp_2 <- cached(bnec(resp ~ crf(sqrt.x, model = c("ecxll3", "nec3param")),
                     data = beta_data, family = Beta(link = "identity")))
```

Each taught module that fits needs `source("fit_cache.R")` in its setup chunk.
Under the render the working directory is `vignettes/`, so the bare filename
resolves; in the bundle it resolves from the root.

## Fit call inventory

Taken from the module sources on 2026-09-16 by parsing chunk bodies for
`bnec(`, `amend(` and `brm(`. Only the calls marked executing need wrapping;
the rest are demonstrations that the render never runs.

| Module | Line | Chunk label | Call | Executes |
|---|---|---|---|---|
| 2 | 141 | `seed` | `bnec(... nec3param ...)` | no |
| 2 | 302 | `fit` | `bnec(... nec3param ...)` | **yes** |
| 4 | 119 | `fit-nec3param` | `bnec(... nec3param ...)` | **yes** |
| 4 | 129 | `fit-ecxll3` | `bnec(... ecxll3 ...)` | **yes** |
| 4 | 161 | `amend` | `amend(bmanecfit, add = "ecxexp")` | **yes** |
| 4 | 272 | `parallel-basic` | `bnec(... decline ...)` | no |
| 4 | 293 | `parallel-nested` | `bnec(... decline ...)` | no |
| 4 | 542 | `ex9-simazine` | `bnec(... decline ...)` | no |
| 4 | 569 | `ex9-nassarius` | `bnec(... decline ...)` | no |
| 4 | 667 | `loo-controls` | `bnec(... decline ...)` | no |
| 5 | 186 | `fit-binom-nec` | `bnec(suc \| trials(tot) ...)` | **yes** |
| 5 | 245 | `fit-binom-ecx` | `bnec(suc \| trials(tot) ...)` | **yes** |
| 5 | 269 | `fit-betabinom` | `bnec(suc \| trials(tot) ...)` | **yes** |
| 5 | 336 | `fit-beta` | `bnec(resp ~ crf(sqrt.x ...))` | **yes** |
| 5 | 417 | `fit-pois` | `bnec(y ~ crf(x ...))` | **yes** |
| 5 | 467 | `fit-negbin` | `bnec(y ~ crf(x ...))` | **yes** |
| 5 | 494 | `fit-gamma` | `bnec(measure ~ crf(x ...))` | **yes** |
| 5 | 520 | `fit-gauss` | `bnec(y ~ crf(sqrt.x ...))` | **yes** |
| 6 | 293 | `fit-base` | `bnec(... nec3param ...)` | **yes** |
| 6 | 325 | `fit-a` | `bnec(... nec3param ...)` | **yes** |
| 6 | 345 | `fit-b` | `bnec(... nec3param ...)` | **yes** |
| 6 | 385 | `fit-c` | `bnec(... nec3param, nec4param ...)` | **yes** |
| 6 | 416 | `fit-d` | `bnec(... nec3param, nec4param ...)` | **yes** |
| 6 | 433 | `fit-e` | `amend(exmp_d, add = "necsigm", ...)` | **yes** |
| 6 | 455 | `constant-prior` | `bnec(... avoid_data ...)` | no |
| 6 | 543 | `sensitivity` | `bnec(... nec4param ...)` | no |
| 7 | 95 | `fit-ametryn-shown` | `bnec(fvfm ~ crf(concentration ...))` | no |
| 7 | 151 | `fit-all-shown` | `bnec(fvfm ~ crf(concentration ...))` | no |
| 8 | 250 | `ogl-fit` | `bnec(... ecxlin ...) + ogl()` | no |
| 8 | 310 | `pgl-fit` | `bnec(... nechorme ...) + pgl()` | no |
| 8 | 445 | `fi-fit-all` | `bnec(... decline ...)` | no |
| 8 | 474 | `fi-brms` | `brms::brm(bf_plain, ...)` | no |

Eighteen calls execute and need wrapping: one in module 2, three in module 4,
eight in module 5, six in module 6.

Modules 7 and 8 fit nothing at render time. Both read objects produced outside
it by `scripts/generate_herbicide_fits.R` and `scripts/generate_grouping_fits.R`,
which have not been run. Those two scripts write their own output and are
outside the `cached()` scheme; their results are added to the bundle as named
`.rds` files and the live scripts for those modules load them by name.

## A key that does not cover its inputs

`cached()` keys on the deparsed text of the call. For a `bnec()` call whose data
is a shipped dataset this is sufficient, because the text names the dataset and
the model set and the seed.

It is **not** sufficient for `amend()`. Module 4 line 161 is
`amend(bmanecfit, add = "ecxexp")`, and the key contains the name `bmanecfit`
and not its contents. If the set that `bmanecfit` holds changes, the key does
not change and a stale object is returned without a message.

Two mitigations, both required.

Pass an explicit key that names what the input was:
`cached(amend(bmanecfit, add = "ecxexp"), key = "m4-amend-ecxexp-from-nec3param-ecxll3")`.
Changing the upstream set then means changing the key by hand, which is visible
in a diff.

Delete `fits/` and regenerate whenever a module's model set changes, rather
than relying on the keys to notice. The regeneration is half an hour, and it is
the same discipline `CLAUDE.md` section 10 requires for `_freeze/`.

Module 6 line 433 is the other `amend()` call and has the same exposure.

## Cache directory resolution

`fit_cache_dir()` reads `CR_FIT_CACHE` and falls back to `file.path("..", "fits")`.

The fallback is correct under the render, where the working directory is
`vignettes/` and the bundle sits at the repository root. It is wrong from the
participant bundle root, where it points outside the bundle. This is why the
generated live scripts set `CR_FIT_CACHE = "fits"` in their preamble, before
sourcing `fit_cache.R`.

Do not change the fallback. Changing it to `"fits"` would break the render,
which is the more frequently exercised path.

## Consequence for generate_taught_fits.R

Once every executing fit call is wrapped, `scripts/generate_taught_fits.R` stops
needing its own `saveRDS()` loop. Rendering the site, or running the script,
populates `fits/` through `cached()` itself. Its own header anticipates this at
lines 14 to 19 and warns that the two routes write different filenames for the
same object, so do not leave both active against one bundle.

Reduce the script to: run each module's executable chunks, then write
`manifest.csv` and `PROVENANCE.txt` from what appeared in `fits/`. The manifest
no longer records the per-object name, because a `cached()` filename is a hash. Record the key alongside the hash instead, which means `cached()` gains
a line writing `<hash>.key` beside `<hash>.rds`.

---

# The runs-here callout

## Purpose

To stop a participant starting a fit that will not finish inside the session.
Stated on the page rather than only aloud, because a person working through the
material alone afterwards needs the same information.

## Wording

Placed immediately after the module title, before the first prose. Second
person is correct in body text here; see `CLAUDE.md` section 7 of this
repository.

```markdown
::: {.callout-note}
## Running this module

The fits in this module take about <N> minutes in total. Open
`live/<file>.R` and run it from the top. The fits are loaded from the
`fits/` folder in the workshop bundle, so nothing is sampled unless you have
deleted it.

Run these yourself: <the short calls>.
Watch rather than run: <the long calls, with their time>.
:::
```

Heading is a noun phrase, per `CLAUDE.md` section 12 rule 2. Do not write
*What you need to run*.

## Per-module content

Fill `<N>` from the render timings in *Timing evidence*, not from an estimate.
The content of the two lists follows from the triage decision in *Material that
does not fit*, which is not yet made, so the callouts are added after it.

---

# Bundle distribution

## Release asset and USB

A GitHub release asset gives a stable URL that can be named in
`vignettes/0Software-setup.qmd` as part of the pre-work. GitHub's per-file limit
for a release asset is 2 GB, well above anything this bundle will reach.

USB sticks cover the case where several dozen people download at a venue at
once. The zip also holds `site/`, an unzipped copy of `docs/`, so the material
opens from the file system if the network fails entirely.

## The published asset

The tag is `fits` and the asset is `cr_modelling_fits.zip`, so the address is

```
https://github.com/open-AIMS/cr_modelling_training/releases/download/fits/cr_modelling_fits.zip
```

and it is held in `FITS_URL` at the top of `vignettes/fetch_fits.R`. The tag
does not change. Rebuild with `Rscript scripts/bundle_fits.R` and replace the
asset with

```bash
gh release upload fits dist/cr_modelling_fits.zip --clobber -R open-AIMS/cr_modelling_training
```

which leaves the address valid, so no module and no script needs editing when
the fits are rebuilt. `scripts/bundle_fits.R` also rewrites `vignettes/fits.sha256`,
which is committed and is what a participant's download is checked against;
commit it in the same change as the upload or the check reports a mismatch that
is not one.

## Measured size

42.3 MB over 25 objects and `PROVENANCE.txt`, measured 2026-09-17. The archive
stores its paths as `vignettes/fits/<object>.RData`, so it extracts correctly at
the root of a `cr_modelling_training` folder and nowhere else. An end-to-end run
of `source("vignettes/fetch_fits.R")` on that date downloaded and unpacked all
25 objects in 32 seconds.

The zip is reproducible. Rebuilding it from unchanged objects on 2026-09-17
produced the same sha256, `0a1186c2...`, as the previous build, so a rebuild
that changes nothing does not invalidate a participant's copy.

For scale, the 2023 objects it replaces are `fitted_model.RData` at 124 MB,
`manfecfit.RData` at 451 MB and `ametryn.RData` at 65 MB.

## The rejected institutional share

A OneDrive share link on the AIMS tenant was used between 2026-09-17 and the
release asset above, and it does not work for participants. Fetched with `curl`
from a client not signed in to the tenant, the `&download=1` form returned 302
to the file path and then 403 with `Access denied. Before opening files in this
location, you must first browse to the web site and select the option to login
automatically`. The form without `&download=1` returns 200 with 58 KB of HTML,
which lands in a file named `.zip` and fails at `unzip` rather than at the
download.

Neither failure is visible to the person who created the share, because their
browser is signed in. Test a hosting change by fetching it unauthenticated:

```bash
curl -sIL '<url>' | grep -iE '^HTTP|content-type|content-length'
```

## Zenodo

Still the destination for a citable identifier once the material is frozen, and
`scripts/zenodo_deposit.R` uploads and sets the metadata without publishing. It
is not the route for the workshop, because the archive is still being rebuilt as
modules are revised and a Zenodo record is versioned rather than replaced.

## Regeneration trigger

Regenerate after any upgrade to `bayesnec`, `brms`, `cmdstanr` or CmdStan, for
the reason `CLAUDE.md` section 10 gives for `_freeze/`: a stored result does not
know that the package which produced it has changed. `PROVENANCE.txt` already
records all four versions, so the check is a comparison against the installed
set.

---

# Material that does not fit

## The measurement

Chunk counts, taken on 2026-09-16 by counting ```` ```{r} ```` openers and
`eval: false` options in each module source.

| Module | Words | Chunks | `eval: false` | Executing |
|---|---|---|---|---|
| 1 The software stack | 1,859 | 3 | 1 | 2 |
| 2 Fitting a single model | 4,938 | 26 | 5 | 21 |
| 3 Toxicity estimation | 3,487 | 18 | 1 | 17 |
| 4 Model averaging | 4,470 | 33 | 6 | 27 |
| 5 Response distributions | 5,406 | 37 | 4 | 33 |
| 6 Priors | 3,660 | 20 | 3 | 17 |
| 7 Case study | 2,431 | 21 | 4 | 17 |
| 8 Groupings | 3,759 | 29 | 8 | 21 |
| **Total** | **29,810** | **187** | **32** | **155** |
| drc reference (not taught) | 3,292 | 22 | — | — |

Word counts are `wc -w` on the `.qmd` source. They count backticked
identifiers, chunk contents and table syntax, so they overstate the prose by a
margin that varies with how much code a module holds. `CLAUDE.md` section 13
records this measurement error from the earlier revision. Use them as a relative
measure of module length only.

## The arithmetic

A nine-to-five day with an hour for lunch and two twenty-minute breaks gives 380
minutes of teaching. Against 155 executing chunks that is 2.45 minutes per
chunk if the whole day went to running code and none to explanation. Allowing
half the day for explanation leaves about 1.2 minutes per chunk, which is less
than the time to introduce a call, run it, and let several dozen people run it.

## Modules against the course outcomes

From the six course outcomes in `index.qmd`.

| Outcome | Served by |
|---|---|
| distinguish linear, non-linear and threshold models | 3 |
| identify the distribution for a response and fit accordingly | 5 |
| fit models in R and plot the results | 2 |
| derive NEC, NSEC and ECx and explain each | 3 |
| apply model averaging and interpret the weights | 4 |
| produce complete workflows for SSD derivation | 7 |

Module 1 serves no outcome directly and enables all of them. Module 8 serves
none of the six. The `drc` reference is already marked as not taught.

## Options

Not decided here. Each is recorded with what it gives up.

### Demoting module 8 to reading

Removes 21 executing chunks and the two
generated fit sets behind it. Gives up the only treatment of comparing
estimates across groups, which is a question participants ask even though it is
not a listed outcome. Cheapest option by a wide margin and the first to
consider.

### Compressing module 1 to ten minutes

It is the shortest module at 1,859 words
and two executing chunks, and its content overlaps the setup pre-work that
participants have already done. Gives up little; the risk is that the
participants who did not do the pre-work are then unrecoverable, which they
already are.

### Splitting module 5 into a taught half and a reference half

It is the longest
module and holds eight of the eighteen executing fits. Teach the binomial and
beta families, which cover most ecotoxicological endpoints, and leave Poisson,
negative binomial, gamma and Gaussian as worked examples on the page. Gives up
the completeness that makes the module useful afterwards, which argues for
keeping the page whole and teaching part of it rather than splitting the file.

### Reducing module 4's parallel fitting section

Added per `CLAUDE.md` section 8
and not yet delivered, so its length is still under the presenter's control.

Whichever combination is chosen, the outcome is a stated number of minutes per
module summing to under 380, recorded in this file, and the runs-here callouts
are written from it.

---

# Timing evidence

Rendering the development site on 2026-09-16, on a four-core WSL2 machine with
`options(brms.backend = "cmdstanr", mc.cores = 4)`:

| Module | Render time |
|---|---|
| 2 | 1.5 minutes |
| 4 | 4.5 minutes |
| 5 | 12.9 minutes |

Recorded in the header of `scripts/generate_taught_fits.R`. Module 6 was not
timed and holds six executing fits, so it should be timed before the callouts
are written. Module 3 fits nothing.

These are render times on the presenter's machine and are not what a
participant's machine will take. They are the argument for the bundle, not an
estimate of what anyone waits.

---

# Rejected alternatives

## revealjs from the same source

Adding `revealjs` as a second format target on the module `.qmd` files was
considered. Rejected because a slide break in Quarto's reveal format falls at a
level-2 heading by default, and the modules have sections running to several
hundred words under one heading. Producing usable slides means inserting breaks
and shortening prose, at which point the two outputs diverge in content and the
single-source argument no longer holds.

## quarto-live and webR

Settled in `CLAUDE.md` section 2 before this note. webR has no C++ toolchain, so
`cmdstanr` cannot compile a model and `brms` cannot run. The Bayesian half of
the course would be non-executable code display. Recorded here so it is not
reopened.

## Reviving the learnr set

The learnr documents are retired by the decision in `CLAUDE.md` section 2. They
also require a running Shiny server, which reintroduces a failure mode during
the session that the static site does not have.

## A deck per module

Rejected because both versions would have to be maintained and would diverge.
Eight modules is roughly two hundred
slides, built against material that is being revised, with three corrections and
two unbuilt modules outstanding.

## Hand-written live scripts

Rejected because a hand-written script drifts from the module the first time
either is edited, and the drift is silent: the script still runs, and produces
something other than what is on the screen.

---

# Open items

1. Module 6 render time, needed for its callout.
2. Bundle size, from the first run of `scripts/generate_taught_fits.R`.
3. The triage decision in *Material that does not fit*, and the per-module
   minutes that follow from it.
4. Whether `scripts/generate_herbicide_fits.R` and
   `scripts/generate_grouping_fits.R` have run, which determines whether modules
   7 and 8 exist in a taught form at all. `index.qmd` currently states that they
   do not.
5. The single published URL for the day, which requires the approved build to
   hold every taught module and the `/dev/` split to be retired.
