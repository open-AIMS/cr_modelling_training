# CLAUDE.md — Concentration-response modelling training course

Project-specific instructions. `C:/Rworking/CLAUDE.md` applies in full; this file adds
context and overrides only where §2 states an override explicitly.

---

## 1. Purpose and audience

A **training course repository**. It is not an R package and not a research compendium:
there is no `DESCRIPTION`, no `R/`, no `NAMESPACE` and no `R CMD check`. Everything in
`vignettes/` is teaching material, despite the folder name.

The course covers the estimation of no-effect toxicity values from concentration-response
(CR) data, by two engines in parallel: frequentist fits with `drc` and Bayesian fits with
`bayesnec`. The learning outcomes are listed in `README.md` and are the specification for
the content — a module that does not serve one of them needs a reason to exist.

**Event:** SETAC Australasia, September 2026. Delivered as a **full-day workshop**: all
eight modules, participants running code on their own machines.

**Audience:** practising ecotoxicologists and risk assessors. Assume competent R users who
are not statisticians and have not used Stan. Most will be meeting Bayesian inference for
the first time, and a substantial number will arrive with a broken C++ toolchain.

**Presenter:** Dr Rebecca Fisher, Australian Institute of Marine Science.

**The repository is public** (`open-AIMS/cr_modelling_training`).

---

## 2. Repository type and the Quarto migration

The course is a Quarto website. The migration decided on 2026-09-11 is complete: every
module is a `.qmd` in `vignettes/`, the `learnr` set is retired, and the parent §9
prohibition on R Markdown applies here without an exception. Do not author a new `.Rmd`
in this repo. Two 2023 `.Rmd` files remain tracked and are not course content,
`0Overview.Rmd` and `Scratch.Rmd`.

**Every module stays in `vignettes/`.** The course is a series of vignettes and that is how
collaborators expect to find it. A converted module replaces its `.Rmd` in place, keeping
the same file stem, so the numbering and any existing reference still resolve. Do not move
a module into a new directory; module 1 was briefly split into a `setup/` folder on
2026-09-11 and moved back the same day.

**Software setup is pre-work and sits outside the numbered sequence.** It is
`vignettes/0Software-setup.qmd`, with `vignettes/check_setup.R` beside it. It is issued
before the workshop and is not taught: a participant who has not completed it cannot fit
models on the day, and there is no troubleshooting slot. Module 1 is
`vignettes/1Software-stack.qmd`, which is taught, and covers what R, Stan, `brms` and
`bayesnec` each contribute and how Positron is used.

**Converting a module touches three places**, and missing one leaves the module
unreachable: the `project: render:` list and the sidebar in `_quarto.yml`, and the module
table in `index.qmd`.

**The runtime is a static Quarto site, and participants run the code on their own
machines.** `quarto-live` and webR were rejected because webR has no C++ toolchain, so
Stan cannot compile a model in the browser and the Bayesian half of the course would
become non-executable code display. A `server: shiny` build was rejected with the
`learnr` set. `notes/delivery-format-human.md` holds the reasoning and what follows from
it for the day.

**The deployment target is GitHub Pages**, served from `docs/` on `main`, at
<https://open-aims.github.io/cr_modelling_training/>. Two builds are published from the
one source. `_quarto.yml` builds the approved site into `docs/` and `_quarto-dev.yml`
builds everything, including modules under revision, into `docs/dev/`. Render the
approved site first, because it clears `docs/` and would take the development site with
it:

```bash
quarto render                 # -> docs/
quarto render --profile dev   # -> docs/dev/
```

As of 2026-09-19 the approved site holds every taught module, the `drc` reference and the
software setup. The modules were promoted from `/dev/` that day: the `render:` list and
the sidebar in `_quarto.yml`, and the module table in `index.qmd`, whose
`content-visible` split into a linked and an unlinked copy was removed, both profiles now
showing the same linked table. The `/dev/` build is still produced and still carries the
not-reviewed banner, and is now a place to put a module under revision rather than a
second copy of the course. `notes/delivery-format-human.md` records that the split should
not survive to the day.

---

## 3. The module set

Numbered modules are the Bayesian thread. The `drc` counterparts that used to sit beside
them, with a `d` suffix, were consolidated into one reference page in `13d386d`.

| # | File stem | Title | Notes |
|---|---|---|---|
| — | `0Software-setup` | Software setup | pre-work, not taught, unnumbered |
| 1 | `1Software-stack` | The software stack | R, Stan, `brms`, `bayesnec`, Positron |
| 2 | `2Fitting-a-CR-model-using-bayesnec` | Fitting a single model | |
| 3 | `3Toxicity_estimation_and_available_models` | Toxicity estimation and the model set | |
| 4 | `4Model_averaging_and_multimodel_inference` | Model averaging and multimodel inference | |
| 5 | `5Response_data_and_statistical_distributions` | Response data and statistical distributions | |
| 6 | `6Priors_and_Bayesian_inference` | Priors and Bayesian inference | |
| 7 | `7Example_case_study` | A worked case study | reads results written by `scripts/generate_herbicide_fits.R` |
| 8 | `8Factor_covariates_and_groupings` | Factor covariates and groupings | reads results written by `scripts/generate_grouping_fits.R` |
| — | `drc-reference` | Fitting the same models with `drc` | reference, not taught |

All of these are `.qmd` in `vignettes/`. Every module holds its executable code as the
`bnec()` call and its `save()` shown but not run, followed by the `load()` that the
render executes, so a render samples nothing and needs `vignettes/fits/` (§5).

`2d Fitting-a-CR-model-using-drc.Rmd`, `3dToxicity_estimation_and_available_models.Rmd`,
`4dModel_averaging_and_multimodel_inference.Rmd`, `8dFactor_covariates_and_groupings.Rmd`
and `test_4dModel_averaging_and_multimodel_inference.Rmd` were deleted in `13d386d` and
their material is in `drc-reference.qmd`. `7dExample_case_study.Rmd` was deleted in 2023
(`8d55cd2`) and is recoverable with `git show 47abfcd:vignettes/7dExample_case_study.Rmd`.

`0Overview.Rmd` is the retired 2023 course outline, still tracked, still carrying module
1's old title. `index.qmd` has replaced it.

`vignettes/Scratch.Rmd` is a 2023 working file, not a course module. `scratch/` and
`ignore/` are likewise not course content; `ignore/` is git-ignored in full.

`vignettes/functions.R` holds three helpers written for module 8 in 2023, before the
equivalent functionality existed in the packages: `pred_out()`, `nec.brmsfit()` and
`nsec.brmsfit()`. Nothing sources it as of 2026-09-17; the revised module 8 uses
`bnec_group()`. `nsec.brmsfit()` calls three unexported `bayesnec` internals
(`bayesnec:::do_wrapper`, `bayesnec:::modify_posterior`, `bayesnec:::min_abs`), which
break without deprecation, so check it against the current `bayesnec` before reviving
anything from it.

---

## 4. Packages in scope

No `DESCRIPTION` and no `packages.R`, so parent §4 gives no canonical list. The list below
was taken from the `library()` calls and `::` usage across `vignettes/*.qmd` and
`scripts/*.R` on 2026-09-17 and is the working set until a `packages.R` is created:

`bayesnec`, `brms`, `cmdstanr`, `rstan`, `drc`, `future`, `future.apply`, `posterior`,
`loo`, `knitr`, `ggplot2`, `dplyr`, `tidyr`, `purrr`, `car`, `digest`, `remotes`,
`usethis`.

`learnr`, `qwraps2`, `bayesplot`, `scales`, `ggpubr` and `cowplot` were on the 2026-09-11
list and no module uses them now. `vignettes/check_setup.R` installs and verifies a
narrower set again, and it is the contract a participant is checked against.

**Create `packages.R` before the course.** A full-day workshop where participants install
packages needs one authoritative list, and the three places above are not one. Ask before
adding anything to it.

**The dead `extractNSEC` dependency is gone.** It was called by the `2d` module, which
`13d386d` deleted. `drc-reference.qmd` derives the *NSEC* by hand from the fitted curve
and depends on no such package. `toxval` exports `nsec()` with an `nsec.drc` method if a
future page needs one.

`bayesnec` is installed from the `dev` branch by `0Software-setup.qmd:149`
(`remotes::install_github("open-AIMS/bayesnec", ref = "dev")`). **Pin this to a release or
a commit before the workshop.** A `dev` branch that changes during the course produces
failures that cannot be diagnosed in the room, and the modules call `bnec_group()`,
`check_sampling()`, `screen_models()`, `failed_models()`, `bnec_record()`,
`curve_params()`, `ecnsec()` and `check_fit()`. The branch is moving: the site was
rendered against 2.1.3.37 (`15c12765`) and `dev` was at 2.1.3.39 later the same day.

**`lum31` is on `dev`.** Pull request #228 merged on 2026-09-17, so module 8 and
`scripts/generate_grouping_fits.R` run against `dev` without a branch install. `6f36a2e`
records the correction.

Where a branch does have to be installed, install it into a private library and set
`R_LIBS_USER` for that session's R calls. On 2026-09-17, before the merge, two sessions
overwrote each other's `bayesnec` twice inside an hour, the second time removing `lum31`
from under a running script.

---

## 5. Data files and their availability

**The fitted objects live in `vignettes/fits/`, which is git-ignored, and are
distributed as a release asset.** 26 objects, 43.1 MB zipped, built by
`scripts/bundle_fits.R` into `dist/` and published under the fixed tag `fits`:

```bash
Rscript scripts/bundle_fits.R
gh release upload fits dist/cr_modelling_fits.zip --clobber -R open-AIMS/cr_modelling_training
```

The tag does not change, so `FITS_URL` in `vignettes/fetch_fits.R` stays valid.
`scripts/bundle_fits.R` also rewrites `vignettes/fits.sha256`, which is tracked and is
what a participant's download is verified against, so commit it with the upload. `N_FITS`
in `fetch_fits.R` is recounted whenever the archive gains an object. A participant runs
`source("vignettes/fetch_fits.R")`, which took 32 seconds on 2026-09-17.

Rebuild and re-upload after any change to a fit call, and after any upgrade to `bayesnec`,
`brms` or Stan, for the reason §10 gives for `_freeze/`: a saved object does not know that
the package which produced it has changed. `vignettes/fits/PROVENANCE.txt` records the
version each object was fitted under and travels inside the archive.

**A fresh clone cannot render any module that fits** until the archive is fetched, because
each fit is shown as a call and a `save()` that do not run, followed by a `load()` that
does. That is the intended arrangement, not a defect: the render samples nothing and takes
minutes rather than half an hour.

Three 2023 objects are still in `vignettes/` and nothing reads them any more:
`fitted_model.RData` (124 MB), `manfecfit.RData` (451 MB) and `ametryn.RData` (65 MB).
They are git-ignored and local to the presenter's machine. Do not commit them; the
repository is public and 640 MB of `.RData` is not a reasonable clone.

Five CSV files are tracked and small. Module 8 reads `example_ogl.csv` and
`example_fi.csv`; modules 5 and 6 read `example_binomial.csv` and `example_proportion.csv`.
`vignettes/data/` holds the results that modules 7 and 8 read back from
`scripts/generate_herbicide_fits.R` and `scripts/generate_grouping_fits.R`. `example_pgl.csv` is no longer
read by anything as of 2026-09-17: it held the Lum-31 bioluminescence data after the
per-plate control division and a further division by the plate maximum, over four plates,
and module 8 now uses the recorded readings from `lum31` over sixteen. It is kept rather
than deleted because the normalised form is what Luter et al. (2025) published and a later module
may want to show the two side by side.

`vignettes/rsconnect/` is 4 KB and stays where it is. It is the shinyapps.io deployment
record for module 1 (appId 9520984, https://open-aims.shinyapps.io/1Getting-started/) and
the only link between this repository and the deployed application. Only module 1 was ever
deployed. Removing it would make `rsconnect` publish a new application rather than update
that one.

`vignettes/images/` holds 48 tracked images. Several are **screenshots of publisher-typeset
journal pages**: `etnc_header.jpg` is the ET&C title block of Fisher and Fox (2023),
checked on 2026-09-11, and `etnc_fig1.jpg`, `etnc_fig2.jpg`, `NSEC_ieam.jpg`,
`ieam_head.jpg`, `necmod_fox2010.jpg`, `modelave_ecol.jpg` and `glmbooks.jpg` are named as
the same kind of capture. Some are the presenter's own papers, which does not by itself
settle the position, because copyright in the typeset version usually sits with the
publisher. The repository is public. Each was checked against its Crossref licence record
on 2026-09-19, and RF decided on 2026-09-20 that every capture stays;
`notes/image-provenance.md` holds the licence of each source and the decision. Do not
reopen it. For a figure added from now on, still prefer redrawing from the underlying
quantities over reproducing a publisher's rendering.

---

## 6. What may be published here

**The repository is public, and nothing from client or unpublished work in this working
directory may reach it.** That includes project and client names, dataset names peculiar to
an engagement, and measurements taken from such work — a measured percentage change is as
identifying as a name once someone has the report.

Other repositories under `C:/Rworking/` are a legitimate source of *understanding*: they
record traps and engine differences that took real work to establish, and reading them
avoids rediscovering the same things. What may cross over is the mechanism, not the
evidence.

**Every number on a page is produced by that page, or comes from a published source.** Where
a claim was learned from unpublished work, either demonstrate it here on data the course
already ships, or verify it against the public package's own source and say so. Both routes
were used on the `drc` reference page: the ECx divergence is demonstrated on `nec_data`, and
the `maED()` limitations are verifiable in `drc::maED` itself.

This was not hypothetical. An earlier version of this file named a client re-analysis
project and its date, in prose about writing style, and that text was committed and pushed.

## 7. Prose and register

Parent §12 applies in full: plain scientific register, no idiom, Australian spelling, and
the ruled corrections for *carry*, *move*, *cost*, *help/hurt* and *cut both ways*. Parent
§12 "Headings and register in long-form documents" applies to every module. What follows is
additional, and is specific to teaching material.

**The 2023 text is the register reference.** The existing modules were written by hand and
are already close to the target. Read the surrounding paragraphs before writing a new one
and match them. Where a revision and the existing text disagree on register, the existing
text wins unless it breaks a parent §12 rule.

**Second person is correct in body text and wrong in headings.** A workshop document
addresses a participant directly: *you will need Rtools before this will run*, *fit the
model and then plot it* are right. A heading names a thing: *Toolchain requirements*, not
*What you need before you start*. This is the one place the teaching register differs from
the report register in parent §12, and it differs only in the body.

**First person plural is the existing convention and is kept.** *We have also allowed a
linear decline in `bayesnec`* is the voice of the package authors, who are presenting. Do
not convert the modules to the impersonal passive of parent §12; that rule is for reporting
results, and these documents are instruction. Use the passive for a measurement within a
module.

### Patterns to remove on sight

These are the ones that make revised teaching prose read as machine-written. Every one has
been corrected elsewhere in this working directory before, the `bayesnec` vignettes
included.

| Pattern | Instead |
|---|---|
| the tricolon — three parallel clauses where two would do | two clauses, or one |
| *not X, but Y* | state Y |
| *It is worth noting that* / *Importantly* / *Crucially* / *Note that* | delete, and state the thing |
| an em-dash aside mid-sentence | a comma, a colon, brackets, or a second sentence |
| an aphoristic one-line summary closing a section | delete; the section already said it |
| meta-commentary — *this module establishes*, *as we saw above*, *the key takeaway is* | delete, or make it a forward pointer with a section name |
| bold sprinkled through running prose | bold only for a defined term at first use, or a figure or table label |
| *simply*, *just*, *straightforward*, *powerful*, *robust*, *seamless* | delete |
| *delve into*, *dive into*, *unpack*, *leverage* | *describe*, *examine*, *use* |
| a rhetorical question as a section opener | the statement it was standing in for |
| *In this section we will…* announcing what follows | do the thing |

**Do not tell a participant that something is easy.** *Simply run `bnec()`* is wrong twice
over: it is filler, and it tells someone whose model has not converged that they have
failed at something easy. State what to run, and what to do when it does not work.

**A claim about a package's behaviour is checked before it is written.** These modules are
teaching material, so an inaccuracy propagates into other people's analyses. Run the code,
or cite the help page or the paper. The traps that are counter-intuitive and easy to state
wrongly are recorded in the repository each belongs to, which parent §11 indexes:
`C:/Rworking/bayesnec/CLAUDE.md` for the identity link forced on every family, estimates
and plot data on different scales, and `dispersion()`; `C:/Rworking/CR_workflows/CLAUDE.md`
for the `drc` standard error that is unavailable over part of the range and the variance
fixed at the mean under `type = "Poisson"`. Neither file loads in a session started here,
so read them before writing about any of those.

**Name model equations precisely.** Parent §12 rule 6, "Name models precisely", applies with
particular force here: this course is where participants learn the vocabulary. Never write
"the nec model" for a specific equation; write `nec3param`, `nec4param`, `nechorme`. Keep
the three senses of `ecx` distinct — the equation group, the `ecx()` function, and an ECx
estimate.

**Exercise and solution text.** An exercise prompt states the task and nothing else. A
solution states what the code does and why that answer, not *Great — you have now
mastered…*. No congratulation, no encouragement written in advance of the participant
doing anything.

---

## 8. Outstanding work

Both items recorded here on 2026-09-11 are done: module 4 covers parallel fitting under a
`future` plan, and the `drc` comparison was rewritten as `vignettes/drc-reference.qmd`.
What is outstanding, as of 2026-09-19:

**`bayesnec` is installed from a moving branch** (§4). Pin it.

**The runs-here callout is not written.** `notes/delivery-format-human.md` step 5 asks each
taught module to open with a note naming the code a participant runs and the code they
watch, so that nobody starts a fit that will not finish before the module ends. No module
has one. The triage it depends on is settled: `notes/workshop-agenda.md` gives every block
a stated number of minutes against the real 285. The other two steps that waited on that
triage are done. `slides/opening.qmd` is the opening deck, and `scripts/live/` holds one
generated `.R` file per taught module.

**The live scripts are generated, not written.** `scripts/generate_live_scripts.R` builds
`scripts/live/` from the module sources, using the chunk parser in `scripts/qmd_chunks.R`
that `scripts/generate_taught_fits.R` also uses. Regenerate after any change to a module;
a hand edit to a file under `scripts/live/` is lost at the next run and leaves the script
disagreeing with the page while the room is looking at both.

Which fit calls run in the room rather than being commented out is decided from
`scripts/fit_times.csv`, the sampling time of each saved object, written by
`scripts/measure_fit_times.R`. **Re-run that after any rebuild of the fits and commit it
with the archive**, for the reason §5 gives: a recorded time does not know that the call
or the package which produced it has changed. It is tracked while `vignettes/fits/` is
git-ignored, so that generating on a clone without the bundle gives the same result. The
file being absent is not an error; every fit call is commented out instead.

The estimate a fit is judged against adds `COMPILE_SECONDS`, 35, for each distinct
equation, because Stan compilation is not recorded in a saved object and is most of the
wall time of a small fit: module 2's fit samples in 2.8 seconds and takes 38.5. Re-measure
it after a toolchain change. Where a fit does run, its `save()` is commented out, or the
script would overwrite the distributed objects and `vignettes/fits.sha256` would no longer
match what a participant downloaded.

**`packages.R` does not exist** (§4).

**The publisher-typeset figures are resolved.** Every source was looked up in its
Crossref record on 2026-09-19, and on 2026-09-20 RF decided that every capture stays,
including the three whose source registers no reuse licence. Nothing is to be removed or
redrawn on this account. `notes/image-provenance.md` is now the provenance record rather
than a list of work outstanding, and holds the licence of each source.

## 9. Known defects in the current material

Found by inspection on 2026-09-11. The learnr-era defects are resolved by the migration:
the missing `css/style.css` (the site uses `styles.css`), the `4d` module that was not a
learnr tutorial, the deleted `7d` source (§3), the "gruopings" title typo, the README
"backage" typo, and the dead `extractNSEC` dependency (§4).

What remains:

- **`0Overview.Rmd` still carries module 1's 2023 title**, "Getting started - Installing
  and running BRMS", and its content is the course outline that `index.qmd` now holds. It
  is tracked and published to nobody. Delete it or leave it; do not revise it.
- **`bayesnec` is installed from `dev`** — §4.
- The 2023 material was last run in November 2023. `bayesnec`, `brms` and `drc` have all
  released new versions since, so assume nothing in `0Overview.Rmd`, `Scratch.Rmd` or
  `scratch/` runs until it has been run.

---

## 10. Working in this repository

**Line endings.** The working tree is CRLF and the index is LF, so `git status` reports all
26 tracked text files as modified while `git diff --ignore-cr-at-eol` reports nothing. The
diff is line endings only, verified 2026-09-11. Do not "fix" this by committing the CRLF
versions — that rewrites every file and destroys the diff of the actual revision. Add a
`.gitattributes` normalising the repository to LF, as `CR_workflows` does, before the
migration begins.

Because of this, parent §7's git-status rule needs care here: a file reported as modified
is not necessarily a file with uncommitted work. Check with `git diff --ignore-cr-at-eol`
before deciding, and warn only if that shows content.

**R.** Parent §1 applies. Confirm which R is active before running anything; the Bayesian
modules need a working C++ toolchain and Stan, which WSL and Windows do not share.

**Rendering.** A learnr document is rendered by `rmarkdown::run()` rather than by knit,
because of `runtime: shiny_prerendered`.

**Build artefacts were moved out of `vignettes/` on 2026-09-11.** The knitr caches
(13.8 GB), the figure `*_files/` directories, the rendered HTML and the diagnostic PDFs are
now under `ignore/render_artefacts_2023/`, together with a `PROVENANCE.md` recording what
each item is, why it was moved and how to reverse the move. The 2023 offline installer
bundle is now `ignore/software.zip`. Nothing was deleted and none of it was ever tracked.
Two consequences:

- Re-rendering a learnr module now recomputes every chunk, which means refitting the Stan
  models. Reverse the move before rendering a 2023 module.
- `scratch/interactive.R` reloads the caches with `qwraps2::lazyload_cache_dir()` and its
  paths no longer resolve. Repoint them at `ignore/render_artefacts_2023/` or reverse the
  move.

**`ignore/` is the local-only holding area** for material that is not course content:
superseded drafts, background reading, the installer bundle and the 2023 build artefacts.
It is git-ignored in full, so nothing placed there is backed up by the remote. This is the
repository's equivalent of the `superceded/` folder in parent §7.

**Review drafts render to `ignore/render_artefacts_2026/`.** When a module or vignette is
rendered ad hoc so a person can read it (as opposed to a `quarto render` of the published
site), put the output HTML there rather than in a system temp directory, since the user
cannot easily open temp paths. Created 2026-09-15 with a render of
`vignettes/0Software-setup.qmd`. Keep this separate from `ignore/render_artefacts_2023/`,
which holds the retired learnr build artefacts, not review drafts.

**Freeze does not notice a package upgrade.** `execute: freeze: auto` keys on the source
document, not on the environment that rendered it. Upgrading `bayesnec`, `brms` or Stan
leaves every frozen result in place, so the site keeps publishing output from the previous
version with no warning and no diff. Clear the stored result for the affected modules and
re-render:

```bash
rm -rf _freeze/vignettes/<module>
quarto render
```

Measured on 2026-09-11: module 2 was rendered under `bayesnec` 2.1.3.7, the package was
upgraded to 2.1.3.33, and a plain `quarto render` reproduced the 2.1.3.7 output unchanged.
The default `resolution` had changed from 1000 to 200 between those versions, so the
published page would have reported a superseded default as current.

After any package upgrade, clear `_freeze/` for every module that fits a model.

**Freeze does not notice a prose change either.** The frozen `result.markdown` in
`_freeze/vignettes/<module>/execute-results/html.json` holds the *whole* rendered document,
prose included, and `freeze: auto` serves it without consulting the source. A correction to
the text of a module can therefore render, report success, and publish the old wording,
with nothing in the log to say so. Measured 2026-09-19: module 2's sampling bullet was
corrected, two renders reported `approved OK`, and both published the superseded sentence;
the frozen markdown held the old text and not the new.

Clear the module's stored result after editing its prose, the same as after a package
upgrade, and check the rendered page for the words that changed:

```bash
rm -rf _freeze/vignettes/<module>
quarto render && quarto render --profile dev
grep -c "<a phrase you added>" docs/vignettes/<module>.html
```

This is the more dangerous of the two freeze traps, because a prose edit gives no reason to
suspect the cache and the rendered page looks finished.

**A `future` multisession plan deadlocked the generation scripts.** `bnec()` fits a model
set under whatever plan is active (`bayesnec` #184), and on 2026-09-17
`scripts/generate_grouping_fits.R` hung under `plan(multisession, workers = 5)` on the
copper plate set: the parent had 439 bytes queued to one worker and that worker 415 queued
back, both in `poll_schedule_timeout`, with no Stan process running and no CPU used by any
of the six processes for fourteen minutes. It is silent, so a script that produces nothing
for half an hour is worth checking with `ss -tnp | grep <port>` before assuming it is
sampling. The scripts here now fit in sequence, with the four chains of each fit in
parallel through `mc.cores`.

**Long fits are fetched rather than re-run where `bayesnec` has already run them.** The
`open-AIMS/grouping-structures` compendium runs the 189 fits behind `vignette("example8")`
as cluster array tasks and keeps the assembled objects in a store on the AIMS HPC, keyed on
the fit call and a digest of the data. `hpc/fetch-store.sh` in that repository brings the
whole store back; a single object is `rsync`ed from
`/export/scratch/rfisher/grouping-structures/store/`, and `store/index.csv` gives the key
for each call. Module 8's copper-against-zinc `bnec_group()` fit comes from there, through
the `LUM31_TOX_FIT` environment variable. A key answers a call only where the call and the
data both match, so a change to either means refitting.

**`autoplot()` draws a single-equation fit at a tenth the resolution of a model
average.** `ggbnec_data()` takes a `bayesmanecfit`'s curve from `w_pred_vals`,
which `bnec()` computed at its own `resolution` of 1000, and a `bayesnecfit`'s
from `brms::conditional_effects()`, whose default resolution is 100. Both grids
are evenly spaced on the **recorded** predictor rather than on the scale the
model was fitted on. A fit written as `crf(log(conc), ...)` and drawn on a log
axis is therefore drawn at very low resolution over its lower decades: measured
on 2026-09-17, over the `lum31` copper range of 0.002 to 10.2 mg/L, the
`bayesnecfit` frame held one grid point below 0.1 mg/L and the `bayesmanecfit`
frame for the same data held ten. The single fit's `ecxll4` curve was drawn as
one straight segment across the whole lower decade and appeared to decline
steadily where the equation is flat. Module 8 draws its single-equation curves
from `fit$pred_vals$data` instead. Worth raising as a `bayesnec` issue.

**The site is built against the head of `bayesnec`'s `dev` branch.** Check `dev` before
rendering, reinstall if it has moved, and record the commit used in the table in
`notes/setup-evidence.md`. The commit is recorded rather than pinned, so a published page
stays attributable without holding the course behind the package. Module 1 is the exception:
it tells participants which version to install, so it names a release rather than a branch.

**Every module that fits sets the Stan backend.** `brms` uses `rstan` unless
`options(brms.backend = "cmdstanr")` is set, and `rstan` compiles a model a good deal more
slowly. Measured 2026-09-19: module 2's two fits took 199 seconds under the `rstan` default
against 73 under `cmdstanr` on the same four-core machine, almost all of the difference
being compilation. Until that day no module set it, so a participant following module 2 and
running `bnec()` waited about two and a half times as long with nothing on the page to say
why. It is now in the `setup` chunk of modules 2, 4, 5, 6, 7 and 8, and stated visibly in
module 2 beside `mc.cores`. Keep it in the setup chunk of any new module that fits. The
live scripts take it from there rather than setting it themselves, so that a line in a
script is the same text as the line on the page.

**Prompt logging.** Parent §10 applies. Course modules are teaching material about the
analyses, so a change to a module's explanation of a method, to its code, or to which model
is fitted **is** logged. Fixing a typo, restyling, or converting a document's format is
not. Log to `prompts/`, which does not exist yet.

---

## 11. The reference library

The papers behind this course and the `bayesnec` package are in
`C:/Rworking/references/`, outside every repository so that they are never
committed, and reached from here through `ignore/references`. Read one by that
path rather than the absolute one, so the read stays inside the working
directory. `C:/Rworking/bayesnec/notes/references.md` lists what is there and
gives the key in `vignettes/bayesnec.bib` for each paper that package cites.
They are publisher PDFs of copyrighted articles, so they stay local.
