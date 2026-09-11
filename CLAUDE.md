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

The course is currently eight `learnr::tutorial` documents with `runtime:
shiny_prerendered`, written in R Markdown and deployed to shinyapps.io under the
`open-aims` account (`vignettes/rsconnect/`).

**Decision (RF, 2026-09-11): the set migrates to Quarto before the course.** The parent
§9 prohibition on R Markdown therefore applies here without an exception, and the learnr
set is retired rather than extended. Do not author a new `.Rmd` in this repo.

**Every module stays in `vignettes/`.** The course is a series of vignettes and that is how
collaborators expect to find it. A converted module replaces its `.Rmd` in place, keeping
the same file stem, so the numbering and any existing reference still resolve. Do not move
a module into a new directory; module 1 was briefly split into a `setup/` folder on
2026-09-11 and moved back the same day.

Module 1 is converted. It is `vignettes/1Getting-started.qmd`, with
`vignettes/check_setup.R` beside it, and it is the only module that is not taught: it is
issued as pre-work, and a participant who has not completed it cannot fit models during
the workshop. There is no troubleshooting slot on the day.

**The runtime is not yet decided, and it is the first thing to settle.** The three routes
differ in what they can execute, and the course fits Bayesian models:

| Route | Executes `drc` | Executes `bayesnec`/Stan | Needs a server |
|---|---|---|---|
| `quarto-live` (webR) | expected yes, untested here | **no** | no |
| Quarto with `server: shiny` | yes | yes | yes, as now |
| Static Quarto, participants run code locally | n/a, participant's machine | yes | no |

webR has no C++ toolchain, so `rstan` and `cmdstanr` cannot compile a model in the
browser and `brms` cannot run there. This has not been tested in this repository; it is
stated from how webR and Stan work, and it should be confirmed before the route is
chosen. A `quarto-live` conversion would leave the Bayesian half of the course as
non-executable code display, which removes the exercises from modules 2, 4, 6, 7 and 8.

Whichever route is taken, the migration converts **all** modules. A part-migrated set —
some modules learnr, some Quarto — is not acceptable for a workshop, because participants
would meet two interaction models in one day.

Record the outcome in this section when it is settled, and note the deployment target.

**Scale of the work.** Eight modules, roughly 170 KB of source, four of them depending on
pre-fitted objects (§5). Converting all of it before the course is a substantial task.
Sequence it by module number so that a partial completion still gives a coherent morning.

---

## 3. The module set

Numbered modules are the Bayesian thread; a `d` suffix is the `drc` counterpart of the
same subject, and the two are taught back to back.

| # | File stem | Title | Notes |
|---|---|---|---|
| 0 | `0Overview` | course outline | title is wrong, see §7 |
| 1 | `1Getting-started` | Software setup | **converted to `.qmd`**, pre-workshop, not taught |
| 2 | `2Fitting-a-CR-model-using-bayesnec` | Fitting a single model using bayesnec | |
| 2d | `2d Fitting-a-CR-model-using-drc` | Fitting a single model using drc | **space in the filename** |
| 3 | `3Toxicity_estimation_and_available_models` | Toxicity estimation and the available models | |
| 3d | `3dToxicity_estimation_and_available_models` | Toxicity estimation and the available models in drc | |
| 4 | `4Model_averaging_and_multimodel_inference` | Model averaging and multimodel inference | |
| 4d | `4dModel_averaging_and_multimodel_inference` | Model averaging in drc | **not learnr**, see §7 |
| 5 | `5Response_data_and_statistical_distributions` | Modelling your response using the right statistical family | |
| 6 | `6Priors_and_Bayesian_inference` | Priors and Bayesian inference | |
| 7 | `7Example_case_study` | Worked example and comparing toxicity | |
| 7d | `7dExample_case_study` | — | **source deleted**, recoverable, see §7 |
| 8 | `8Factor_covariates_and_groupings` | Factor covariates and gruopings | title typo |
| 8d | `8dFactor_covariates_and_groupings` | Are these ECx values different? | |

`vignettes/Scratch.Rmd` and `vignettes/test_4dModel_averaging_and_multimodel_inference.Rmd`
are working files, not course modules. `scratch/` and `ignore/` are likewise not course
content; `ignore/` is git-ignored in full.

`vignettes/functions.R` holds three helpers written for module 8 before the equivalent
functionality existed in the packages: `pred_out()`, `nec.brmsfit()` and `nsec.brmsfit()`,
all for `brmsfit` objects with a grouping variable. `nsec.brmsfit()` calls three
unexported `bayesnec` internals (`bayesnec:::do_wrapper`, `bayesnec:::modify_posterior`,
`bayesnec:::min_abs`). Check against the current `bayesnec` before relying on it; triple-colon
calls break without deprecation.

---

## 4. Packages in scope

No `DESCRIPTION` and no `packages.R`, so parent §4 gives no canonical list. The list below
was taken from the `library()` calls and `::` usage across `vignettes/*.Rmd` on 2026-09-11
and is the working set until a `packages.R` is created:

`bayesnec`, `drc`, `brms`, `cmdstanr`, `rstan`, `posterior`, `bayesplot`, `learnr`,
`knitr`, `tidyverse` (`ggplot2`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`), `scales`,
`ggpubr`, `cowplot`, `car`, `qwraps2`, `parallel`.

**Create `packages.R` as part of the migration.** A full-day workshop where participants
install packages needs one authoritative list, and eight scattered `library()` blocks is
not one. Ask before adding anything to it.

**`extractNSEC` no longer exists.** `2d Fitting-a-CR-model-using-drc.Rmd:17` calls
`library(extractNSEC)`, and the last commit on the repository (`4298c3d`, "use new
package", 2023-11-07) added `nsec()` calls for `drc` objects on the strength of it. That
functionality is now in **`toxval`**, which exports `nsec()` with an `nsec.drc` method
(confirmed in `C:/Rworking/toxval/NAMESPACE`). The module must be repointed, and the
commented-out `abline()` at line 137 of the same file resolved rather than left commented.

`bayesnec` is installed from the `dev` branch by module 2
(`remotes::install_github("open-aims/bayesnec", ref = "dev")`). Pin this to a release or a
commit before the workshop. A `dev` branch that changes during the course produces
failures that cannot be diagnosed in the room.

---

## 5. Data files and their availability

`.gitignore` excludes `*.RData`, so three large files that the material depends on are
present only on the presenter's machine:

| File | Size | Used by |
|---|---|---|
| `vignettes/fitted_model.RData` | 124 MB | modules 2, 4, 6 |
| `vignettes/manfecfit.RData` | 451 MB | module 7 |
| `vignettes/ametryn.RData` | 65 MB | modules 4, 7 |

These stay in `vignettes/`, because the modules `load()` them by relative path at render
time.

**A fresh clone cannot render modules 2, 4, 6 or 7.** Anyone rebuilding the course from
GitHub (a co-presenter, a participant, a future session) is blocked at this point. Settle
the distribution route as part of the migration: regenerate the fits from a script held in
the repository, or publish the objects as a release asset, or reduce them to the posterior
draws the modules actually use. Do not commit them; the repository is public and 640 MB of
`.RData` is not a reasonable clone.

The three CSV files are tracked and small: `example_ogl.csv`, `example_pgl.csv` and
`example_fi.csv`, read by module 8.

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
publisher. The repository is public. Check each before the course, and prefer redrawing a
figure from its underlying quantities over reproducing a publisher's rendering.

---

## 6. Prose and register

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
been corrected in this working directory before, in the `bayesnec` vignettes (August 2026)
and in the `ndorsatus` reanalysis report (2026-08-14).

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
or cite the help page or the paper. Parent §11 records six traps in `bayesnec` and `drc`
that are counter-intuitive and easy to state wrongly — read it before writing about links,
transformed predictors, `dispersion()`, `drc` standard errors or `type = "Poisson"`.

**Name model equations precisely.** Parent §12 "Naming models precisely" applies with
particular force here: this course is where participants learn the vocabulary. Never write
"the nec model" for a specific equation; write `nec3param`, `nec4param`, `nechorme`. Keep
the three senses of `ecx` distinct — the equation group, the `ecx()` function, and an ECx
estimate.

**Exercise and solution text.** An exercise prompt states the task and nothing else. A
solution states what the code does and why that answer, not *Great — you have now
mastered…*. No congratulation, no encouragement written in advance of the participant
doing anything.

---

## 7. Known defects in the current material

Found by inspection on 2026-09-11 and not yet fixed. These are recorded so that the
migration does not reproduce them. The revision plan itself is decided separately.

- **`css/style.css` does not exist.** Every learnr module's YAML names it. The whole set
  is running unstyled, or partly so.
- **`0Overview.Rmd` has module 1's title**, "Getting started - Installing and running
  BRMS". Its content is the course outline.
- **`4dModel_averaging_and_multimodel_inference.Rmd` is not a learnr tutorial.** It is
  `html_document`/`pdf_document`, so it behaves differently from every other module.
  `test_4dModel_averaging_and_multimodel_inference.Rmd` appears to be a learnr version of
  the same material. Decide which is the module.
- **`7dExample_case_study.Rmd` was deleted from the repository.** The rendered
  `7dExample_case_study.html` is still present. The source was added in `0f53a3d`, last
  modified in `47abfcd`, and deleted in `8d55cd2` ("Minor updates") with no explanation.
  Recover it with `git show 47abfcd:vignettes/7dExample_case_study.Rmd` and establish
  whether the deletion was deliberate before reinstating the module.
- **Title typo:** "gruopings" in `8Factor_covariates_and_groupings.Rmd`.
- **README typo:** "the backage bayesnec".
- **`extractNSEC` dependency is dead** — §4.
- **`bayesnec` is installed from `dev`** — §4.
- The material was last touched in November 2023. `bayesnec`, `brms` and `drc` have all released
  new versions since. Assume nothing runs until it has been run.

---

## 8. Working in this repository

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

**The site is built against the head of `bayesnec`'s `dev` branch.** Check `dev` before
rendering, reinstall if it has moved, and record the commit used in the table in
`notes/setup-evidence.md`. The commit is recorded rather than pinned, so a published page
stays attributable without holding the course behind the package. Module 1 is the exception:
it tells participants which version to install, so it names a release rather than a branch.

**Prompt logging.** Parent §10 applies. Course modules are teaching material about the
analyses, so a change to a module's explanation of a method, to its code, or to which model
is fitted **is** logged. Fixing a typo, restyling, or converting a document's format is
not. Log to `prompts/`, which does not exist yet.
