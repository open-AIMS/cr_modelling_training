# Revising the course for intermediate practitioners

Written 2026-09-16, before any module was edited.

## Purpose

This note decides what is added to modules 2 to 8 and what is corrected in them.
The course was written in 2023 and converted to Quarto in September 2026. The
conversion preserved the content. Since 2023 `bayesnec` has gained a set of
functions the course does not mention, changed three defaults the course still
states in their old form, and acquired vignettes covering five subjects the
course does not raise at all. Four papers published or revised since 2023 change
what should be taught about normalisation, guideline acceptance of the NSEC, and
experimental design.

The audience is unchanged: practising ecotoxicologists who are competent R users
and are not statisticians. The structure, the module order and the learning
outcomes in `README.md` are unchanged. What changes is the depth at which each
existing subject is treated.

## Sources compared

The comparison covered three bodies of material.

The `bayesnec` `dev` branch vignettes were read in full: `example1` (single
model usage), `example2` (multi-model usage), `example2b` (model details),
`example3` (priors), `example4` (comparing posterior predictions), `example6`
(hurdle and zero-inflated models), `example7` (growth data and negative
responses) and `example9` (a complete analysis workflow). `example8`, the
grouping and factor covariate vignette, was read from the open pull request
`open-AIMS/bayesnec#228`, as was `#372` on the dispersion sub-model. The
`NEWS.md` entry for version 2.2.0 was read for the changes to defaults.

The reference library at `ignore/references` supplied nine papers that bear on
the additions: Fisher et al. (2023, *IEAM*), Fisher and Fox (2023, *ET&C*), Ritz
et al. (2026), Helsel (2006), Weimer et al. (2012), Krull (2020), Mebane (2015),
OECD Test Guideline 201 (2026), and Warne et al. (2025).

The eight course modules were read as they currently stand.

## Corrections required

Three statements in the current material are now wrong. They are separated from
the additions because they must be fixed whether or not anything is added.

### The model weighting default

Module 4 states that stacking is the default weighting method and quotes
`Method: stacking_weights` from a summary. The default is `pseudobma`, and has
been since `bayesnec` issue #320. A set assembled by `c()` or `amend()` was
being weighted by stacking while the same set fitted by `bnec()` was weighted
by pseudo-BMA, and the fix made pseudo-BMA the default everywhere. The package
does not recommend stacking at the sample sizes concentration-response
experiments produce.

### Guideline status of the NSEC

Module 7 states that use of the NSEC is recent and not yet formally adopted.
Warne et al. (2025), the current Australian and New Zealand method, lists the
NSEC alongside the NEC as one of the two most preferred statistical estimates
of toxicity. That document adds a condition the course should teach: an NSEC
whose measured effect size exceeds 20% should probably not be used to derive a
guideline value.

### The ECx reference

Module 3 states that `bayesnec` measures every ECx from the control, which is
correct, but does not say that this changed in version
2.2.0, that `type = "relative"` now means something different from what it meant
before, or that the `hormesis_def` argument has been removed. Anyone returning
to an analysis written against an earlier version needs all three.

## Gaps in the present material

The gaps fall into five groups. Each is a subject the sources treat at length
and the course either omits or states in one sentence.

### Preparation of the data before fitting

The course does not say how to decide the scale of the predictor, and the
`bayesnec` vignette states that this decision generally matters more than the
choice of equation. The course does not mention control normalisation at all.
Ritz et al. (2026) measured what normalising to a control mean does to an ED10:
over 1000 simulations with six control replicates, bias was 6.8% against 2.1%
for the same quantity estimated from the unnormalised response, the coefficient
of variation was 26.4% against 12.7%, and nominal 95% intervals covered at 90%
against 95%. `bnec()` now detects a normalised response and says so, so a
participant will meet this message and the course does not explain it.

### Values that are not measurements

Concentration-response data contain censored values, structural zeros and
values pushed off a distribution boundary, and the course treats none of them.
Helsel (2006) is explicit that substituting a fraction of a detection limit is
fabrication rather than estimation. `bayesnec` passes through the `cens()` term
from `brms`, and `bnec_record()` reports every value the package itself shifted
to bring it inside a family's support.

### The fit as distinct from the sampler

The course checks convergence and stops there. Convergence says whether the
sampler explored the posterior. It says nothing about whether the fitted model
reproduces the variability in the data, and an NSEC depends on that assumption
more than on any other, because its reference is a quantile of the control
posterior. `check_fit()` and a grouped `pp_check()` answer the second question
and neither appears in the course.

### Screening of a model set

Module 4 says a model that has not converged should be removed, and shows no
way of finding one other than inspecting plots. Three functions now do this:
`check_sampling()` reports R-hat, effective sample size and divergent
transitions per candidate, `screen_models()` drops what failed and reports why,
and `failed_models()` names the equations that did not fit at all. The
distinction between a failure of the run, a failure the design cannot fix, and
a model whose shape does not suit the data is the distinction an intermediate
practitioner most needs here, and the sources set it out.

### Reporting

No module states what a methods section has to contain. The `bayesnec` workflow
vignette gives a list, `bnec_record()` supplies most of it
from the fitted object, and `curve_params()` supplies the parameter estimates
that `summary()` does not report.

## Design guidance from the papers

Two findings are worth teaching directly because they change how an experiment
is laid out, and neither appears in the course.

The resolution of the model weights responds more to the number of test
concentrations than to the replication within each concentration (Fisher et al.
2023). A design with twelve concentrations, five replicates and ten trials
gave a weight of 0.987 to the generating model; a design with eight
concentrations and double the within-concentration replication gave 0.984 from a
third more trials.

The NSEC becomes less conservative as the control becomes more variable, because
its reference is a lower quantile of the control (Fisher and Fox 2023). Poor
experiments therefore produce higher no-effect estimates. Mebane (2015) states
the same point for the ECx: at a 20% standard error in the control mean, an EC10
is not a meaningful quantity.

## Order of work

The modules are revised in numerical order, so that a partial completion still
gives a coherent morning. Module 2 is first because the data preparation
material belongs there and three later modules refer back to it.

1. Module 2, which gains the scale of the predictor, control normalisation, the
   seed and sampler settings, and the first sight of `check_sampling()`.
2. Module 3, which gains the ECx vocabulary in full, the censored-estimate
   outcome, `ecnsec()` defined at its first use, and the guideline position.
3. Module 4, which gains the weighting correction, the screening workflow,
   parallel fitting, and the three kinds of failure.
4. Module 5, which gains the dispersion screen, censoring, non-constant
   dispersion, zeros, and boundary substitution.
5. Module 6, which gains `get_priors()`, fixing a parameter, and prior
   sensitivity.
6. Module 7, which gains the screening workflow in place of its own helper,
   and a reporting section.
7. Module 8, which is converted to Quarto and rewritten around `bnec_group()`,
   reading pre-computed results rather than fitting anything.

## The participant's position after the revision

A participant currently leaves the course able to fit a model set, average it,
and extract an estimate. After the revision they also leave able to say why
their predictor is on the scale it is on, to recognise a response that should
not have been normalised, to declare a value below a detection limit rather
than substitute for it, to screen a model set and state which equations were
removed and why, to test whether the fit reproduces the variability their
estimate depends on, and to write the methods section the analysis requires.

## Constraint on new code

New material that needs a fitted model reuses the fit the module already makes,
or is shown with `eval: false`. No module gains a new Stan fit. This was
decided because the site stores computed output under `_freeze/` and a new fit
would have to be sampled on the presenter's machine before the site could be
rebuilt.

The consequence is that some additions describe output rather than displaying
it. Where that happens the text names the function and the argument, so that a
participant running the code sees what the module describes.

The limit is a few minutes per chunk, and it applies to post-processing as well
as to fitting. Extracting **EC~x~** and *NSEC* values across a large model-
averaged set is itself slow: the `bayesnec` grouping vignette takes 38.5
minutes to render against a store in which every model is already fitted, and
none of that time is sampling. Anything above the limit is produced by a
script, saved, and read back by the module. A participant who has to wait ten
minutes for a chunk during a workshop has stopped learning.

Module 7 already works this way. `scripts/generate_herbicide_fits.R` fits the
seven herbicides and writes compact results, and the module reads those. Module
8 adopts the same arrangement.

That covers the published site and those two modules. It does not cover a
participant running modules 2, 4, 5 and 6 alongside the presenter, who fits
everything those modules fit. Rendering the development site on 2026-09-16
measured the wait: 1.5 minutes for module 2, 4.5 for module 4 and 12.9 for
module 5, with module 6 longer again. Thirteen minutes for one module is not a
workshop.

The modules answer it by showing the workflow rather than hiding it. Each fit
appears as three steps, which are the steps an analysis follows:

```r
# shown, not run
set.seed(333)
bnec_fit <- bnec(y ~ crf(x, model = "nec3param"), data = nec_data, seed = 333)
save(bnec_fit, file = "fits/m2_bnec_fit.RData")
```

```r
# run
load("fits/m2_bnec_fit.RData")
```

A first attempt wrapped each call in a `cached()` helper that fitted or loaded as
needed. It worked, and it was the wrong answer: it invented a course-specific
abstraction for something the participant should be learning to do. Fitting
once, saving, and loading is what a real analysis does, and a module that shows
it teaches something beyond the model. The helper was withdrawn.

Two consequences follow. A module no longer samples when the site is built, so a
render takes minutes rather than half an hour. And rebuilding a module from
scratch now needs the saved objects, in the way module 7 already did, which is
what makes the distribution below part of the method rather than a convenience.

## Distributing the fitted models

The archive is a release asset on the course repository, and a USB stick holds
the same file on the day. The two cover different failures: the download is
pre-work and can be done on a connection that works, and the stick covers anyone
whose download failed or who arrives without having done it. Several hundred
megabytes over conference wifi on the morning is not a plan.

The archive is **42.3 MB** over 25 objects, measured on 2026-09-17. That is much
smaller than the 100 to 200 MB estimated before it was built, and small enough
that the download is no obstacle: an end-to-end run of
`source("vignettes/fetch_fits.R")` on that date took 32 seconds. A `bayesnec`
fit of these teaching datasets is a few megabytes at most. The figure was 23.2
MB over eighteen objects when this note was written, and it grows as modules are
revised, so take it from `dist/cr_modelling_fits.txt` rather than from here.

A release asset was chosen over the alternatives because it needs no account to
download, sits beside the material it belongs to, and can be replaced in place
under a tag that does not change, so the address in `fetch_fits.R` stays valid
when the archive is rebuilt. An institutional OneDrive share was used first and
was withdrawn on 2026-09-17: fetched by a client not signed in to the AIMS
tenant, the share link returned 403, so no participant outside AIMS could have
run `fetch_fits.R`. A file-transfer service is unsuitable because its links
expire.

Zenodo remains the destination for a citable identifier once the material is
frozen. `scripts/zenodo_deposit.R` uploads the archive and sets the metadata,
and deliberately stops before publishing, because a published Zenodo record
cannot be withdrawn. The last step is a person looking at the draft. It is not
the route for the workshop, because the archive is still being rebuilt as
modules are revised and a Zenodo record is versioned rather than replaced.

Four pieces support the distribution. `scripts/bundle_fits.R` builds the archive
with a checksum. `vignettes/fetch_fits.R` downloads it, verifies it and unpacks
it, and states that fitting the models instead is fine where it fails.
`check_setup.R` gained a seventh stage that reports whether the fits are
present, as a warning rather than a failure because they are optional. And the
software setup module gained a section covering both routes.

The address is held in one constant at the top of `fetch_fits.R`, and the setup
module names the command rather than the address, so hosting it elsewhere changes
one line and requires no re-rendering of anything that fits a model.

## Shipped datasets

The course holds five data files of its own. Three of them teach something a
dataset shipped with `bayesnec` now covers. Using the shipped data instead
removes undocumented files from a public repository.

One change was made immediately, and then corrected. `manec_example` ships with
the package and fails its sampler screen, so module 4 briefly used it to show a
failing screen. It is a test fixture built at 100 retained draws, so it fails
because it is small rather than because of anything an analysis would meet, and
teaching from it would have taught the wrong lesson. Module 4 now reports the
three outcomes the `bayesnec` workflow article measures on `herbicide` and
`nassarius`, both of which ship: a screen that removes five of eleven equations
at reduced settings, one of ten once the settings are raised, and eight of ten
on a design whose threshold falls in an untested gap.

The larger substitution waits on `open-AIMS/bayesnec#228`, which adds
`coral_colour`, `coral_pam` and `lum31`. Each was built for the grouping
vignette and each matches one of module 8's own files: a chamber that sits at
one concentration, a plate that spans the whole series, and a factor whose
levels are of interest. `lum31` covers all three in one dataset and has
censoring columns besides. Module 5 can use `nassarius` and `herbicide` without
waiting for anything.

The specification gives the mapping file by file, and the three cautions: the
prose has to be rewritten against the new output rather than reused, the new
datasets are larger and so slower to fit, and nothing can be written against
them until the pull request merges.

## The dispersion sub-model

`open-AIMS/bayesnec#372` measures what a dispersion sub-model changes on one
plate of the `lum31` copper series, and module 5 now reports that rather than
arguing the direction in the abstract. Holding the dispersion constant made the
fit simulate about eleven times the spread the data show at the control, and the
resulting *NSEC* was twice the value the sub-model returns. The sampler
diagnostics of that fit were clean throughout, which is the reason this check is
needed: nothing in the usual convergence output says the control is wrong.

The connection to Ritz et al. (2026) makes this a change to the teaching rather
than an addition. A response on an arbitrary scale is often divided by its
observed maximum so that a beta distribution can be used, because a beta
accommodates variability that changes with the mean and a Gaussian does not.
That answers a real problem indirectly and pays for it with the bias Ritz et
al. measured. Modelling the response on its own scale, with the family that
suits it and a dispersion term where the variability changes, addresses the
same problem directly. Module 2 now says so, because dividing by the maximum
and fitting a beta is the practice this group used to follow.

## Prerequisites for building the site

Two modules cannot be rendered as they stand, and this is the most
consequential thing in the note for anyone picking the work up.

`vignettes/data/` does not exist. Module 7 reads its results from there and
`scripts/generate_herbicide_fits.R` has never been run, so module 7 has been
unbuildable since it was converted, which is why it appears in no render list.
Module 8 now joins it, reading from `scripts/generate_grouping_fits.R`. Both
scripts have to be run, deliberately and outside a render, before either module
can be added to `_quarto.yml` or `_quarto-dev.yml`.

Nothing in the project sets the Stan backend. `notes/setup-evidence.md` records
the decision that it is `cmdstanr`, set with
`options(brms.backend = "cmdstanr")`, but there is no `.Rprofile` and no module
sets it, so a plain `quarto render` builds the site under `rstan` instead. That
should be fixed before the next render, either with a project `.Rprofile` or in
each module's setup chunk, so that the published pages are built by the backend
the setup instructions tell participants to install.

`_freeze/` has to be cleared for every module whose chunks changed, which is
modules 2 to 6. Module 4 needs it whether or not its chunks changed, because the
summary output it quotes was produced under the previous weighting default and
`freeze: auto` keys on the source document rather than on the package that
produced the output.

## Upstream changes to watch

`open-AIMS/bayesnec#225` adds `hurdle_poisson` and `hurdle_negbinomial`, the
count analogues of `hurdle_gamma`. A count response with structural zeros is the
common case in this field, a fecundity count where the animal died, so module 5
names the families it can name now and prints the family list rather than fixing
it in prose.

`open-AIMS/bayesnec#373` was raised during this work. `dispersion()` fails with
an internal error on a model-averaged fit rather than returning the per-equation
values or refusing the class by name. Module 5 uses it only on a single fit, so
nothing here depends on the outcome.

## Specification

`notes/intermediate-revision-claude.md` holds the full specification: every gap,
the source that establishes it, the module and section it belongs in, the
evidence behind each claim, and the rejected alternatives.
