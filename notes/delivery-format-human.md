# Delivery format for the one-day workshop

Written 2026-09-16, before any delivery material was built.

## Purpose

This note decides what is displayed on the projector, what each participant has
open on their own machine, and how code is run in a room where a single model
fit takes minutes. It does not decide the content of any module; that is decided
in `intermediate-revision-human.md`.

It is written before that revision finishes because the answer changes what a
module has to contain. A module taught from its own page needs its code
runnable in the room, and that is a property of the module source rather than
of the presentation.

Two terms are used throughout and are defined here. The **fits bundle** is the
set of saved model objects written to `fits/` by
`scripts/generate_taught_fits.R`, which a participant downloads so that nothing
has to be sampled on the day. A **live script** is the plain `.R` file holding
one module's executable code, extracted from the module source, which the
presenter types from and the participants run.

## The projected material

The projector shows the rendered module page. No parallel slide deck is built
for modules 1 to 8.

The reasons, in order of weight.

A deck and a page covering the same material diverge as soon as either is
edited, and the revision in `intermediate-revision-human.md` will edit every
module. Section 13 of the working-directory `CLAUDE.md` makes this argument for
planning documents, and it holds here for the same reason: a decision stated in
two places is a decision that will shortly be stated two ways.

The participant has the same page open on their own screen, so an instruction
to look at the section headed *The model weights* resolves for everyone. A
participant who looks away for two minutes can rejoin a page and cannot rejoin
a deck.

Chunk output is frozen under `_freeze/`, so the figure on the screen is the one
their own script reproduces. A redrawn version of the same figure on a slide
is a second thing that can disagree with the code.

Eight modules is roughly two hundred slides of new material to build before
September, on top of a revision that has two modules still unbuilt and three
corrections still to apply.

## The participant's copy of the code

Each taught module gets a live script at `scripts/live/<module>.R`, generated
from the module source rather than written by hand. The presenter runs it in
Positron beside the browser, and every participant has the identical file open.

The generation already exists in part. `scripts/generate_taught_fits.R` holds
`executable_chunks()`, which pulls every `{r}` block whose options do not say
`eval: false` out of a `.qmd`. That function is moved to a shared file and used
for both purposes, so a live script cannot drift from the module it came from.
Specification section *Live script generation*.

Generating the scripts by hand is rejected for the same reason the deck is
rejected.

## Fits loaded rather than sampled

Every fit in a taught module is shown as three steps: the `bnec()` call, the
`save()` that keeps the result, and the `load()` that reads it back. The first
two are displayed and not run; the third is what the render and the participant
execute. A participant holding the fits bundle has the fit in seconds, and a
participant without it runs the `bnec()` call instead and waits, so the material
still works for someone reading it at home.

The measured fit times are the reason. Rendering the development site on a
four-core WSL2 machine on 2026-09-16 took 1.5 minutes for module 2, 4.5 minutes
for module 4 and 12.9 minutes for module 5. A room that
waits thirteen minutes for one module has stopped learning, and there is no
version of the day in which those fits are sampled live.

A `cached()` helper wrapping each fit call was written first and then withdrawn
on 2026-09-17. It worked. It also invented a course-specific abstraction for the
thing a participant should be learning to do, which is to fit once, save the
object, and load it next time. `vignettes/fit_cache.R` does not exist, and the
specification sections describing it are superseded; the reasoning is in
`intermediate-revision-human.md`, section *Constraint on new code*.

## The boundary between running and loading

Each taught module opens with a **runs-here callout**, a note naming the code a
participant runs and the code they watch. The purpose is to stop a participant
starting a fit that will not finish before the module ends, which is the failure
that takes a person out of the room for the rest of a session.

The boundary is stated on the page rather than only said aloud, because someone
working through the material alone afterwards needs it as much as someone in
the room. Specification section *The runs-here callout* gives the wording and
the per-module content.

## Slides

One deck of ten to fifteen slides is built, covering the opening of the day
only: the presenter, the shape of the day, and what a concentration-response
model is before any code appears. It is not extended to cover module content.

Conceptual figures go into the module pages rather than into slides. The ones
worth drawing are what NEC, NSEC and ECx each measure on a fitted curve, what
the sampler does, and what a model weight means. All three belong in the
permanent material, because that is where someone will look for them in six
months.

## The fits bundle and its distribution

The distribution of the bundle is settled before the course rather than during
it. Under the format decided here a participant without the bundle cannot run
any of the code, so it is no longer an optional download.

Two routes are used together. The bundle is a release asset on the course
repository, at
<https://github.com/open-AIMS/cr_modelling_training/releases/download/fits/cr_modelling_fits.zip>,
which `vignettes/fetch_fits.R` downloads; USB sticks hold the same file for a
venue where several dozen people download at once. The tag `fits` does not
change, so rebuilding the archive replaces the asset and the address in
`fetch_fits.R` stays valid.

The archive is 42.3 MB over 25 objects, measured on 2026-09-17, which is well
inside the 2 GB limit on a release asset and downloads in about half a minute.

An institutional OneDrive share was used first and does not work. Fetched on
2026-09-17 by a client not signed in to the AIMS tenant, the share link
redirected to the file and then returned 403, so no participant outside AIMS
could have run `fetch_fits.R`. A share link is served through a sign-in, and
that cannot be tested by opening it in the browser of someone who is already
signed in. Specification section *Bundle distribution*.

The same script must be re-run after any upgrade to `bayesnec`, `brms` or Stan,
for the reason `CLAUDE.md` section 10 gives for `_freeze/`: a saved result does
not know that the package which produced it has changed.

## The time budget

The material does not fit a day at its present length, and no choice of
presentation format changes that. A decision is required about what is taught in
full, what is demonstrated briefly and what becomes reading. It is the one
decision in this note that cannot be deferred to the specification.

Measured on 2026-09-16, modules 1 to 8 hold 187 code chunks, of which 32 are
marked `eval: false`, leaving 155 that execute. Assuming a nine-to-five day
with an hour for lunch and two twenty-minute breaks, there are 380 minutes of
teaching. If the whole day went to code and nothing else, that is 2.45 minutes
per chunk. Half the day will go to explanation, which leaves a little over a
minute per chunk. That is less than the time to introduce a call, run it, and
let several dozen people run it.

The source word counts are 4,938 for module 2, 4,470 for module 4 and 5,406 for
module 5, against 1,859 for module 1. These were taken with `wc -w`, which
counts backticked identifiers and table syntax, so they overstate the prose;
they are given as a relative measure of module length and not as a reading
time.

The triage is not made here because it needs the presenter's judgement about
which outcomes in `README.md` matter most. Specification section *Material that
does not fit* records the options and what each one gives up.

## The room

Practical points to settle before the day, each of which is difficult to fix
once a session has started.

The venue network may not support several dozen simultaneous downloads, so the
built site in `docs/` is shipped as a zip on the same USB as the fits bundle
and opens from the file system without a server.

The docked sidebar takes width that a projected code chunk needs. Collapse it
when presenting, or set `page-layout: full` for the session.

Code legibility at the back of a room is a property of the browser zoom rather
than of the theme. Check a code chunk from the back row before the day starts.

The site is published in two builds, the approved one at the repository root
and the development one at `/dev/`. Participants are given one URL, and that
must be the approved build with every taught module in it. The
`_quarto-dev.yml` split is a development convenience and should not survive to
the day.

## Order of work

Each step is complete when the stated condition holds.

1. Show each fit in the taught modules as its call, its `save()` and its
   `load()`. Done when a render with `fits/` populated samples nothing.
2. Move `executable_chunks()` to a shared file and add the live script
   generator. Done when `scripts/live/` holds one script per taught module and
   each runs start to finish against the bundle.
3. Generate the bundle and record its size. Done when `fits/PROVENANCE.txt`
   exists and the size decides the distribution route.
4. Make the triage decision on the time budget. Done when every module has a
   stated number of minutes and the total is under 380.
5. Add the runs-here callout to each taught module, with content that follows
   from step 4. Done when each taught module opens with one.
6. Build the opening deck.

Steps 1 and 2 are independent of the revision and can proceed alongside it.
Step 5 depends on step 4 and on the module content being settled.

State on 2026-09-17: steps 1, 3, 4 and 6 hold. Every taught module shows call,
save and load, a render of the development site samples nothing, the bundle is
built, measured and published as a release asset, `workshop-agenda.md` gives each
module its minutes, and the opening deck is drafted at `slides/opening.qmd`.
Steps 2 and 5 are outstanding: `scripts/live/` does not exist, and no module
opens with a runs-here callout.

The deck is a `revealjs` document with `embed-resources: true`, so it renders to
one self-contained `slides/opening.html` that opens from a USB stick with no
other files beside it. It is not in either `render:` list, so a project render
neither builds nor deletes it; build it with `quarto render slides/opening.qmd`.
The rendered file is git-ignored.

## Rejected routes

`quarto-live` and webR cannot compile Stan, so the Bayesian half of the course
would become non-executable code display. This was already settled in
`CLAUDE.md` section 2 and is recorded here so it is not reopened.

A `revealjs` format target on the same `.qmd` files was considered and
rejected. Long-form teaching prose does not partition into slide-sized sections
without being rewritten, and rewriting it for the second target returns the
maintenance problem the single-source decision was made to avoid.

The retired `learnr` set and its shinyapps.io deployment are not revived.

## The live scripts as built

Step 2 of the order of work is done, and this section records what was built
and what it was measured at. The parsing that both generators share is in
`scripts/qmd_chunks.R`, and `scripts/generate_live_scripts.R` writes
`scripts/live/`, holding one `.R` file per taught module and a README. The
shared parser was checked against the function it replaced in
`generate_taught_fits.R`: for all eight modules the extracted code is
identical, so the fits bundle does not need rebuilding.

Each script was sourced start to finish in a fresh R session on 2026-09-19,
under R 4.6.1 with `bayesnec` 2.1.3.39, `brms` 2.23.0 and `cmdstanr` 0.9.0,
against the 26 objects in `vignettes/fits/`. All eight completed without error.
The elapsed times were 0, 91, 6, 246, 248, 56, 416 and 4 seconds for modules 1
to 8 in order, the ten live fits included. The 26 saved objects were checksummed
before and after and none changed, which is what the commented `save()` lines
are there to ensure.

Module 2 is the check on the estimates. It ran in 91 seconds against a
prediction of 94, being the 76 seconds estimated for its two fits on the 18
seconds the script took before any fit was run live. Module 7 spends most of
its 416 seconds on the posterior predictive checks and the model weights, which
are computed from the loaded objects rather than read from them, and module 8
takes 4 seconds because its results are read from `data/` as comma-separated
files.

Those times are for a whole module sourced in one go. A demonstration block
runs a few lines at a time, so the figures above are a ceiling rather than a
wait anyone will sit through.

A chunk is commented out in the script when running it would take too long or
would fail, and `eval: false` on the page is neither of those. The page
withholds a chunk for whatever suits the render. Most of those chunks are fits,
and a few are a single `options()` call that returns instantly and that a
participant should run. The first version of the generator commented out all 57
of them, which teaches the wrong thing about `options(mc.cores = 4)` and spends
on it the one mechanism available for saying *do not run this*.

### The time a fit takes

A fit is commented out where it would take more than `MAX_LIVE_SECONDS` in the
generator, which is 50. Estimating that needs two quantities, because neither
alone is the wall time a participant waits.

The sampling time is in the saved object, and `scripts/measure_fit_times.R`
reads it out with `rstan::get_elapsed_time` into `scripts/fit_times.csv`. A fit
takes as long as its slowest chain, the chains being run at the same time, and
a model set takes the sum over its models, `bnec()` fitting them in turn. Across
the 26 saved objects the sampling times run from 0.8 seconds to 241 seconds,
and 21 of the 26 sample in under 15 seconds.

Compilation is the rest, and is the larger part of a small fit. It is not
recorded in the object, so it was measured directly on 2026-09-19 on the
four-core WSL2 machine: a `nec3param` fit on `nec_data` took 38.5 seconds of
wall time against 2.9 seconds of sampling, and an `ecxll3` fit took 33.5
seconds against 2.4. Refitting a model already compiled in the same session took
6.6 seconds, so the charge falls once per distinct equation rather than once per
call. `COMPILE_SECONDS` is set to 35 from those three measurements. A
participant's machine is not that machine, and Windows with Rtools is generally
slower, so the estimates are a floor.

The estimate is therefore 35 seconds for each distinct equation plus the
measured sampling. Module 2's first fit comes out at 38 seconds, which is what
the module tells the participant when it says the fit took a little over two
minutes on a slower machine. An `amend()` call is charged for the equations
named in its `add` argument rather than for every equation in the result, since
it refits only what it adds.

A fit with no saved object cannot be estimated and stays commented. That covers
the calls made only for illustration, such as module 4's `future` plans, and
every fit in module 8, which are made on the cluster and never saved here.

### The threshold and its effect on each module

The threshold was set at 120 first and lowered to 50 on 2026-09-19 (RF). Two
minutes is the right rule for a single fit and the wrong one for a module. At
120 the nine fits of module 5 ran, about 10 minutes of a block the agenda gives
30 minutes as a demonstration, and module 5 is about recognising which
distribution suits an endpoint rather than about executing a particular call.

At 50, ten fits run and twenty-eight stay commented. The time this adds is 76
seconds for module 2, about 2 minutes for module 4, 72 seconds for module 5, 83
seconds for module 6 and 40 seconds for module 7, so a little over five minutes
across the day. Module 5 keeps its two single-equation fits and its seven
two-equation sets are commented.

The value is 50 rather than 45 because the estimates are not evenly spread. The
ten fits that run estimate between 36 and 45.1 seconds and the next is 73, so
any threshold from 46 to 72 gives the same ten. A threshold of exactly 45 would
turn on the 0.1 second by which module 6's `fit-a` exceeds it, which is well
inside the error of an estimate built on a compile time measured three times.

A per-module budget was considered as the alternative to lowering the per-fit
threshold, and is not implemented: it adds a second constant to reason about,
and which fits it drops depends on their order in the module rather than on
what they teach.

### The Stan backend

A script that fits sets `options(brms.backend = "cmdstanr")` in its preamble,
guarded on `cmdstanr` being installed. `brms` uses `rstan` unless told
otherwise, and `rstan` compiles a model far more slowly.

This was found by running the scripts rather than by reading them. Module 2's
two fits took 199 seconds against an estimate of 76, and the log held
`SAMPLING FOR MODEL 'anon_model'`, which is `rstan` output. `COMPILE_SECONDS`
had been measured in a session with the backend set, so the estimate was right
for `cmdstanr` and wrong for the backend the script used.

The guard matters because a participant who has not finished the software setup
has no `cmdstanr`, and setting the option unconditionally would leave every fit
failing rather than merely slow.

The same gap is in the modules. Module 1 tells a participant to set the backend
and no module sets it, so someone who runs `bnec()` while following module 2
compiles under `rstan` and waits about two and a half times as long, with
nothing on the page to say why. The live scripts no longer have that problem and
the modules still do.

### The save() inside a fit that runs

The `save()` line inside a fit chunk is commented out wherever the fit itself
runs. It writes into `vignettes/fits/`, which holds the distributed objects, so
running the script would replace them with fits made under whatever versions the
participant has and the checksum in `vignettes/fits.sha256` would no longer
match what they downloaded. The fit still happens and the object is still in the
session. The `load()` below it then restores the distributed copy, with a line
saying so, because otherwise it reads as the participant's own fit being
discarded for no reason.

### Objects an earlier module fitted

A chunk that refers to objects the script never creates raises an error and is
commented, with the comment naming them: module 2's `ecx(fit, xform = ...)`
refers to `fit`, where the script has `bnec_fit`. Module 3 is the exception.
Its `ecnsec()` demonstration uses `bnec_fit` and module 3 never loads it,
because the page reads on from module 2 where it was fitted. The generator
resolves a name like that against the objects the modules save, and emits the
`load()` that makes the block run, marked as an addition. The file chosen is the
one saved by the latest module at or before this one, so module 3 gets module
2's fit rather than module 6's, which is a different fit under the same name.

A chunk whose body is a bare formula is commented as showing the shape of a
term, its names being placeholders.

### The free-variable check

Whether a chunk refers to something absent is decided by
`codetools::findGlobals` against the names assigned by the chunks emitted before
it, rather than by `all.vars`, which reports a function's own arguments and its
loop variables as references. Three corrections were needed before that was
reliable.

Formulas are removed before the free variables are counted, because a formula is
not evaluated when the call around it is and its names are not references the
script has to satisfy. Module 8's `brms::bf(y ~ bot + ..., nl = TRUE)` runs
against an empty session, and counting `bot`, `top` and the rest would have
commented out a call that cannot fail. The substitution puts `TRUE` in the
formula's place rather than `NULL`, because assigning `NULL` into a call removes
that element and shrinks the call underneath the loop walking it.

The target of a replacement assignment is counted as a read, because module 6's
`fixed_prior$prior[i] <- "beta(6, 6)"` needs `fixed_prior` to exist and
`findGlobals` reports it as a local assignment.

The map from a saved file to the names it restores is built across all eight
modules rather than per module, because module 3 loads module 2's object. It is
read from the `save()` call in the source rather than from the file, so
generation works without the bundle on disk.

The commented blocks are commented rather than wrapped in `if (FALSE) { ... }`.
The wrapper was written first and rejected. It keeps the syntax highlighting,
and in Positron and RStudio a cursor on a line inside the block still sends
that line to the console, so a participant working down the file line by line
starts the fit anyway. That is the failure the runs-here callout exists to
prevent. The projector shows the rendered page, which is where the call is read
from, so the script gains nothing from being the prettier of the two.

A chunk whose body is nothing but `knitr::include_graphics()` is dropped. It
places a figure on the page and does nothing in a console, and no later chunk
depends on it. Thirty-two were dropped across the eight modules, eleven of them
in module 3.

Each script sets the working directory to `vignettes/` before anything else.
Every path in a module resolves from there, and setting the directory keeps a
line in the script identical to the line on the page. Rewriting the paths
instead would give a participant comparing the two a difference to explain.

## Specification

`delivery-format-claude.md` holds the file layout, the generator, the callout
wording, the per-module inventory of fit calls, and the evidence behind each
decision above.
