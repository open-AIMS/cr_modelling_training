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

State on 2026-09-17: steps 1 and 3 hold. Every taught module shows call, save
and load, a render of the development site samples nothing, and the bundle is
built, measured and published as a release asset. Steps 2, 4, 5 and 6 are
outstanding, and `scripts/live/` does not exist.

## Rejected routes

`quarto-live` and webR cannot compile Stan, so the Bayesian half of the course
would become non-executable code display. This was already settled in
`CLAUDE.md` section 2 and is recorded here so it is not reopened.

A `revealjs` format target on the same `.qmd` files was considered and
rejected. Long-form teaching prose does not partition into slide-sized sections
without being rewritten, and rewriting it for the second target returns the
maintenance problem the single-source decision was made to avoid.

The retired `learnr` set and its shinyapps.io deployment are not revived.

## Specification

`delivery-format-claude.md` holds the file layout, the generator, the callout
wording, the per-module inventory of fit calls, and the evidence behind each
decision above.
