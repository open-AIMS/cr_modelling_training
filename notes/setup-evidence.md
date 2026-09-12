# Setup requirements: evidence and decisions

Record of what the pre-workshop setup instructions rest on, and how strong each
part of it is. Written 2026-09-11 while rewriting module 1, which had been
unchanged since 2023 and directed participants to Rtools43 and R 4.3.x.

## Decisions

**Native Windows is the primary route. WSL is the documented fallback.** The
`rstan`-era difficulties on Windows were toolchain mismatches between
`StanHeaders` and `RcppEigen`, hand-edited `Makevars` and `-march=native`. The
`cmdstanr` backend removes most of them. WSL means a second R installation, a
second package library, path translation under `/mnt/c/`, and Positron remote
configuration, which is a second operating system to debug alongside the
modelling. It is the right answer only where compilation is blocked or crippled
by endpoint security, where a Linux server or HPC environment has to be matched,
or where the participant already works in a shell.

**The sampling backend is `cmdstanr`, set once with
`options(brms.backend = "cmdstanr")`.** `rstan` still has to be installed,
because `brms` imports it and loads its namespace, but it does not have to be
the backend.

**Two platforms are tested: native Windows and WSL.** macOS cannot be tested on
the available hardware. WSL doubles as the Linux path and as the Windows
fallback, so testing it serves two purposes.

**The verification script is the contract, not the instructions.**
`vignettes/check_setup.R` defines what "ready" means. Platform instructions are
best-effort routes to passing it. A participant on an untested platform who
passes is ready; one who fails sends the output back before the day. This is
what makes "no troubleshooting during the workshop" achievable, since
instructions alone cannot establish who is ready.

## Verified on Linux, 2026-09-11

`vignettes/check_setup.R` was run end to end under WSL2 Debian and passed every
stage: R 4.6.1, `bayesnec` 2.1.3.7, `brms` 2.23.0, `cmdstanr` 0.9.0, `rstan`
2.32.7, CmdStan 2.39.0 at `~/.cmdstan/cmdstan-2.39.0`,
`check_cmdstan_toolchain()` clean, the bundled `bernoulli` example compiled and
sampled, and a small `brms` model fitted through the `cmdstanr` backend. The
`bayesnec` version check reported the expected warning against the placeholder
target of 2.2.0.

This verifies the end state on Linux. It does not verify the instructions,
because the machine was already configured.

## Observed on a fresh Windows machine

RF assisted an installation on a fresh Windows machine shortly before
2026-09-11. One failure occurred and one fix resolved it.

`cmdstanr::install_cmdstan()` failed. `cmdstanr::install_cmdstan(overwrite =
TRUE, quiet = FALSE)` then succeeded. The likely cause is that the first attempt
created the version directory and failed partway through the build; on a second
attempt without `overwrite`, `cmdstanr` declines to clobber the existing
directory, and the resulting error resembles the original failure without being
it. A transient file lock, of the kind antivirus scanning produces on a
freshly written file, is a less likely alternative.

**Consequence for the instructions:** direct participants to
`install_cmdstan(overwrite = TRUE, quiet = FALSE)` on any retry, and tell them
that `quiet = FALSE` only changes how much build output is printed. The default
output is truncated and does not identify which of the three stages —
download, extraction, build — failed.

## Windows verification, 2026-09-11

Two rounds were run on Windows 11 (build 26200). The full account is in
`ignore/windows-checks/windows-setup-report.md`, which is not tracked; the report
was rewritten after the first round and supersedes itself in three places.

### The specified route passes

The combination module 1 instructs was run in full for the first time — R 4.5.1,
Rtools45, `cmdstanr` 0.9.0, CmdStan 2.39.0 at `C:/cmdstan` — and every stage
passed, with only the `bayesnec` 2.2.0 placeholder warning outstanding. No step
beyond the module's own instructions was needed.

### Two defects in module 1, both of which stop a participant

`install_cmdstan(dir = "C:/cmdstan")` fails when the directory does not exist:
`cmdstanr` does not create it. Module 1 now calls `dir.create()` first.

The CmdStan download exceeds R's default `options(timeout)` of 60 seconds. The
archive is 48.6 MB; on a connection delivering about 500 KB/s the download failed
at the default and completed in 4.9 minutes at `options(timeout = 1800)` with
nothing else changed. The `timeout` argument of `install_cmdstan()` does not
govern the download, and its default of 1200 was in force for the failed attempt.
Conference and hotel networks are commonly slower, so this will affect more
participants than it did there. Module 1 now sets the option and says why.

### Changing Rtools invalidates an existing CmdStan

CmdStan builds its own `tbb.dll`, linked against the toolchain that built it. A
CmdStan 2.36.0 built under Rtools44 requires `libgcc_s_seh-1.dll` and
`libstdc++-6.dll` from `C:\rtools44\ucrt64\bin`, which Rtools45 does not
provide. After installing Rtools45 the model still compiled and every chain
exited with `STATUS_DLL_NOT_FOUND`, reported as "No chains finished
successfully". Rebuilding CmdStan resolved it. Anyone upgrading R and Rtools
between installing CmdStan and the workshop will meet this; it is now in module
1's troubleshooting section.

### Two claims retracted

**The TBB PATH entry is not required.** `CmdStanRun` places that directory on the
PATH of the chain process itself, and the passing run above had no TBB entry in
the user PATH. The original evidence came from a test in Git Bash, which supplies
`libwinpthread-1.dll` — that DLL, not `tbb.dll`, was what the stale executable
could not find. The `tbb_path` check has been removed and a comment in stage 3
records why.

**`check_cmdstan_toolchain(fix = TRUE)` is not a required step.** It was reported
as required because `C:\rtools44\ucrt64\bin` was empty on that machine. That
applies only to `cmdstanr` 0.8.0, or to CmdStan older than 2.35.0. In 0.9.0,
`rtools4x_toolchain_path()` selects `x86_64-w64-mingw32.static.posix`, which the
Rtools installer ships populated, and the check passed against a stock Rtools45.
A stock Rtools does have an empty `ucrt64\bin`; current `cmdstanr` does not use
it. Module 1 has been reverted to describing the step as a check.

### The cmdstanr version requirement

`cmdstanr` 0.9.0 maps R 4.5.x to Rtools45 and refuses Rtools44 outright.
`cmdstanr` 0.8.0 maps every R from 4.4 onward to Rtools44 and never looks for
`RTOOLS45_HOME`, so a participant with an older copy who follows module 1 and
installs Rtools45 is told Rtools44 is missing. Module 1 installs from the Stan
r-universe, which serves 0.9.0, so a clean machine is unaffected and the exposure
is to anyone with an existing installation.

`check_setup.R` reports the `cmdstanr` version on Windows as a **warning below
0.9.0, not a failure**. A warning is the right level: an older `cmdstanr` paired
with the Rtools it expects is a working configuration, as the presenter's machine
shows, so the mismatch is with module 1's instructions rather than with the
machine.

## Advice not yet observed here
## Advice not yet observed here
## The bayesnec version the site is built against

**The site tracks the head of `dev`** (RF, 2026-09-11). The 2.2.0 development
work is close enough to settled that further API changes are expected to be few,
so building against a frozen commit would hold the course behind the package for
no benefit.

```r
remotes::install_github("open-AIMS/bayesnec", ref = "dev")
```

The commit each build actually used is recorded here rather than pinned, so a
published page remains attributable to a state of the package.

| Date | commit | Version | Modules rendered |
|---|---|---|---|
| 2026-09-11 | `dev` `fc33ca84` | 2.1.3.33 | 2, 3, 5 |
| 2026-09-11 | PR #321 `3226749` | 2.1.3.34 | 4, 6 |

The second row departs from `dev`. Module 4 demonstrates model averaging by
combining single fits with `c()` and `amend()`, which is the path on which the
pseudo-BMA weighting default was lost (#320); PR #321 supplies the default
wherever a set is assembled, and RF approved building module 4 against whatever
version works. Module 6 was rendered in the same pass.

Modules 2, 3 and 5 remain frozen against `fc33ca84` and are unaffected: #321
changes only how weights are defaulted when a set is assembled outside `bnec()`,
and none of them does that. Re-render everything once #321 merges to `dev`, so
the whole site reports one version.

Re-check `dev` before a render, reinstall if it has moved, and clear `_freeze/`
for every module that fits a model, because freeze does not notice a package
upgrade. Add a row above when the commit changes.

Module 1 is a separate matter. Participants are told to install a specific
version, and that instruction cannot track a moving branch; its placeholder is
filled from the release once the 2.2.0 work concludes.

Two 2.2.0 changes were confirmed present at `fc33ca84` by rendering module 2
against it: the default `resolution` in `ecx()` and `nsec()` is 200 rather than
1000, and `hormesis_def` is absent from both while `ecnsec()` exists. The
`ne_posterior` and `ne_type` elements of a `bayesnecfit`, and the `ecx()` default
of EC10 and `nsec()` default `sig_val` of 0.01, are unchanged from 2.1.3.7.

## Instruction corrections for module 1

- Rtools must match R: **Rtools45 for R 4.5.x**, Rtools44 for 4.4.x. Module 1
  currently names Rtools43 and R 4.3.x and links to
  `rtools43-5550-5548.exe`. A mismatched Rtools is the most common single
  failure.
- `cmdstanr` is installed from **`https://stan-dev.r-universe.dev`**. Module 1
  uses `https://mc-stan.org/r-packages/` and states that `cmdstanr` "is not on
  CRAN yet".
- Module 1 describes CmdStan with reference to version 2.29.2 and includes a
  build-from-git-clone route. `cmdstanr::install_cmdstan()` replaces all of it.
- Add `options(cmdstanr_cmdstan_path = "C:/cmdstan/cmdstan-x.xx.x")` to
  `.Rprofile` where CmdStan was installed to a non-default path, or every new
  session reports CmdStan as absent.
- Give participants a time box. Thirty minutes for the native route, after
  which the fallback is WSL rather than further configuration.

## Provenance

The Windows material comes from a conversation in which RF assisted the
installation. The observed failure and its fix are first-hand. The four
conditions in "Advice not yet observed here" are from the same conversation but
are general guidance rather than observations from that machine, and are
labelled as such above.

## Module 7 is on hold

**RF, 2026-09-12.** Conversion is written and committed but is not registered on
the site, and the fit-generation run was stopped after a few minutes.

The reason is not the fits. Module 7 fits seven herbicides by splitting the data
and mapping `bnec()` over the pieces, which is what the 2023 analysis did and
what `scripts/generate_herbicide_fits.R` reproduces. `bayesnec` now exports
`bnec_group(formula, data, group_var, family, ...)`, which fits the levels of a
grouping variable in one call and returns one object, together with
`crossed_group_weights`, `compare_estimates` and `compare_fitted`. The
split-and-map pattern is superseded, and teaching it would teach the old way of
working.

The new grouping vignette on PR #228 is the reference for how this should now be
written. Module 7 should be rewritten around `bnec_group()` once that vignette
has settled, and the generation script rewritten with it or dropped.

Module 8 is the grouping and factor-covariate module and is the same PR's
subject, so 7 and 8 should be planned together rather than separately. Module 8
should also no longer need the raw `brms` helpers in `vignettes/functions.R`,
which were written before this functionality existed.
