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

`vignettes/check_setup.R` was run on Windows 11 (build 26200) with R 4.5.1,
Rtools44, `cmdstanr` 0.8.0, CmdStan 2.36.0, `brms` 2.23.0 and `bayesnec`
2.1.3.0. It halted at stage 5 with "No chains finished successfully" and never
printed its summary. The full account is in
`ignore/windows-checks/windows-setup-report.md`; that folder is not tracked.

Two defects in the script were found and are fixed. `$sample()` does not throw
when every chain fails to start, so stage 5 recorded PASS and then died reading
the result outside its `tryCatch` — meaning a participant in exactly the state
the script exists to detect had nothing to send back. And `cmdstan_model()`
reuses any executable newer than the `.stan` file, so stage 4 reported PASS
without invoking the compiler; the executable on that machine was dated June 2025
and had been built by a toolchain since removed.

**A TBB PATH check was added and then removed, and the diagnosis behind it was
wrong.** The sampling failure was first attributed to CmdStan's TBB directory
being absent from `PATH`. The test supporting that was run in Git Bash, which
carries `C:\Program Files\Git\mingw64\bin` and therefore supplies
`libwinpthread-1.dll` — that DLL, not `tbb.dll`, was the one the stale executable
could not find. `CmdStanRun` places the TBB directory on the PATH of the chain
process itself, and sampling was measured to succeed with it removed from the
user PATH entirely. The check warned on correctly installed machines. A comment
in stage 3 records this so it is not added again.

Two machine faults that do not generalise: a User-scope `RTOOLS44_HOME` pointing
at a directory that does not exist, overriding a correct Machine-scope value, and
the stale executable above. One that does: `C:\rtools44\ucrt64\bin` was empty,
and `check_cmdstan_toolchain(fix = TRUE)` populated it. `cmdstanr` ships
`install_toolchain()` for this, which indicates a stock Rtools does not provide
what it needs. Module 1 now describes that step as required rather than as a
repair. A fresh Rtools install was not tested.

### The Rtools pairing, resolved

The report left open whether module 1's R 4.5.x with Rtools45 had ever been run:
the test machine worked on Rtools44 because `cmdstanr` 0.8.0 mapped every R from
4.4 onward to `RTOOLS44_HOME`.

Resolved by reading `cmdstanr:::rtools4x_version()` in 0.9.0, the version on the
Stan r-universe: it returns "44" for an R minor below 5.0 and "45" otherwise, so
R 4.5.x resolves to `RTOOLS45_HOME`. Module 1's pairing is correct **provided
`cmdstanr` is 0.9.0 or later**. `check_setup.R` now reports the `cmdstanr`
version on Windows, and module 1 states the requirement. This was established by
inspecting the installed 0.9.0 source, not by running it on Windows.

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

| Date | `dev` commit | Version | Modules rendered |
|---|---|---|---|
| 2026-09-11 | `fc33ca84` | 2.1.3.33 | 2, 3 |

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
