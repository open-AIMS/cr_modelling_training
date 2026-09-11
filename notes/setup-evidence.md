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
`setup/check_setup.R` defines what "ready" means. Platform instructions are
best-effort routes to passing it. A participant on an untested platform who
passes is ready; one who fails sends the output back before the day. This is
what makes "no troubleshooting during the workshop" achievable, since
instructions alone cannot establish who is ready.

## Verified on Linux, 2026-09-11

`setup/check_setup.R` was run end to end under WSL2 Debian and passed every
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

## Advice not yet observed here

The following are the recurring Windows failure modes. None was encountered on
the fresh machine above, so they are recorded as anticipated rather than
measured, and the instructions should not imply otherwise.

| Condition | Effect | Response |
|---|---|---|
| `HOME` inside OneDrive | CmdStan build fails; the default target is `~/.cmdstan` | `install_cmdstan(dir = "C:/cmdstan")` |
| Space or non-ASCII character in the username | build fails | as above |
| Windows Defender real-time scanning | each compile is inspected; a 40-second compile can take several minutes | exclude the CmdStan directory, if permitted |
| Policy blocking freshly compiled executables | compilation cannot proceed | WSL; this is the case WSL exists for |

`setup/check_setup.R` reports on the first two directly, by inspecting `HOME`
before anything is compiled.

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
