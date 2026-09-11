# Concentration-response modelling training material

Code and worked examples for estimating no-effect toxicity values from
concentration-response data in R, using Bayesian methods via
[`bayesnec`](https://open-aims.github.io/bayesnec/).

**The material is published as a website:
<https://open-aims.github.io/cr_modelling_training/>**

Start there. This repository holds the source; the site is what the material is
written to be read as, and the modules are ordered and cross-linked on it.

## Software setup

Running the code requires R, a C++ compiler, CmdStan and several R packages. The
[software setup
module](https://open-aims.github.io/cr_modelling_training/vignettes/1Getting-started.html)
covers the installation on Windows, macOS and Linux, and ends with a script that
verifies it.

Workshop participants must complete it beforehand. No time is set aside during
the workshop for installing or repairing software.

## Repository layout

| Path | Contents |
|---|---|
| `index.qmd` | the site landing page |
| `vignettes/` | the course modules, and `check_setup.R` |
| `_quarto.yml` | the website definition, including the list of modules to render |
| `docs/` | the built site, served by GitHub Pages |
| `notes/` | decisions and the evidence behind them |
| `ignore/` | local working material, not part of the course |

Modules 2 to 8 are being converted from an earlier version of this course and
are added to `_quarto.yml` as they are completed. A frequentist treatment of the
same material using [`drc`](https://cran.r-project.org/package=drc) is supplied
as reference rather than as part of the taught sequence.

## Building the site

```bash
quarto render
```

The site is rendered locally and the built output in `docs/` is committed.
It is not built by continuous integration, because fitting the Bayesian models
requires a Stan toolchain and takes longer than any runner allows. Computed
chunk output is stored in `_freeze/` and committed for the same reason, so the
site rebuilds without re-fitting.

Render deliberately rather than habitually. Each render rewrites every page, so
a render commit is large; keeping it separate from source changes keeps those
changes reviewable.

## Licence and citation

Released under [CC0 1.0 Universal](LICENSE). The no-significant-effect
concentration is introduced in Fisher and Fox (2023), *Environmental Toxicology
and Chemistry*.

Developed at the [Australian Institute of Marine
Science](https://www.aims.gov.au/). Corrections and questions are welcome as
[issues](https://github.com/open-AIMS/cr_modelling_training/issues).
