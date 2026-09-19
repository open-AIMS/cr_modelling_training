# Live scripts

One `.R` file per taught module, holding the code from that module's page.
Open the file beside the page and run it as the module is worked through.

| Script | Module |
|---|---|
| `1Software-stack.R` | The software stack |
| `2Fitting-a-CR-model-using-bayesnec.R` | Fitting a single model |
| `3Toxicity_estimation_and_available_models.R` | Toxicity estimation and the model set |
| `4Model_averaging_and_multimodel_inference.R` | Model averaging and multi-model inference |
| `5Response_data_and_statistical_distributions.R` | Response data and statistical distributions |
| `6Priors_and_Bayesian_inference.R` | Priors and Bayesian inference |
| `7Example_case_study.R` | A worked case study |
| `8Factor_covariates_and_groupings.R` | Factor covariates and groupings |

## Running them

Open `cr_modelling_training.Rproj` first. Each script sets the working
directory to `vignettes/`, because that is where the data files and the
saved model objects sit, and every path in the module is written relative
to it.

The model fits are loaded rather than sampled. Fetch them once with
`source("vignettes/fetch_fits.R")`, which takes about half a minute; a
script that loads a fit stops with an instruction if they are absent.

Every block runs unless the line above it says otherwise. A fit is commented
out where it would take more than 50 seconds, and so is a block that refers
to objects the script does not create. The line above each one says which,
and gives the estimate for a fit. Remove the leading `#` to run a commented
block yourself.

The fits that do run are the short ones, and most of their time is Stan
compiling the model rather than sampling it. Their `save()` is commented out
so that the fits you downloaded are not overwritten, and the `load()` after
each one restores the distributed copy so that what follows matches the page.

## Regenerating

These files are generated from the module sources and are not edited:

```bash
Rscript scripts/generate_live_scripts.R
```

Run it after any change to a module, or the script and the page will
disagree while the room is looking at both.
