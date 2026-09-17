# Specification for the intermediate revision

Companion to `notes/intermediate-revision-human.md`, which gives the decisions and
the order of work. This document holds the evidence behind each one, the exact
destination of each addition, and the alternatives that were rejected.

Written 2026-09-16 against `bayesnec` 2.1.3.35 as installed under WSL, `dev` at
2.1.3.36, `brms` 2.23.0, R 4.6.1.

Each item below has three fields. *Source* names where the material came
from. *Destination* names the module and the section it belongs in. *Code* states
whether the chunk is executed at render or shown with `eval: false`, following
the constraint recorded under *Constraints that apply throughout*.

---

# Verified state of the package

Everything in this section was checked against the installed package rather than
taken from a vignette. A vignette on `dev` may describe a function that the
installed version does not export, and the course must not name one.

## Functions present in the installed package

`check_sampling()`, `screen_models()`, `pull_best()`, `check_fit()`,
`curve_params()`, `bnec_record()`, `dispersion()`, `get_priors()`,
`failed_models()`, `bnec_group()`, `ecnsec()`, `compare_estimates()` and
`average_estimates()` are all present, confirmed by
`exists(f, where = asNamespace("bayesnec"))`.

Signatures, read from the installed package:

```
check_sampling(x, rhat_cutoff = 1.01, ess_cutoff = 400, divergence_cutoff = 10)
screen_models(x, rhat_cutoff = 1.01, ess_cutoff = 400, divergence_cutoff = 10,
              quiet = FALSE)
check_fit(x, group = NULL, ndraws = 1000, seed = 10, ...)
curve_params(object, summary = TRUE, xform = identity, ...)
dispersion(model, summary = FALSE, seed = 10)
get_priors(object, ...)
pull_best(object, ...)
bnec_record(x)
nsec(object, sig_val = 0.01, resolution = 200, x_range = NA, xform = identity,
     prob_vals = c(0.5, 0.025, 0.975), ..., dpar = NULL)
ecnsec(object, nsec, resolution = 200, x_range = NA, type = "absolute",
       xform = identity, prob_vals = c(0.5, 0.025, 0.975), ...)
```

`bnec()` formals: `formula`, `data`, `x_range`, `resolution`, `sig_val`,
`loo_controls`, `x_var`, `y_var`, `trials_var`, `model`, `random`, `random_vars`,
`prior`, `prior_type`, `timeout`, `model_survival`, `...`.

## Absent from the installed package

`models()` takes `object` only. The `max_pars` argument that `example2b` on `dev`
describes is not in the installed version. Do not write
`models("decline", max_pars = 3)` into a module until the installed package
exports it. State the underlying point in prose instead: a design with few
distinct concentrations cannot support a five-parameter equation.

The `models()` groups in the installed version are `nec`, `ecx`, `all`,
`bot_free`, `zero_bounded`, `decline` and `hormesis`. Module 4's list of these is
correct as it stands.

## The weighting default

`bayesnec:::define_loo_controls` in the installed package:

```r
if (missing(loo_controls)) {
  loo_controls <- list(fitting = list(), weights = list(method = "pseudobma"))
} else {
  loo_controls <- validate_loo_controls(loo_controls, family_str)
  if (is.null(loo_controls$weights$method)) {
    loo_controls$weights$method <- "pseudobma"
  }
}
```

The default is `pseudobma`. The `is.null()` branch exists because `pull_out()`
and `update.bnecfit()` pass the method read off the object being operated on,
which is `NULL` for an object that recorded none. `loo::loo_model_weights()`
resolves a `NULL` method through `match.arg()` to its own first choice, which is
stacking. That is `bayesnec` issue #320 reached by a second route.

Module 4 currently asserts the opposite. See *Correcting the weighting default*.

---

# Module 2, fitting a single model

## The scale of the predictor

Source: `example1`, section "The scale of the predictor".
Destination: a new section after "The formula" and before "Calling `bnec`", so
that the decision is made before the first fit.
Code: executed. The diagnostic needs no fit.

The vignette measured one algal dataset on both axes. Fitting `nec4param` on
`log(dose)` rather than on the raw dose axis reduced the maximum lack of fit in
the mean from 0.2115 to 0.0676, on a response spanning about 0.7, while five
equations fitted on the log axis all landed within 1.2 `elpd_loo` of one another.
The vignette draws the conclusion that the scale generally matters more than the
choice of equation.

The diagnostic is the spacing of the design rather than the shape of the
response. The coefficient of variation of the gaps between adjacent distinct
concentrations is small on whichever scale the series was laid out on:

```r
spacing_cv <- function(x) sd(diff(sort(x))) / mean(diff(sort(x)))
conc <- c(0.01, 0.3, 1, 3, 10, 30, 100, 300)
round(c(linear = spacing_cv(conc), log = spacing_cv(log(conc))), 2)
```

Four further points belong in this section.

A transformation is written inside `crf()`, as `crf(log(x), model = ...)`.

The transformation is evaluated before `bayesnec` inspects the predictor, so a
control recorded as `0` arrives as `log(0)` and the fit stops. `bayesnec` does
not correct a zero on the predictor axis, because a concentration of zero is an
ordinary control and no response distribution constrains the values a predictor
may take. The offset is the analyst's to choose and to add to the data first.

An offset applies to everything read off the fit, which is what `xform` does:
`ecx(fit, xform = function(x) exp(x) - 1)`.

Estimates are returned on the fitted scale. `example9` makes the point that a
logged value can land in the same numeric range as a concentration and pass
unremarked, and that reporting a logged value as a concentration understates the
threshold by a factor of e for every unit on the logged scale.

Rejected: putting this in module 5 beside the other data-preparation material.
Module 2 fits the first model, and a participant who meets the decision only in
module 5 has already fitted three models without it.

## Preparing the response

Source: `example1`, section "Preparing the response"; Ritz et al. (2026).
Destination: a new section after "The example dataset".
Code: prose and one `eval: false` chunk. No fit.

Ritz et al. (2026), Table 4, scenario 2, six control replicates, ED10:
normalisation gave a bias of 6.77% and a coefficient of variation of 26.38%, with
coverage of 0.90; the growth-rate-based approach gave 2.07%, 12.65% and 0.95.
Table 3, scenario 1, one control replicate, ED10: coverage was 0.57 under
normalisation against 0.93 unnormalised. The mechanism is Jensen's inequality
applied to the function taking y to 1/y, which biases the normalised inhibition
trend downwards and the effective doses read off it upwards.

Ritz et al. state the exemption explicitly, and the sentence should be quoted:
"If the average growth rate in the control group was known beyond uncertainty,
perhaps from historical data, there would be no problem."

Six points belong in this section.

Supply the response as measured. Do not convert it to percent inhibition, percent
of control, or percent of the observed maximum before calling `bnec()`.

Nothing is lost by not normalising, because the concentration at which inhibition
rises by x per cent is the concentration at which the response falls by x per
cent. `ecx(type = "absolute")` already returns the second, evaluated separately
within every posterior draw, so uncertainty about the control level propagates
into the interval instead of being discarded.

Dividing by the observed maximum is worse than dividing by the control mean, on
three counts. An extreme order statistic is more variable and more biased than a
mean of three to six values. The divisor then depends on every treatment rather
than on the controls alone. And it forces exactly one observation to exactly 1,
which is outside the open support of the Beta distribution.

`bnec()` detects both practices and issues a message. The maximum check is
suppressed where the response lies on a rational grid, because 19 of 20 surviving
is a genuine count-derived proportion and is not evidence of anything.

Where a divisor is unavoidable, as it is for a `Beta` response that exceeds 1, it
must be a constant fixed in advance of the analysis. ECx is invariant to the
choice of a constant divisor, because it is a relative decline from the fitted
`top`.

Weimer et al. (2012) supports the last point from the frequentist side and should
be cited beside it. A linear transformation of the data changes neither the EC50
nor the Hill slope estimate, provided the fitting function is kept and no
parameter is fixed.

## Seeds

Source: `example9` setup chunk; `example1` chunk `disp-fits`.
Destination: the existing "Reproducibility" section, which names `set.seed()`
alone.
Code: executed. The existing fit already sets a seed.

`seed` passed to `bnec()` reaches Stan's sampler. The initial values `bayesnec`
generates are drawn in R, so `set.seed()` is required as well. Neither alone
makes a fit reproducible.

## Sampler settings

Source: `example9`, step 2.
Destination: extends the existing paragraph stating that `bayesnec` raises
iterations from 2,000 to 10,000.
Code: prose only.

`bayesnec` takes `chains = 4` and `iter = 1e4`, and sets `warmup` to
`floor(iter / 5) * 4`, discarding four fifths and retaining 8000 draws. `brms`
discards half, so its own defaults of `iter = 2000` over four chains retain 4000.
That is more than a `bayesnec` call at `iter = 4000, chains = 2`, which retains
1600, from a quarter of the iterations. Compare retained draws rather than
`iter`, because the fraction discarded is a choice each package makes for itself.

Name the two arguments worth reaching for on a marginal fit. Raise `iter` and
`warmup` where R-hat is marginal or the effective sample size is low. Set
`control = list(adapt_delta = 0.99, max_treedepth = 12)` where there are
divergent transitions or the sampler reports hitting the tree depth. Neither
corrects a misspecified model.

## Sampler diagnostics

Source: `example9`, step 3; Vehtari et al. (2021).
Destination: "Checking the fit", after the trace plots.
Code: executed. `check_sampling()` runs on the existing `bnec_fit`.

The `check_sampling()` defaults are `rhat_cutoff = 1.01`, `ess_cutoff = 400` and
`divergence_cutoff = 10`. The first two follow Vehtari et al. (2021). The third
does not, and `example9` says so: Stan's own guidance is that any divergence
means the sampler failed to explore the posterior, and the value of ten is a
package convention reflecting that these non-linear models commonly produce a
small number of divergences near a parameter boundary.

The effective sample size floor of 400 is absolute, corresponding to 100 draws
per chain at four chains. Read `min_ess` rather than `min_ess_ratio` when
deciding. The ratio is read alongside it, and separates a model that cleared 400
because many draws were taken from one that cleared it efficiently.

Both the effective sample size and the divergence cutoffs are absolute counts,
and the number of retained draws changes what each demands, in opposite
directions. At 1600 draws a model must reach an efficiency ratio of 0.250 to
clear 400, against 0.050 at 8000, so that cutoff is harder to clear at reduced
settings. Ten divergent transitions is 0.63 per cent of 1600 and 0.125 per cent
of 8000, so the divergence cutoff is harder to clear at the defaults. A screen
result is a property of the settings as much as of the model, and the settings
have to be reported with it.

Define a divergent transition in plain words. The sampler's trajectory through
parameter space goes numerically unstable, which happens where the posterior has
curvature too sharp for the step size adapted during warmup. Taking more draws
samples the same surface for longer and does not smooth it.

## The record of what the fit did to the request

Source: `example9`, step 2; `bayesnec` NEWS for `bnec_record()`, issues #261 and
#93.
Destination: a short section after "The summary".
Code: executed on the existing fit.

`bnec_record()` returns `requested`, `attempted`, `excluded` with a reason per
equation, and `substitutions`. `requested` is partitioned exactly by `attempted`
and `excluded$model`. An equation that was attempted and failed to fit stays in
`attempted` and appears in `failed_models()`.

Before this existed, both were reported by `message()` and then discarded, so
neither survived a knitted document or a call wrapped in `suppressMessages()`.

`substitutions` reports the changes `check_data()` made to the response because
the family cannot represent a value as recorded. Under `Beta` a zero becomes one
tenth of the smallest non-zero value and a one is reduced by an absolute 0.001.

## Checking the fit rather than the sampler

Source: `example9`, "The assumptions the diagnostics test" and step 4.
Destination: a new section after the sampler diagnostics.
Code: executed on the existing fit, subject to the caveat below.

Convergence diagnostics ask whether the sampler explored the posterior, and `loo`
and `waic` rank equations against one another. None asks whether the variability
the fitted model implies matches the variability in the data.

The NSEC depends on that assumption most directly. `nsec()` sets its reference at
a quantile of the posterior of the control mean, and the width of that posterior
is set by the dispersion the fit estimated at the control. A model simulating
more variability at the control than the data show places the reference lower, so
the curve crosses it further along and the reported concentration is higher and
less protective. A model simulating less variability returns a lower
concentration. The direction is not known in advance and no convergence
diagnostic reveals it.

One statistic for the whole fit does not detect this. Every family with a free
dispersion parameter estimates that parameter from the same residuals such a
statistic summarises, so the parameter absorbs the discrepancy and the fit
reports a ratio near 1 while simulating the wrong spread at the concentrations
that determine the estimate. The check is therefore made within concentration
groups, which is what `check_fit()` does.

`check_fit()` returns, per group, `n`, `obs_sd`, `sim_sd`, `sd_ratio`,
`mean_ratio`, `ppp_mean`, `ppp_sd`, and a logical `control` flagging the lowest
value of the predictor. The posterior predictive p-value is `mean(sim >= obs)`,
so a value near 0.5 says the observed statistic sits in the middle of what the
model simulates and a value near 0 or 1 says it sits in a tail. Read `ppp_mean`
and `ppp_sd` rather than the ratios.

`plot()` on the returned object shows the observed statistic against the 95% span
of what the model simulates, with location and scale in separate panels because
they fail independently. The printed message flags a p-value outside [0.05, 0.95]
and the bars are the 95% span, [0.025, 0.975], so a group can fall beyond the
message's criterion and inside the bar. Read the message as a prompt to look at
the group rather than as a finding that the fit is inadequate.

`pp_check(fit, type = "dens_overlay_grouped", group = "x")` is the display for a
continuous response and `type = "bars_grouped"` for a discrete one. State why
grouping matters: pooled over the whole dataset a model can match the pooled
distribution while being wrong at every concentration.

Caveat to check before writing the chunk. `nec_data` has a continuous predictor
with 100 distinct values rather than a designed series, so grouping by
concentration may give one observation per group. Where the grouping is
degenerate, use `check_fit()` alone and show `pp_check()` with `eval: false`, or
group on a rounded predictor and say in the text that this is what was done.

## The control convention

Source: `example9`, step 4.
Destination: module 3, beside the ECx reference, cross-referenced from module 2.
Code: prose only.

The control is the lowest value of the predictor, which is `bayesnec`'s
convention throughout. Both estimators read their reference from the posterior of
the fitted mean there. `nsec()` sets its significance reference from it, and
since version 2.2 every ECx is measured from it as well.

Identifying the control by position rather than by value keeps the convention
intact under the common preparation in which a small constant is added so that a
zero control can be plotted on a log axis. An experiment whose lowest treatment
is not a control is being analysed on that assumption and should say so.

---

# Module 3, toxicity estimation and the model set

## The undefined ECNSEC

Source: the module's own hormesis section.
Destination: a definition at first use, in the same section.
Code: prose, with `ecnsec()` named.

Module 3 currently reads "`bayesnec` measures every **EC~x~**, *NSEC* and ECNSEC
from the control". ECNSEC appears nowhere else in the course and is not defined.
Parent CLAUDE.md §13 rule 4 requires every coined term to be defined at first
use.

Define it. The ECNSEC is the effect size, as a percentage, that the NSEC
concentration corresponds to on the fitted curve. It answers the question the
NSEC leaves open, which is how much effect the no-significant-effect
concentration represents. `ecnsec()` computes it, taking `type` with the
same four values as `ecx()` and defaulting to `"absolute"`.

This is needed now rather than later. Warne et al. (2025) states that an NSEC
with a measured effect size above 20% should probably not be used to derive a
guideline value, and the ECNSEC is how that effect size is obtained.

## The type vocabulary and its history

Source: `bayesnec` NEWS 2.2.0, "Breaking changes to ECx, NSEC and ECNSEC";
`example2b`.
Destination: the existing `ecx()` subsection, which describes the four values
correctly and does not describe the history.
Code: prose only.

Four values, each measured from the control and differing in what they measure
towards:

| `type` | target |
|---|---|
| `"absolute"`, the default | 0 |
| `"relative"` | the equation's theoretical asymptote, `bot` where it has one, 0 otherwise |
| `"range"` | the lowest response the curve predicts over the predictor range |
| `"direct"` | a response value supplied as `ecx_val` |

Three changes that a returning analyst needs. `"range"` is what `"relative"`
computed up to 2.1.3, and supplying `type = "relative"` explicitly now warns and
names `"range"`. `"relative"` is refused where the bound is infinite, which is an
equation with no `bot` under a family unbounded below, because there is then no
denominator. And `hormesis_def` has been removed from `ecx()`, `nsec()`,
`ecnsec()`, `compare_estimates()`, `compare_posterior()` and
`average_estimates()`; a call still passing it is refused by name rather than
absorbed by `...`.

## The censored estimate

Source: `bayesnec` NEWS 2.2.0, issues #39 and #325; Fisher and Fox (2023).
Destination: a new subsection after `nsec()`.
Code: prose only.

A target the curve never reaches within the predictor range returns `NA`, with a
warning naming how many draws were affected. Earlier versions returned the grid
point nearest the target, which for such a curve is the highest concentration in
the series, reported as an estimate with nothing said.

An ECx summarised over the remaining draws is a censored estimate and has to be
reported as one, naming how many draws were excluded and stating that the value
is bounded below by the highest concentration tested. Reported without that label
it reads as an ordinary estimate and understates the true value.

The opposite case returns the control concentration. The NSEC reference is the
`sig_val` quantile of the control posterior, so `sig_val` of the draws have a
control at or below the reference and reach it at the control itself. This is
what produces the lower credible bound of zero that Fisher and Fox (2023) report
in their Table 3 at every significance level above the 0.025 quantile the bound
is read at.

Give the rule from the paper, because a participant meeting a lower bound of zero
will otherwise read it as a failure. Where the significance value for testing
against the control is higher than the value used to calculate the lower bound of
the NSEC, typically 2.5% for a 95% credible interval, the estimated lower bound
will be zero.

## The evidence for the N(S)EC

Source: Fisher et al. (2023, *IEAM*).
Destination: expands the existing "Choosing a toxicity estimate" section, which
gestures at the paper through a screenshot and states its conclusion without a
number.
Code: prose only.

Findings from the simulation study, over four scenarios and eight designs:

- where data were generated from a threshold model the model-averaged N(S)EC was
  close to the true NEC, because the threshold equation took the high weight;
- where data were generated from a smooth model the N(S)EC was close to the NSEC,
  and the NEC estimated from those data was higher than even the EC10;
- fitting a threshold equation to smooth data implied an actual effect of 15% to
  40% depending on the design, and the implied effect did not fall as sampling
  effort increased;
- the estimated effect at the N(S)EC and at the NSEC approached 10% on average
  for the worst-designed experiment simulated, and about 1% for the most
  replicated;
- ECx estimates were themselves biased upwards at low sampling effort;
- even for the poorest design simulated, the N(S)EC was lower than the EC10.

Two case study results are usable directly. For *Cryothecomonas armigera* under
copper the weights spread across three smooth equations at 0.334, 0.377 and
0.154, with the threshold equations at 0.001 and 0; the N(S)EC was 7 (1 to 11)
µg/L against an NEC of 20.9 (14.9 to 23.5), and the estimated effect at the
N(S)EC was 2.3%. For *Stomopneustes variolaris* the NEC was definitively higher
than the EC10, the N(S)EC was 140 µg/L, and a NOEC computed on the same data was
31.8 µg/L.

That last comparison is the one worth teaching. It answers the objection that the
N(S)EC is overly conservative, by showing that a NOEC on the same data was more
than four times lower.

## Experimental design and the estimate

Source: Fisher et al. (2023), Figure 2B and Supporting Table S1; Krull (2020).
Destination: a new section, placed after the NSEC and the NOEC are both defined,
since the argument is about which design each needs.
Code: prose only.

For the NEC 1 scenario the weight on the generating model was 0.987 for twelve
concentrations with five replicates of ten trials, which is 600 trials, and 0.984
for eight concentrations with five replicates of twenty trials, which is 800.
The paper's conclusion is that replication within treatments, which is a
requirement of the one-way analysis of variance used to generate NOECs, should be
reduced in favour of increasing the number of treatments, and that redistributing
replication this way can often be achieved without spending more on the
experiment.

Krull (2020) gives the limitation that belongs beside it. Simulating threshold
estimation by maximum likelihood, by MCMC and by piecewise regression across
slopes, background mortalities and designs, all methods performed poorly on
shallow and intermediate curves, and accuracy increased with the slope of the
curve. Information criteria weights usually did not provide strong evidence for
the true model.

## The guideline position

Source: Warne et al. (2025), sections 3.2.4 and 3.4.2.1.
Destination: the existing "NSEC against ECx" section, which cites the 2008 and
2018 guidelines only.
Code: prose only.

The hierarchy is now: NEC and NSEC as the most preferred, directly estimated
negligible-effect concentrations; EC/IC/LCx where x is at most 10, and BEC10, as
other appropriate estimates; EC/IC/LCx above 10 and up to 20, and the NOEC, as
less preferred; and LOEC, MATC and EC/IC/LC50 as requiring a conversion factor.

The condition to teach is stated in section 3.4.2.1: "any NOECs and NSECs with an
effect size > 20% should probably not be used to derive GVs".

## ECx from threshold equations

Source: `example2b`; `example4`.
Destination: the discussion of what each equation group supports.
Code: prose only.

ECx estimates can be obtained from both `nec` and `ecx` equations, and will
usually be lower, meaning more conservative, for `ecx` equations fitted to the
same data. The package recommends the `all` set where ECx estimation is required,
because `nec` equations fit some datasets better and the averaging places the
greatest weight on whichever suits the data.

`example4` demonstrates the opposite direction on data simulated from a smooth
curve. There the threshold equations produced higher EC10 values, because a
broken-stick curve fitted to a smooth decline is flat before the break and then
falls sharply.

---

# Module 4, model averaging and multi-model inference

## Correcting the weighting default

Source: the installed `bayesnec:::define_loo_controls`; `example2`; issue #320.
Destination: the existing "Weights in `bayesnec`" section.
Code: executed. The chunk already runs; its printed output changes.

The module currently states: "Stacking is the default. `bayesnec` passes no
method of its own, so `loo`'s own default applies, and the summary above names
it: `Method: stacking_weights`."

`pseudobma` is the default and `bayesnec` passes it explicitly. Keep the history
in one sentence, because it is the reason a published analysis may report either
method: earlier versions left the method unset in some paths, and
`loo::loo_model_weights()` resolves an unset method to stacking, so a set
assembled by `c()` or `amend()` was weighted by stacking while the same set
fitted by `bnec()` was weighted by pseudo-BMA.

Keep the existing description of what the two methods do, which is accurate. Add
the package's position from `example2`: stacking is not recommended at the sample
sizes concentration-response experiments produce, because the motivation for
averaging here is to capture model uncertainty rather than to reduce prediction
error.

The summary output must be re-checked. The module quotes `Method:
stacking_weights` from a rendered summary, and the chunk that produced it will
now print `pseudobma`. The quoted text in the prose has to be changed to match
and the `_freeze` entry for module 4 cleared.

## Screening a model set

Source: `example2`; `example9`, steps 3 and 5.
Destination: replaces and extends "Diagnostics across a set".
Code: executed. Both functions run on the existing three-model `bmanecfit_more`.

`rhat()` on a `bayesmanecfit` returns a list keyed by model name, each element
holding the R-hat values and a logical `failed`. `rhat(fit)$failed` is `NULL` on
that structure and silently drops nothing, which module 7 already documents and
module 4 does not. The default cutoff is 1.01, following Vehtari et al. (2021);
the older 1.05 is available through `rhat_cutoff` and is useful mainly for
comparison against analyses predating the current recommendation.

`check_sampling()` reports all three diagnostics per candidate. `screen_models()`
screens on the three together, drops what failed, and reports what went and why,
one line per model, naming every threshold it missed rather than only the first.
That message is the record of the exclusion, and `quiet = TRUE` should not be
used in an analysis that will be reported.

`amend()` does the work, and is the only route by which a model leaves a
`bayesmanecfit`. The weights of a screened set are therefore read from the
screened fit rather than rescaled by hand, because `amend()` recomputes the
weights across what remains.

Where every candidate fails, `screen_models()` stops with an error listing the
reasons rather than returning a meaningless average. Where none fails, it says so
and returns the set unchanged.

The module also needs a correction. Its current text reads "A poor model with a
negligible weight does no harm, because it cannot influence the averaged
prediction." That holds for shape and does not hold for convergence. A model that
has not sampled contributes to the averaged prediction on the strength of a
posterior the sampler never explored, and its weight is computed from that same
posterior, so the weight is no evidence that keeping it is safe.

## Kinds of sampler failure

Source: `example9`, "Kinds of failure" and "Summary".
Destination: its own subsection within the screening material.
Code: prose, with the `example9` outcomes quoted and attributed.

This is the most valuable single addition in the revision.

An equation that misses the effective sample size or R-hat threshold with no
divergent transitions was sampled from a posterior the sampler could traverse. It
was not sampled for long enough. Refitting with a higher `iter` is what it needs,
and excluding it on this evidence is premature.

An equation that fails on divergent transitions is describing the posterior
rather than the run. More draws sample the same surface for longer. Raising
`adapt_delta` reduces the step size and does not remove a ridge in the
likelihood.

An equation that fails to fit at all, and is the most heavily parameterised in
the set on a design with few distinct concentrations, has a posterior that is
flat in some direction. Sampling a flat direction for longer does not locate a
maximum that is not there.

The `example9` outcomes are the evidence. On `simazine` at deliberately reduced
settings the screen removed five of eleven equations. At the defaults on the same
data it removed one, with a further equation failing to fit at all. On
`nassarius` at those same defaults it removed most of the set, on divergent
transitions in the hundreds and thousands. The `nassarius` result is a statement
about the experiment: across the untested interval from 1.25 to 2.5 the threshold
is unidentified and the likelihood has a flat ridge there.

`failed_models()` retains the priors and initial values each failed equation was
given, which are built inside `bnec()` and are otherwise unrecoverable, and they
are the starting point for an adjustment. State the warning that goes with it: a
prior wide enough to let an over-parameterised equation fit would be doing the
work the data cannot, and the estimate it produced would be a property of that
prior.

## Parallel fitting

Source: `example2`, "Fitting the model set in parallel"; `bayesnec` issue #184.
Required by repository CLAUDE.md §8, which records that module 2 forward-
references module 4 for this and that module 4 does not yet deliver it.
Destination: a new section after the model-set material.
Code: `eval: false` throughout. Running any of it would refit the set.

`bnec()` and `amend()` fit their models under whatever `future` plan is set.
There is no argument, because the plan already records how much of the machine
you are willing to use. `future` and `future.apply` are Suggests, so a session
without them takes the sequential path.

The two levels nest. Under a plan of more than one worker `bnec()` passes
`cores = 1` to `brm()`, because `workers x chains` would otherwise be requested,
which is sixteen processes for four workers and four chains. Passing `cores`
yourself is left alone and is how the two levels are nested deliberately.

`plan(multisession)` is the portable choice. `plan(multicore)` forks, which is
unavailable on Windows and inside RStudio and Positron, where the plan resolves
to a single worker. `bnec()` reports the number of workers it was given, so a
collapsed plan is visible rather than merely slow.

Use the `cmdstanr` backend under a forked plan. `rstan` compiles into the
session's temporary directory, which a forked worker shares with its parent.

Check a plan before committing a long run to it. Fit two models with a small
`iter` and confirm it returns.

The timings, attributed as `example2` attributes them. The whole `decline` set,
measured 2026-09-12 on an otherwise idle 22-core workstation under WSL2, R 4.6.1,
`cmdstanr` and `plan(multicore)`, each run from a fresh session with the compile
cache emptied first: no plan 299.5 s, two workers 226.4 s, four workers 174.2 s,
eight workers 149.9 s. Against a warm cache the same five configurations took 99,
114, 108, 101 and 94 seconds, so two workers then took longer than no plan at
all.

The point to teach is which run is being timed. A first fit on a machine where
the Stan programs are not built is compilation-bound and parallelises well. A
repeat against a warm cache is sampling-bound, and for models this small nothing
is gained. Where a single model takes minutes rather than seconds the sampling
dominates in both cases.

State the reproducibility caveat. Each fitted model reproduces exactly under a
plan for the same `seed`. The model-averaged quantities do not, because `bnec()`
draws the seed for the weighted posterior from the session's random number stream
after the models are fitted, and a sequential run has advanced that stream while
a parallel run has not. `set.seed()` in the calling session before `bnec()` fixes
it.

## Reading the weights

Source: `example9`, step 2; the `pull_best()` entry in `bayesnec` NEWS.
Destination: the existing summary and `pull_out()` material.
Code: executed on the existing set.

`wi` is the share of the model-averaged prediction each candidate contributes,
summing to one across the set. Where a model reports Pareto k values above 0.7
for a few observations, `brms` warns; the weights are usable and the estimates
behind them are flagged as unreliable at those points, which belongs in a methods
statement alongside the weights themselves.

`pull_best()` returns the highest-weighted candidate as a `bayesnecfit`, and
returns a `bayesnecfit` unchanged. It reports the weight with the number of
candidates it was selected from, because the two together say whether the
selection means anything: the highest weight of a flat set of ten is little more
than the equal share of 0.1.

This corrects the module's existing treatment, which tells the reader to read the
weights out of `mod_stats` and pass a name to `pull_out()`.

## Curve parameters

Source: `bayesnec` NEWS, issue #297; `example9`, step 5.
Destination: "Elements of the fitted object", beside `mod_stats`.
Code: executed on the existing set.

`summary()` reports the model weights, the per-equation dispersion, the weighted
no-effect estimate and the per-equation Bayesian R-squared, and no parameter
estimates. For a model average nothing returned them before `curve_params()`.

The estimates are per equation and are not averaged across the set, because the
equations do not share a parameter list. `ecxexp` has no `bot`, the
three-parameter equations have no `d`, and only the `nec` group estimates `nec`.
Averaging a parameter over whichever equations estimate it would average over a
different subset for each parameter. The model weight is reported beside each row
instead.

`xform` reaches `nec` and `ec50` and no other parameter, because those two are
measured on the predictor axis while `top`, `bot` and `beta` are not.

---

# Module 5, response data and statistical distributions

## The dispersion screen

Source: `example1`; `example9`, "The dispersion screen"; `bayesnec` NEWS, issue
#262.
Destination: the existing "Binomial" and "Beta-binomial" sections, replacing the
judgement by eye.
Code: executed. `dispersion()` runs on the existing `exp_1nec` fit.

The module currently judges over-dispersion by eye, in the sentence "The bounds
above are implausibly narrow for data of this kind", and never names
`dispersion()`.

`dispersion()` is defined for `poisson` and `binomial` only, the two families
whose variance is fixed by the mean, and returns an empty vector otherwise. For
`Beta`, `gaussian`, `Gamma`, `negbinomial` and `beta_binomial` there is no
over-dispersion decision to take, because the family already has a free parameter
for the spread.

`dispersion(fit, summary = TRUE)` returns, for each posterior draw, the ratio of
the observed to the simulated sum of squared Pearson residuals, summarised to a
median, an equal-tailed 95% interval, and `P(>1)`, the posterior probability that
the ratio exceeds 1. The statistic centres near 1 when the family's variance
assumption holds. Read `P(>1)` rather than the point estimate, which says nothing
about how well determined the ratio is. It is symmetric, so `1 - P(>1)` answers
the under-dispersion question.

Each draw compares one observed residual sum against a single simulated
replicate, so the spread of the posterior is dominated by replicate-to-replicate
simulation noise rather than by uncertainty about a dispersion parameter. At a
small design its power to separate a genuine departure from noise is low, and a
rule that asks only whether the interval excludes 1 will detect gross
over-dispersion and nothing else.

The screen is model-conditional. A dispersion ratio cannot separate
over-dispersion from a shape that fits badly, so the screening fit has to sample
adequately before the dispersion it reports is used.

`beta_binomial` mixes the success probability over a beta distribution, which can
only add variance. It addresses over-dispersion and does not address the opposite
case. A screen pointing below 1 indicates something about the fit, the design or
the data-generating process rather than about the family.

`summary()` on a model set reports `dispersion_Estimate` and
`dispersion_P_over_1` per equation.

## Counts observed over an exposure

Source: `example1`, "Counts observed over an exposure".
Destination: the existing "Poisson" section.
Code: `eval: false`.

A count is often observed over an exposure that differs between replicates:
offspring per replicate where the number of fecund adults differs, cells over
different areas, events over different lengths of time. The `rate` term declares
that exposure and is available for `poisson` and `negbinomial`:

```r
y | rate(n_females) ~ crf(x, model = "nec4param")
```

Raw counts are passed unchanged, so a replicate with twenty females contributes
more information than one with eight, which is what fitting `y / n_females` as a
continuous response would discard.

Because `bnec()` uses an identity link, `brms` applies the denominator
multiplicatively on the response scale rather than as a log offset on the linear
predictor. The fitted mean therefore is the rate, and `top`, `bot` and `nec` are
directly interpretable as offspring per female. The prediction grid holds the
denominator at 1, and `autoplot()` divides the observed counts through to match.

The caveat is the teaching point. Where the exposure varies for reasons unrelated
to the treatment, `rate` does what is wanted. Where the exposure is itself an
effect of the treatment, as it would be were the denominator the number of
females surviving to the point of counting, the endpoint has changed from
reproductive output to fecundity conditional on survival. Both are legitimate and
they are different endpoints.

## Censoring

Source: `example1`, "Censoring" and "Limits of the censored likelihood"; Helsel
(2006).
Destination: a new section after the family-by-family examples.
Code: `eval: false`. The `alga` results are quoted as a table with their source
named.

A value is censored when the truth is known to lie in an interval rather than at
a point: a concentration below a limit of detection, a count below the resolution
of the method, a value rounded at the recording precision. The number written
down is a bound.

Helsel (2006) is the citation for why substitution is not acceptable:
"Substituted values using a fraction anywhere between 0 and 0.99 times the
detection limit are equivalently arbitrary, equivalently precise, equivalently
wrong." He is equally clear that deleting non-detects is worse, producing a
strong upward bias in every subsequent measure of location.

`bayesnec` passes through the `cens()` aterm from `brms`:

```r
y | cens(censoring) ~ crf(x, model = "nec4param")
```

`censoring` is a variable in the data rather than a constant, because only some
rows are censored. It takes `"none"`, `"left"`, `"right"` or `"interval"`,
equivalently `0`, `-1`, `1` and `2`. Interval censoring takes a second argument
giving the upper bound.

The response column holds the bound. A left-censored row asserts only that the
truth lies between the family's lower support and that bound.

A censored row is exempt from the boundary shifts `bnec()` otherwise applies,
because altering the bound would restate what is known. Declaring a row censored at a value
the family excludes is an error rather than a shift.

The mechanism. A censored row contributes F(bound), the probability of falling at
or below the bound, rather than a density at a point. As the fitted curve
descends past the bound, F rises toward 1 and saturates: the curve is not
rewarded for passing further below and is not penalised. A substituted value
contributes a density, which falls away as the curve descends, so substitution
actively pulls the curve back up toward the value that was substituted.

The limit has to be taught beside the mechanism. Saturation is what confines the
censored likelihood to what is known, and it is also its limit. Once the curve sits a few
residual standard deviations below the bound the likelihood is flat there, so a
parameter whose only expression is in that region is not identified, and what
gets reported for it is the prior. The lower asymptote is the usual casualty.

The worked evidence, from the `alga` dataset in `example1`. Fitting `nec4param`
with a Gaussian family and the same prior throughout, varying only how the
below-limit rows are declared:

| below-limit rows declared as | `bot`, median [95% CI] | posterior SD |
|---|---|---|
| the bound, treated as observed | −1.985 [−2.013, −1.958] | 0.014 |
| left-censored at the bound | −5.64 [−11.0, −2.51] | 2.23 |
| interval-censored from extinction to the bound | −2.45 [−2.68, −2.10] | 0.17 |

Substitution is wrong about the value and states it precisely: an interval a
tenth the width of any other. Left-censoring withdraws the false assertion and
reveals that the data contain almost nothing about `bot`; the chains mix, R-hat
is 1.00, and five divergent transitions in 8000 draws is nothing a routine check
would stop on. Interval censoring bounds the flat region from below and restores
a usable posterior, which then rests on the declared extinction floor.

Two points follow for practice. Read identification from the prior-to-posterior
contraction rather than from the interval: across those three fits it is 0.996,
0.379 and 0.953, and the middle value is the diagnostic the interval hides. And
the problem propagates, so an `ecx()` value whose target falls near the asymptote
deserves the same scepticism as `bot` itself.

## Zeros

Source: `example6`, "Three kinds of zero"; Martin et al. (2005); Warton (2005);
Blasco-Moreno et al. (2019).
Destination: expands the single sentence the module currently gives.
Code: prose only.

Three paragraphs rather than a full treatment.

A zero that is a measurement below the resolution of the method is censored, and
belongs in the censoring section above.

A zero that marks a distinct event, such as an organism that died so that growth
was undefined rather than merely small, is structural and belongs in a hurdle
model. A hurdle fits two blocks: one for whether the event occurred, and one for
the response given that it did not.

A zero that is an ordinary draw from a distribution admitting zeros is not
evidence of anything. Warton (2005) is the citation: many zeros does not mean
zero inflation. Test the simpler family before reaching for a two-block one.

State that `bnec_hurdle()`, `bnec_joint()` and the hurdle families exist and
point at the `bayesnec` article, rather than teaching them. A full-day workshop
cannot cover a two-block fit.

## Boundary substitution

Source: `example9`, step 1; `bnec_record()`.
Destination: the existing "Beta" section, cross-referenced from module 2.
Code: prose only.

A `Beta` response is supported on the open interval from 0 to 1, so a value
recorded as exactly 0 or 1 lies outside the family's support. `bayesnec`
substitutes for both. A zero becomes one tenth of the smallest non-zero value,
which is relative to the scale of the response, and a one is reduced by an
absolute 0.001. A series with several concentrations at a boundary is being
modelled through that substitution rather than as recorded.

`bnec_record(fit)$substitutions` names each altered value, so a methods section
can state what was modelled.

## Non-constant dispersion

Source: `example1`, "Non-constant dispersion".
Destination: a new section, kept short.
Code: `eval: false`.

By default a fit holds the dispersion parameter constant across the whole curve.
Where that is wrong the consequences fall on the credible intervals of the
toxicity estimates, and they fall unevenly, because the assumed variability is
too large at one end of the series and too small at the other.

```r
y ~ crf(x, model = "nec4param") + disp(~x)         # on the predictor
y ~ crf(x, model = "nec4param") + disp("power")    # on the fitted mean
```

Diagnose before modelling. A variance function will absorb any cause of apparent
heteroscedasticity and report a confident slope regardless. Three explanations to
exclude first: substituted or censored values clustered at one end of the series;
a misspecified family, since a constant coefficient of variation is what a
`Gamma` already implies; and lack of fit in the mean, since a variance function
is free to describe a region as noisy rather than the curve as wrong.

What it changes, and the direction, is the teaching point. On the worked example
in `example1` the `ecx` intervals narrowed while the `nsec` interval widened.
Under a variance function the tightly grouped high-concentration replicates take
more weight and the variable control replicates less, so `bot` and `beta` are
estimated more precisely and `top` less so. ECx depends on the curve as a whole
and inherits the better-determined shape. The NSEC is referenced against the
control and inherits the one thing the variance function made less certain.

An NSEC reported from a constant-dispersion fit rests on an assumption about
control variability that the fit has not checked, and it is the estimate most
exposed if that assumption is wrong.

## Growth rates and the ErCx notation

Source: OECD Test Guideline 201 (2026); `example7`.
Destination: the existing "Gaussian" section, in two sentences with a pointer.
Code: prose only.

OECD TG 201 distinguishes two response variables and two corresponding estimates.
Average specific growth rate is the logarithmic increase in biomass over the
exposure period, and the estimate derived from it is the ErCx. Yield is the
biomass at the end of the period minus the biomass at the start, and its estimate
is the EyCx. The guideline states that the ErCx will generally be higher than the
EyCx, that growth rate is the scientifically preferred basis, and that yield is
retained to satisfy current regulatory requirements.

A growth rate can be negative, which rules out families bounded below at zero and
rules out the equations whose lower asymptote is fixed at zero. `example7` covers
this in full and is the pointer.

## Model suitability rules

Source: `example2b`, "Model suitability for response types"; `example9`.
Destination: the existing "Model suitability by response type" section, which
says `bayesnec` discards unsuitable models without saying which or why.
Code: executed. `bnec_record()` runs on an existing fit.

State the three rules.

Equations with an exponential decay and no `bot` are zero-bounded and are
unsuitable for a Gaussian response or any log or logit link, because they cannot
predict negative values.

Equations with a linear decay, whose names contain `lin`, are unsuitable for
zero-bounded families under an identity link, because their fitted mean is
unbounded below.

Equations that raise the predictor to a fractional power cannot be evaluated
where the predictor is negative, which a logged concentration series below 1
routinely is.

`bnec_record(fit)$excluded` gives the reason per equation, so none of this has to
be inferred.

---

# Module 6, priors and Bayesian inference

## Priors without a fit

Source: `example3`.
Destination: "Inspecting the priors" and "Specifying priors".
Code: executed. `get_priors()` needs no fit.

`get_priors()` builds the default priors for a formula and dataset without
fitting anything, so the set can be inspected and edited before the first run.
This is the answer to the module's existing observation that `pull_prior()`
reports the priors of a model that has already been fitted.

## The predictor-scaled prior in full

Source: `example3`, "Priors for predictor-scaled parameters"; `bayesnec` issue
#302.
Destination: the existing "Predictor-scaled parameters" section.
Code: prose only.

The module describes the lognormal correctly and states the truncation. It does
not give the rule that sets the width, and it does not give the consequence for
how a control is recorded.

The width is stated rather than chosen. The central 95% interval of the
untruncated prior covers every concentration tested: the standard deviation is
the larger of the two half-widths from the mean to the ends of the logged series,
divided by `qnorm(0.975)`. A prior on a threshold should not exclude a
concentration the experiment applied, and that criterion is the whole of the
rule.

Because the width is set by the two extreme concentrations rather than by the
spread of the series between them, it responds to how a control is recorded. On
the `nassarius` contaminant A series the standard deviation is 2.30 with the
control recorded as `0`, 2.59 with it recorded as `0.001` and 6.11 with it
recorded as `1e-6`. Record a control as `0`.

The location is the median of the distinct positive predictor values, logged.
Distinct values are used rather than the observations so that the prior describes
the series of concentrations tested rather than how many replicates each
received.

Keep the module's existing statement that `nec` and `ec50` are truncated to the
observed range of the predictor, and that an estimate sitting at the edge of the
tested range should be treated with suspicion rather than reported. Give it the
numbers above.

## The regularizing set

Source: `example3`, "Selecting the default prior set"; `bayesnec` NEWS 2.2.0,
issue #305.
Destination: the existing "Two prior sets" subsection, which gives two sentences.
Code: prose only.

The `"regularizing"` set is stated once, as a location and a spread, and each
family's entry is derived from it.

The location is the mean response at the end of the concentration series where
the parameter is the level of the curve: the lowest concentrations for `top`,
which for a design with an unexposed control is the control group, and the
highest for `bot`.

The spread is 0.4 of the standard deviation of the `"uninformative"` prior for
the same parameter, floored at the standard error of the location and capped at
the `"uninformative"` width. The floor exists because a prior narrower than the
noise in its own anchor states a precision the data do not supply.

The measurement, attributed. Over 420 `top` and `bot` entries built from one
simulated response each, across five designs, three predictor transforms, seven
families and two links, the ratio is 0.4 exactly in 385, the standard-error floor
binds in 35, and the prior is wider than the `"uninformative"` one in none.

The assumption that can fail is the teaching point. The rule assumes the lowest
concentration tested is at the level `top` describes and the highest is at the
level `bot` describes. Where the highest concentration has not reached the lower
asymptote, the location for `bot` sits above the true value, and no choice of
spread corrects a bias. Use the `"uninformative"` set on such a design.

`prior_type` is ignored whenever the user supplies their own `prior`.

## Fixing a parameter

Source: `example3`, "Fixing a parameter"; Weimer et al. (2012).
Destination: a new section after "Specifying priors".
Code: the `get_priors()` call is executed; the fits are `eval: false`.

A parameter is held at a known value by giving it a `constant()` prior.
`bayesnec` has no `fixed` argument of the kind `drc` provides, because a
point-mass prior expresses the same thing through machinery that already exists.

```r
fixed_prior <- get_priors(y ~ crf(x, model = "nec3param"), data = avoid_data,
                          family = Beta(link = "identity"))
fixed_prior$prior[fixed_prior$nlpar == "top"] <- "constant(0.5)"
```

Stan does not declare a parameter whose prior is constant, so it contributes
nothing to the posterior and is reported with no error.

Two things are given up, and both are easy to miss because the fit looks better
rather than worse.

The first is a diagnostic. In a two-choice avoidance assay the expected
proportion in the treated container at zero dose is 0.5 by design. Fitting with
`top` free gives a credible interval that either covers 0.5 or does not, and that
is positive evidence about whether the assay was unbiased before the contaminant
was involved. No other output of the model supplies it. Where a real baseline
preference is present, fixing `top` at 0.5 produces a fit that looks exactly as
sound; the disagreement is absorbed into `beta` and `nec`
instead, and so into the toxicity estimate.

The second is uncertainty. A constant asserts that the value is known exactly, so
none of it propagates into ECx or the NEC, and the credible intervals narrow.
That reads as a better analysis and is the opposite.

A tight informative prior is the better tool in most cases. `beta(6, 6)` has a
mean of 0.5 and a standard deviation of 0.14, concentrating `top` around the
design expectation without asserting it.

Weimer et al. (2012) reaches the same conclusion from the frequentist side and
should be cited. Fixing parameters in a log-logistic model may produce erroneous
estimates, and parameters should be fixed only with compelling reason. This
matters most where data are background-corrected and normalised and a reduced
model fixes the response range between 100% and 0%.

One case must be refused outright. A response divided by its control mean has a
maximum of 1 by construction, and fixing `top` at 1 is a common reflex.
Normalising to an estimated control already discards the sampling variability of
the divisor, and fixing the asymptote afterwards makes that loss invisible rather
than merely present. Cross-reference module 2.

## Prior sensitivity

Source: `example3`; Depaoli et al. (2020) and Gelman et al. (2017), both already
cited in the module.
Destination: extends the existing paragraph saying priors "still need
interrogating".
Code: `eval: false` for the refit; prose for the contraction.

A sensitivity analysis refits under a different prior set and compares the
posteriors. `prior_type = "regularizing"` gives a second set without writing any
priors by hand, so the comparison takes one argument. The module already compares
default against user priors by histogram; extend that comparison to name what a
difference would mean and what to do about it.

The numeric diagnostic comes from `example1`. Prior-to-posterior contraction, one
minus the ratio of the posterior variance to the prior variance, separates a
parameter the data determined from one the prior determined. It is computed from
`pull_prior()` and the posterior draws. A value near 1 means the data determined
the parameter, and a low value means the reported posterior is the prior. The
interval alone does not show this.

## Group-level priors

Decision: do not add to module 6. This belongs in module 8, where grouping is
taught, and module 6 is already the longest module in the set. Add one sentence
pointing forward.

---

# Module 7, the worked case study

## Replacing the bespoke screening helper

Source: `example9`, step 3.
Destination: the existing "Dropping models that did not converge" section.
Code: `eval: false`, as the module's fitting code already is.

The module defines its own helper over `rhat()` and drops models with `amend()`.
Use `screen_models()` instead. Keep the existing explanation of why
`rhat(fit)$failed` is `NULL` on a `bayesmanecfit`, because it is correct and
useful, but state it as the reason `screen_models()` exists rather than as a
reason to write the helper.

Check whether `scripts/generate_herbicide_fits.R` needs the same change. Where
the stored results were produced by the R-hat helper, either regenerate them or
state in the provenance block which screen produced them. Do not change the prose
to describe a screen that did not produce the numbers on the page.

## A reporting section

Source: `example9`, step 5.
Destination: a new closing section after "Probabilities of difference".
Code: prose, with the functions named.

List what a methods section has to state:

- the data as analysed, including the response, its bounds, and any values
  altered to bring them inside the family's support;
- the family and how it was chosen;
- the candidate set as fitted, meaning the set requested, the equations excluded
  as invalid with the reason, and any model that failed to fit at all, with a
  note of whether a retry was attempted, all of which `bnec_record()` returns;
- the sampler settings, including the number of retained draws, since both the
  effective sample size and divergence cutoffs are absolute counts read against
  them;
- the diagnostic thresholds used, all three, since they are choices and they
  determine what is excluded;
- which models were excluded and why, naming the threshold each failed, and
  distinguishing a model that would pass at more draws from one the design does
  not support;
- the estimate settings, meaning `sig_val`, `type`, `resolution` and any `xform`;
- the estimates with their credible intervals, and the weights of the retained
  models;
- the parameters of the fitted curve, from `curve_params()`.

This is the natural close of the course's Bayesian thread and is the section an
intermediate practitioner will return to.

## The guideline position

Source: Warne et al. (2025).
Destination: replaces the sentence stating that use of the NSEC is "recent and
not yet formally adopted".
Code: prose only.

Warne et al. (2025) lists the NSEC alongside the NEC as one of the two most
preferred statistical estimates. Add the effect-size condition and point at
`ecnsec()` as the way to check it.

## Reporting the scale

Source: `example9`, step 5.
Destination: the existing table of estimates.
Code: prose only.

The module fits on `log(concentration)` and back-transforms for the table, which
is correct. Add the warning. A logged value can land in the same numeric range as
a concentration and pass unremarked, and only a negative logarithm announces the
error. A reported value has to state the scale it is on.

---

# Module 8, factor covariates and groupings

## Status

`vignettes/8Factor_covariates_and_groupings.Rmd` is a 29 KB `learnr` tutorial,
last touched in 2023 and not converted. Repository CLAUDE.md prohibits authoring
a new `.Rmd`, so revision means conversion.

It depends on `vignettes/functions.R`, which holds `pred_out()`, `nec.brmsfit()`
and `nsec.brmsfit()`, written for module 8 before the equivalent functionality
existed in the package. `nsec.brmsfit()` calls three unexported internals:
`bayesnec:::do_wrapper`, `bayesnec:::modify_posterior` and
`bayesnec:::min_abs`. Check every one against the installed package before
running the module. Triple-colon calls break without deprecation.

## The compendium behind the fits

The `example8` fits are not re-run here, and they are not re-run by the vignette
either. `C:/Rworking/grouping-structures` is the research compendium that
produces them.

What it records, from its own `README.md`: `example8` fits 189 models, being five
model-averaged sets over the whole candidate list, two two-equation contrasts,
and two `bnec_group()` calls fitting a complete set at each of seven and two
levels. Measured on the AIMS high-performance cluster on 2026-09-16 with 40 tasks
resident, those 189 units hold 22.2 hours of sampling between them; run as an
array they completed in 52 minutes, and the resulting store is 484 MB. The
vignette then renders against the store, reading each fit rather than making it.

Three consequences for this course.

The course does not fit these models. A 484 MB store cannot be committed to a
public repository and cannot be fetched during a workshop, so module 8 reads
compact results in the pattern module 7 established. A script derives those
results from the compendium's store, and only the results are committed.

The store is not in the working tree. It is fetched with
`grouping-structures/hpc/fetch-store.sh` and located by the `BAYESNEC_FIT_STORE`
environment variable. The generating script for module 8 therefore has a
prerequisite that must be stated in its header.

The 38.5 minute render is the warning that matters for the course. That time is
`ecx()` and `nsec()` extraction across the 14- and 18-equation sets, with nothing
sampled. A module that reads a saved `bayesmanecfit` and then calls `ecx()` on it
has not avoided the run time; the estimates themselves have to be saved.

`negative-response-conventions` is the equivalent compendium behind `example7`,
and the same reasoning applies to anything drawn from that vignette.

## Rewriting around `bnec_group()`

Source: `example8` on the open pull request `open-AIMS/bayesnec#228`.
Destination: the fitting sections of the converted module.
Code: `eval: false` for every fit, with results read from `vignettes/data/`.

`bnec_group()` fits the model set separately within each level of a factor, which
is what module 8's helpers were written to approximate. It is exported in the
installed version.

Because that branch is not merged, anything taken from it must be checked against
the installed package before it is written into a module, and the module must not
name a function the installed version lacks.

## Group-level term syntax

Source: `example8`; `example3`, "Priors for group-level parameters".
Destination: a new section before the fitting material.
Code: prose, with the three forms shown.

Three forms, and what each does. `ogl(group)` displaces the whole curve, adding
an offset parameter as well as a standard deviation. `(par | group)` places a
deviation on one named parameter. `pgl(group)` places a deviation on every
parameter of the equation.

The mechanism has to be explained, because it is why a grouped `bayesnec` fit
behaves differently from a grouped `brms` model. `brms` declares a group-level
deviation unconstrained. Where the response distribution restricts the range of
the mean, which is every family except Gaussian under an identity link, a
proposal that puts the mean outside that range is rejected by Stan and counted
as a divergent transition. `bayesnec` therefore applies the deviation
multiplicatively where it applies to the whole curve or to `top` or `bot`. On a 0
to 1 response it scales the odds of the mean, and on a positive response it
scales the mean itself. A deviation of zero leaves the curve exactly as it was,
so the parameters keep their meanings, while the mean can no longer leave its
support.

The evidence. On a simulated `beta_binomial` design across three seeds the
multiplicative form gave no divergent transitions at all, against 70 to 95 per
cent for the additive form on the same data, and it sampled about twice as fast.
On the `herbicide` data with `Beta(link = "identity")` and `nec4param`, a
`(bot | herbicide)` term gave 51 divergent transitions of 2000 at
`adapt_delta = 0.95` under the additive form, and gives none at Stan's default of
0.8 under the multiplicative one.

A deviation on `nec`, `ec50`, `beta`, `slope`, `d` or `f` still adds to that
parameter directly, because none of those is bounded by the likelihood.

## Comparing groups

The module's existing subject is whether ECx values differ between groups. Keep
it. `compare_posterior()` and `compare_estimates()` are the current route and
module 7 already introduces the first, so cross-reference rather than repeating.

## Scope limit

Module 8 is the last module of a full day. Convert it, rewrite the fitting around
`bnec_group()`, explain the three term forms and the multiplicative deviation,
and keep the comparison material. Do not import the three case studies from
`example8`; they are longer than the whole current module.

---

# Constraints that apply throughout

## Limits on chunk run time

Decided with the user on 2026-09-16. New code either reuses the fit the module
already makes, or is written with `eval: false`. The reason is `_freeze`: the
site stores computed output, and a new fit would have to be sampled on the
presenter's machine before the site could be rebuilt.

The user added a second limit on the same day, and it is the binding one. Any
chunk taking longer than a few minutes is pre-computed and saved, whether it
fits a model or not. A workshop in which a participant waits ten minutes for a
chunk is not teaching.

Post-processing crosses that limit more easily than it appears. `ecx()` and
`nsec()` rebuild the component predictions of every equation in a set on each
call, so extraction across a large model-averaged set is slow independently of
sampling. The `bayesnec` grouping vignette renders in 38.5 minutes against a
store in which all 189 fits already exist, and none of that time is sampling;
`bayesnec` issue #306 measured the same run time on `example7`. Two mitigations
apply where a module must run such a call: reduce `resolution`, whose default is
200 and whose precision saturates well below that, and read the estimate once
into an object rather than calling `ecx()` repeatedly.

Functions that run on an existing fit quickly, and are therefore executed:
`check_sampling()`, `screen_models()`, `check_fit()`, `curve_params()`,
`bnec_record()`, `pull_best()`, `failed_models()`, `dispersion()` and
`get_priors()`. Functions that need no fit at all, and are therefore executed:
`show_params()`, `models()` and the `spacing_cv()` helper. Each was timed on the
module's own fit before being written in; the timings are in
`prompts/` alongside the session log.

`compare_estimates()`, `ecnsec()` and `average_estimates()` are not executed on a
model-averaged set, because each rebuilds the component predictions.

## Saved fits for participants

The run-time limit above was read too narrowly at first. It was applied to the run time of a
render, and the user's instruction covered the room as well: a participant
running the code alongside the presenter should not wait for the sampler
either. `_freeze/` does nothing for that person, because it keeps the
published page from refitting and not their session.

Measured on 2026-09-16, rendering the development profile on a four-core WSL2
machine with `cmdstanr`:

| module | wall time | fits |
|---|---|---|
| 2 | 1.5 min | 1 |
| 3 | 7 s | 0 |
| 4 | 4.5 min | 3 |
| 5 | 12.9 min | 8 |
| 6 | over 15 min | 6, at `iter = 1e4` and `adapt_delta = 0.99` |

Thirteen minutes for one module is the failure the instruction described.

The modules show the workflow instead of hiding it. Each fit is three steps:

```r
# shown, eval: false
set.seed(333)
exp_1nec <- bnec(suc | trials(tot) ~ crf(raw_x, model = "nec3param"),
                 data = binom_data)
save(exp_1nec, file = "fits/m5_exp_1nec.RData")
```

```r
# executed
load("fits/m5_exp_1nec.RData")
```

A `cached()` helper was written first, wrapping each call so that it fitted or
loaded as needed, and was withdrawn on the user's judgement. It worked. It was
the wrong answer, because it invented a course-specific abstraction for the one
thing the participant should be learning to do for themselves, and because
fitting once and loading thereafter is what the presenter does in their own
work. A module that shows it teaches something beyond the model.

Naming. A file is `fits/m<N>_<object>.RData`, the module number being needed
because `bnec_fit` is the object name in modules 2, 4 and 6 and module 2's call
differs from the other two in passing a `seed`. Without the prefix those would
collide with different content under one name.

Reproducibility. A `set.seed()` inside a fit chunk stops running once that chunk
is `eval: false`, which would leave later random steps unseeded --
`sample_priors()` in module 6 is not internally seeded and was the case that
surfaced this. The seed stays in the fit chunk where it belongs in the workflow,
and each module's setup chunk gained a `set.seed(333)` so that the document is
reproducible without any line appearing twice. Module 5 already had one.

Module 4 also gained `data(nec_data)` in its setup, that call having lived only
inside a fit chunk.

Two consequences. A module no longer samples at render, so a rebuild is minutes
rather than half an hour. And re-executing a module from scratch now requires the
saved objects, exactly as module 7 already required its generated results, which
makes the distribution below part of the method rather than a convenience.

### The wiring as applied

All eighteen executing fit calls across modules 2, 4, 5 and 6 are wrapped, and
`source("fit_cache.R")` was added to each setup chunk. Module 2 explains
`cached()` in a subsection beside the existing `saveRDS()` material, which is
where a reader first meets it and where the argument for saving a fit is already
made.

The wrapping was applied by a script rather than by hand, because a regular
expression cannot do it safely: a fit call may span several lines, and five of
the twenty-three matches sit inside `eval: false` demonstrations that must not
be touched. The script parses the chunks, skips the ones marked `eval: false`,
and finds each call's closing parenthesis by counting depth across lines.

Running it dry first was what made it correct. The first version inserted the
closing parenthesis but not the opening `cached(` wherever a call sat on a
single line, because the second edit was computed against the original string
and overwrote the first. That is invisible in a diff summary and would have
produced files that do not parse.

### Distribution

A release asset on the course repository, with a USB stick as the fallback on
the day. The reasoning, and the alternatives rejected, are in the human
document. The tag is `fits`, the asset is `cr_modelling_fits.zip`, and the
address is

```
https://github.com/open-AIMS/cr_modelling_training/releases/download/fits/cr_modelling_fits.zip
```

Rebuild with `Rscript scripts/bundle_fits.R`, then
`gh release upload fits dist/cr_modelling_fits.zip --clobber -R open-AIMS/cr_modelling_training`.
The tag does not change, so the address holds and nothing else needs editing.
Commit the rewritten `vignettes/fits.sha256` in the same change as the upload,
since that file is what a participant's download is verified against.

Four supporting pieces. `scripts/bundle_fits.R` zips `fits/` and writes a
SHA-256 and a provenance file to `dist/`; the checksum matters because a
truncated archive has a plausible size and fails later and obscurely.
`vignettes/fetch_fits.R` downloads, verifies and unpacks, and falls back to fitting
the models on any failure. `check_setup.R` gained Stage 7, recorded
through `record("fits", NA, ...)` so that a missing archive is a warning rather
than a failure, matching its optional status. And `0Software-setup.qmd` gained
a section beside the existing clone and ZIP instructions.

Zenodo is deferred rather than rejected. It gives a citable identifier and is
the right destination once the material is frozen; it is not the route for the
workshop, because the archive is rebuilt as modules are revised and a Zenodo
record is versioned rather than replaced in place.

`scripts/zenodo_deposit.R` creates the deposition, uploads the archive and sets
the metadata, then stops. It does not publish, because a published Zenodo record
cannot be withdrawn, only superseded. It takes `--sandbox` for rehearsal against
`sandbox.zenodo.org`, which uses separate tokens. `httr` and `jsonlite` are both
installed.

The address lives in one constant, `FITS_URL` at the top of
`fetch_fits.R`, with a `CR_FITS_URL` environment override for testing before any
deposit exists. `0Software-setup.qmd` names the command rather than the address
and has no R chunks at all, so setting the real identifier is one line and a
render measured in seconds.

## Pre-computed results

Module 7 established the pattern and module 8 adopts it.
`scripts/generate_herbicide_fits.R` fits the seven herbicides deliberately,
outside any render, and writes compact results to `vignettes/data/`. The module
shows the fitting code with `eval: false` and reads the saved results. The fitted
objects are never committed: the 2023 equivalent was 451 MB and the repository is
public.

`vignettes/data/` does not exist in the working tree, so module 7 cannot be
rendered as it stands. The script has not been run. This has to be resolved
before module 7 is added to either render list, and it is the reason module 7 is
absent from `_quarto.yml` and `_quarto-dev.yml`.

A generating script states, in its header: what it fits, why it is not run at
render time, how long it took, and the command that runs it. It writes only what
the module reads.

## The Stan backend is not set by the project

Found on 2026-09-16 while rendering. `notes/setup-evidence.md` records the
decision that the sampling backend is `cmdstanr`, set with
`options(brms.backend = "cmdstanr")`. Nothing in the project applies it. There is
no `.Rprofile` at the repository root, and the only module mentioning
`brms.backend` is module 1, where it appears in an `eval: false` illustration.
Every module's setup chunk sets `options(mc.cores = 4)` and none sets the
backend.

A plain `quarto render` therefore builds the site under `rstan`, which is the
`brms` default, and not under the backend the setup module tells participants to
install. The render made on 2026-09-16 forced it through `R_PROFILE_USER`
pointing at a one-line profile, which is a workaround for one render rather than
a fix.

Two ways to fix it. A project `.Rprofile` holding
`options(brms.backend = "cmdstanr", mc.cores = 4)` applies to every render and
every interactive session in the project, and is the smaller change. Adding the
option to each module's setup chunk beside the existing `mc.cores` line is more
explicit and matches the style already there, at the price of repeating it in
seven files.

The first is preferable, subject to one caveat that belongs in the file itself.
A participant who opens a module outside the project does not get the
`.Rprofile`, so module 1's instruction to set the backend stays where it is.

## Clearing the stored output

`execute: freeze: auto` keys on the source document rather than on the
environment. Any module whose chunks change needs its stored result cleared:

```bash
rm -rf _freeze/vignettes/<module>
```

Module 4 needs this whether or not its chunks change, because the summary output
it quotes was rendered under the previous weighting default.

## Bibliography additions

To be added to `vignettes/bayesnec.bib`, taken from
`C:/Rworking/bayesnec/vignettes/bayesnec.bib` where an entry exists there:

| key | paper | needed by |
|---|---|---|
| `Vehtari2021` | Vehtari et al. (2021), on R-hat | sampler diagnostics, screening |
| `Ritz2026` | Ritz et al. (2026), on control division | preparing the response |
| `helsel2006` | Helsel (2006), on non-detects | censoring |
| `Weimer2012` | Weimer et al. (2012), on transformations | preparing the response, fixing a parameter |
| `oecd2026tg201` | OECD Test Guideline 201 (2026) | growth rates |
| `warton2005` | Warton (2005), on excess zeros | zeros |
| `martin2005` | Martin et al. (2005), on sources of zeros | zeros |
| `Krull2020` | Krull (2020), on threshold estimation | experimental design |
| `warne2025` | Warne et al. (2025), the ANZ method | guideline position |

`warne2025` and `Krull2020` have no entry in the `bayesnec` bibliography and must
be written. `fisheretal2023`, already in the course bibliography, is the IEAM
paper and is the key to use for the N(S)EC evidence; it duplicates
`fisher2023ieam` in the `bayesnec` bibliography.

## Register

Repository CLAUDE.md §7 governs. Second person in body text, first person plural
where the package authors are speaking, noun-phrase headings, no congratulation
in solution text, and no statement that anything is easy. The 2023 text is the
register reference, so read the surrounding paragraphs before writing a new one.

Every addition above is a claim about package behaviour, so repository CLAUDE.md
§7 applies in full: run the code, or cite the help page or the paper. Where a
claim came from a vignette on `dev` rather than from the installed package, it is
marked as such here and must be checked before it is written into a module.

## Revising the course onto shipped datasets

Decided with the user on 2026-09-16, after `open-AIMS/bayesnec#228` was raised.
The course holds five data files of its own, three of which duplicate the
teaching purpose of a dataset that either already ships with `bayesnec` or will
once that pull request lands. Using the shipped data instead removes
undocumented files from a public repository and gives every example a manual
page a participant can read.

### Available in the installed package

`data(package = "bayesnec")` in 2.1.3.35 returns `alga`, `herbicide`,
`manec_example`, `nassarius` and `nec_data`.

`manec_example` was used for a failing-screen demonstration and then withdrawn,
on the user's correction. It is a test fixture rather than an analysis: 100
retained draws, so its R-hat of 1.213 and effective sample size of 10 are
properties of how it was built and not diagnostics anyone would meet. Teaching a
screen from it would say that a screen fails when the object is small, which is
the wrong lesson.

The replacement is the pair of datasets the `bayesnec` workflow article uses,
both of which ship: `herbicide`, of which `simazine` is one series, and
`nassarius`. That article reports three outcomes from the same screen, and module
4 now gives all three, attributed and shown with `eval: false` because fitting
three model sets is well past the run-time limit. The figures are 5 of 11
equations removed on `simazine` at reduced settings; 1 of 10 at `bnec()`'s
defaults on the same data, three of the five having passed once they had more
draws, with `ecxll5` then failing to fit at all and `ecxwb2` still failing on
divergent transitions; and 8 of 10 on `nassarius` at those same defaults, on
divergence counts of 120 to 5,589. The `nassarius` estimate over the two
survivors is 1.46 with an interval of 0.31 to 1.92, which is the untested gap
between 1.25 and 2.5 in that design.

That set of three is what the section needed. It shows all three kinds of failure
on real data, and the `nassarius` result makes the point that a screen outcome
can be a statement about the experiment.

`alga` holds the below-limit counting records that module 5's censoring section
currently quotes from the `bayesnec` documentation. The section could compute
that three-way comparison instead of quoting it. It is three fits, so it belongs
in a generating script rather than in the module.

`nassarius` is survival counts out of a known number exposed, which is the
binomial and beta-binomial case module 5 currently teaches from
`example_binomial.csv`. `herbicide` is a continuous proportion on the unit
interval, which is the `Beta` case module 5 teaches from
`example_proportion.csv`. Both shipped datasets are documented, published and
already used elsewhere in the course.

### Arriving with the grouping pull request

`open-AIMS/bayesnec#228` adds `coral_colour`, `coral_pam` and `lum31`. Each was
built for the grouping vignette and each matches one of module 8's own files.

`coral_colour` is 180 rows: three climate scenarios by five diuron
concentrations by three chambers, four fragments per chamber. A chamber sits at
one concentration, so it admits a displacement and nothing else, which is
exactly the within-concentration grouping module 8 teaches from
`example_ogl.csv`. Its manual page also records that five concentrations is at
the conventional minimum and that four-parameter equations do not sample
reliably against it while three-parameter equations do, which is the design point
module 3 makes in prose and cannot presently demonstrate.

`coral_pam` is 414 rows over 54 chambers with repeated readings, and its manual
page states it is the better conditioned of the two for an `ogl(chamber)` term.
Sixty-three readings are exactly 0, measured rather than substituted, which is
the boundary case module 5 describes.

`lum31` is 2904 rows and replaces `example_pgl.csv`, which holds the same assay
in undocumented form. It supplies replication at three nested scales in one
dataset: four wells within a concentration, 33 plates each spanning a whole
dilution series, and two toxicants read at two exposure times. One dataset
therefore covers the within-concentration term, the across-concentration term
and the factor covariate, which module 8 currently needs three files for. It
also has `censoring` and `rlu_cens` columns, because 386 of its readings are
blank-corrected negatives floored at zero, so it is a worked censoring case as
well.

### The order of the substitutions

Module 8 is the substantial one. Replace `example_ogl.csv` with `coral_colour`
or `coral_pam` for the displacement example, and `example_pgl.csv` with `lum31`
for the whole-curve example. The factor-covariate section can then use `lum31`'s
`toxicant` or `minutes`, or `coral_pam`'s `climate`, rather than
`example_fi.csv`. Doing so removes the one remaining dataset in the course whose
source cannot be cited, since `example_fi.csv` is identified only by a DOI and
the name of the person who supplied it.

Module 5's binomial and proportion sections can use `nassarius` and `herbicide`
independently of the pull request, and its censoring section can compute the
`alga` comparison rather than quoting it.

Three cautions. Every substitution changes the numbers the prose reads, so the
prose has to be rewritten against the new output rather than reused. The
three new datasets are larger than the files they replace, and `lum31` at 2904
rows will fit more slowly than `example_pgl.csv`, which matters under the
run-time limit above. And none of the three is in the installed package, so
nothing may be written against them until the pull request merges and the
package is upgraded.

## The dispersion sub-model and the case against a divisor

`open-AIMS/bayesnec#372`, "Confirm the dispersion sub-model with a posterior
check", adds a third worked response to `example9`: one plate of the `lum31`
copper series at 15 minutes, 44 wells over 11 concentrations at four replicate
wells. It targets `#228` rather than `dev` and cannot merge before it, because
`lum31` is only on that branch. The user referred to it as 376, which does not
exist; 372 is the open pull request matching the description.

Module 5's non-constant dispersion section now gives its measurement rather than
arguing the direction in the abstract. Under constant dispersion `check_fit()`
flags 7 of 11 groups on spread and the `sd_ratio` at the control is 0.09, so the
fit simulates about eleven times the spread the data show there. Under
`disp("power")` no group is flagged and the ratio is 1.02. The *NSEC* halves,
from 0.105 to 0.055 mg/L, because a control posterior that is too wide places the
significance reference too low and a reference too low is crossed further along
the curve. The **EC~x~** changes from 0.083 to 0.070, being measured at a fixed
fraction of the control mean rather than at a quantile of its posterior, though
its interval narrows from 0.055 to 0.011. The estimated exponent is 0.78 against
0.77 from the replicate wells with no model involved, which is the evidence that
the sub-model describes measured variation rather than absorbing a bad curve.

The connection to @Ritz2026 is the part that needed writing down, and module 2
now states it. A response on an arbitrary scale is often divided by its observed
maximum so that it falls between 0 and 1 and a `Beta` can be used, because a
`Beta` accommodates variability that changes with the mean where a `gaussian`
does not. That identifies a real problem and answers it indirectly: it pays for a
variance assumption with a bias in the mean, and it forces one observation onto
the boundary of the `Beta` support by construction. Keeping the response on its
own scale, choosing the family that suits it, and modelling the changing
variability with `disp` addresses the same problem directly, keeps the parameters
in measured units, and states the change in variability as something estimated
rather than absorbed.

This is the practice the group used to follow, so module 2 states it as a
correction rather than as a general principle.

## The source of the saved fits

`scripts/generate_taught_fits.R`, and nothing else. This follows from the design, and it
is recorded here because it inverted a rule that had been correct.

Every fit chunk is now `eval: false`, so a render samples nothing. The objects
therefore cannot come from a render, and the script is the only thing that makes
them. Its chunk-selection rule had been "run the executed chunks, skip every
`eval: false` chunk", which was right under the previous design and, after the
change, excluded precisely the chunks that fit. It would have run to completion
and reported success having fitted nothing.

The rule is now: run an executed chunk, and run an `eval: false` chunk where it
saves into `fits/`. Those are the fit chunks and they are the point of the
script.

One refinement was needed beyond that. Module 2 illustrates the fit-save-load
pattern with a chunk that repeats its own fit call verbatim, including the
`save()`, so the rule selected nineteen fits where the modules contain eighteen
and module 2's model would have been fitted twice. The selection now skips a
chunk whose save targets have all been written already in the same run, which
counts 1, 3, 8 and 6 across the four modules.

## The location of the saved fits

`vignettes/fits/`, not the repository root. The modules execute with
`vignettes/` as the working directory, which is why module 5 reads
`example_binomial.csv` and module 7 reads `data/herbicide_sampling.csv` without
a prefix. A `load("fits/...")` in a module therefore resolves to
`vignettes/fits/`, and putting the objects at the repository root would have
failed at the first render. `.gitignore` holds `vignettes/fits/`.

The archive was 23.2 MB over eighteen objects when first built, against an
estimate of 100 to 200 MB made before anything was built. The estimate was wrong
by roughly an order of magnitude, and the consequence is that the distribution
question is much less constrained than it appeared. It is 42.3 MB over 25
objects as of 2026-09-17, and still an unremarkable download: an end-to-end run
of `source("vignettes/fetch_fits.R")` took 32 seconds on that date. The count
grows as modules are revised, so read it from `dist/cr_modelling_fits.txt` and
reset `N_FITS` in `fetch_fits.R` after every rebuild.

The zip is reproducible. Rebuilding on 2026-09-17 from unchanged objects gave
the same sha256 as the previous build, so a rebuild that changes nothing does
not invalidate a copy a participant already holds.

`scripts/bundle_fits.R` had a path defect on its first run, caught by the
overnight chain rather than by review. It changed into the parent of the fits
directory before zipping, in order to store the archive's paths relative to the
repository root; but that parent *is* the repository root, which was already the
working directory, so the archive was written to `../dist/` outside the
repository and the run failed with `zip I/O error`. The directory change was
unnecessary and is removed.

## Pending upstream changes

`open-AIMS/bayesnec#225` adds `hurdle_poisson` and `hurdle_negbinomial` as
two-block families, the count analogues of `hurdle_gamma`. It was a draft parked
behind `paul-buerkner/brms#1923`, was unblocked on 2026-09-16, and is expected to
merge. It is not in the installed 2.1.3.35, whose `mod_fams` holds `gaussian`,
`Gamma`, `poisson`, `negbinomial`, `bernoulli`, `binomial`, `beta_binomial`,
`beta`, `hurdle_gamma`, `zero_inflated_beta`, `zero_inflated_poisson` and
`zero_inflated_negbinomial`.

Module 5's zeros section was written against the installed set, and prints
`names(bayesnec:::mod_fams)` rather than listing the families in prose, so the
list updates itself on the next render once the package is upgraded. Two
paragraphs of that section are version-independent and were written from the pull
request, because they are conceptual rather than a claim about what is exported:
the observed-against-latent distinction between a hurdle and zero inflation, and
the zero-truncation trap.

What to add once it merges. Name the two count families where the section names
`hurdle_gamma`, since a count response with structural zeros is the common
ecotoxicological case, a fecundity count where the animal died. State that
`bnec_hurdle()` refuses `poisson` and `negbinomial` growth families, and why: it
fits the positive part to `data[y > 0, ]` with an untruncated count family, which
estimates `mu / (1 - exp(-mu))` rather than `mu`, and the bias grows as the mean
falls towards zero, which is the end an *NEC* is read from. The joint route,
`bnec(family = "hurdle_poisson")`, needs no truncation because `brms` writes the
zero-truncated positive part itself. The pull request reports a verification on
simulated data with 47% structural zeros: a true *NEC* of 3 was estimated at
3.049 with a 95% interval of 2.536 to 3.607.

The stack `#226` to `#227` to `#228` and `#238` sits above `#225`, and `#228` is
the grouping vignette module 8 draws on, so module 8 should be re-checked against
`example8` after that stack lands.

## Prompt logging

Repository CLAUDE.md §10 applies. These changes alter what the modules teach
about method, so they are logged to `prompts/`, which does not yet exist.
