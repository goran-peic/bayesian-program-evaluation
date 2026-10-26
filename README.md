# Bayesian methods for program evaluation: two R code samples

Two self-contained R scripts, each built around a question a program office or payer actually asks, answered with a hierarchical Bayesian model fit in Stan through `brms`. Both scripts simulate their data from known parameters, so each one checks that the model recovers the truth before interpreting anything. Both end with a decision, stated in plain language and generated from the posterior rather than typed by hand.



|Script|Question|Method|Runtime|
|-|-|-|-|
|`01\_factorial\_information\_design.R`|Which way of presenting program information most increases the share of parents who choose the higher-quality option?|2 x 2 x 2 factorial experiment across 24 sites; hierarchical logistic model with two-way interactions and site-varying effects; prior and posterior predictive checks; deployment recommendation with its uncertainty|about 5 minutes|
|`02\_practice\_ranking\_shrinkage.R`|Which primary care practices are really in the top decile on a noisy quality measure, and how sure can we be?|Hierarchical binomial model with risk adjustment; ranking by posterior probability of being top-decile; scorecard against two naive rankings, with the truth known by simulation|about 2 minutes|

## 

## What to look at

**Sample 1.** Section 4 is the parameter-recovery table: six of seven population-level parameters fall inside their 90% intervals, with the one miss being an interaction whose true value is zero (about 0.7 misses are expected by chance). Section 6 shows why averaging over the other factors matters when interactions are present: monthly price framing adds about six points when programs are shown one at a time and nothing once they are side by side. The final paragraph is the recommendation a program office could act on, including the posterior probability (0.84) that the recommended format is in fact the best of the eight.

**Sample 2.** Section 6 is the scorecard. Ranking practices on their raw rate rewards 15 practices of which 8 belong in the top decile; the mistaken awards go to small practices with lucky years. A minimum-size rule, the usual fix, does worse, because it throws out small practices that are genuinely excellent. The hierarchical ranking gets 10 of 15 right and, more usefully, says which awards are confident and which sit in a grey zone where another year of data is worth more than a decision now. `output/02\_pr\_top\_decile.png` is the figure a payer would use to set the reward threshold.

## 

## Running the scripts

The scripts use `brms` with the `cmdstanr` backend, which compiles models with CmdStan. On Windows this needs Rtools45.

```r
install.packages(c("brms", "dplyr", "tidyr", "ggplot2", "posterior"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install\_cmdstan()      # one-time build of CmdStan

source("01\_factorial\_information\_design.R")
source("02\_practice\_ranking\_shrinkage.R")
```

Each script sets `options(brms.backend = "cmdstanr")` at the top and writes its figures to `output/`. Tested with R version 4.6.1, Windows 11, brms version 2.23.0, CmdStan version 0.9.0.



## Design choices worth a word

* **Simulated data, on purpose.** Real evaluation data cannot be shared, and a simulation is the only setting in which a method can be checked against a known answer. Each script states its true parameters at the top.
* **Weakly informative priors, stated and checked.** A prior predictive check in sample 1 confirms the priors imply plausible choice rates rather than piling up at 0 or 1.
* **Convergence is asserted, not eyeballed.** Each script stops if any R-hat exceeds 1.01 or any bulk effective sample size falls below 400.
* **The deliverable is a decision.** Posterior probabilities that a factor helps, that a format is the best, that a practice belongs in the top decile. Point estimates and p-values are not what a policymaker needs.
* **Models are deliberately simple.** The point of both samples is the inferential framework and the honest accounting of uncertainty, not model complexity.

