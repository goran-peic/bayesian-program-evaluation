## =============================================================================
##  Bayesian Analysis of a Factorial Experiment on Program Information Design
##  -----------------------------------------------------------------------------
##  Question. A state agency wants parents to choose higher-quality child-care
##  programs. It can present program information in different ways. Which
##  combination of presentation choices most increases the share of parents who
##  pick the higher-quality option?
##
##  Design. A 2 x 2 x 2 factorial experiment run across 24 sites. Each parent is
##  randomized to one of eight information formats and makes a binary choice.
##    cost_format : "annual total" vs "monthly"            (how price is shown)
##    quality_cue : "star rating" vs "star rating + text"  (how quality is shown)
##    comparison  : "single program" vs "side-by-side"     (how options are shown)
##
##  Analysis. A hierarchical logistic model with all main effects and two-way
##  interactions, site-varying intercepts, and site-varying slopes for the
##  comparison factor. Because the data are simulated from known parameters, the
##  script checks that the model recovers them before interpreting anything.
##
##  Output. Posterior probabilities that each factor helps, the expected effect
##  of each of the eight formats, and a deployment recommendation with its
##  uncertainty: the deliverable a program office can act on.
##
##  Reproducible: set.seed() is called once; no files are read.
##  Requires: brms (>= 2.20), rstan or cmdstanr, dplyr, tidyr, ggplot2, posterior.
##  Last Tested on R-4.6.1.
## =============================================================================

## =============================================================================
# Quick Install of Any Missing Packages
# ------------------------------------------------------------------------------
# install.packages("brms") # Uncomment if your local R run is missing brms.
# install.packages("tidyr") # Uncomment if your local R run is missing tidyr.
# install.packages("cmdstanr",
#                  repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
# library(cmdstanr)
# check_cmdstan_toolchain(fix = TRUE)   # verifies Rtools45 is wired up correctly
# install_cmdstan(cores = 4)            # downloads and builds CmdStan once
# cmdstan_version()                     # prints a version number if all is well
## =============================================================================

rm(list=ls())

suppressPackageStartupMessages({
  library(brms)
  options(brms.backend = "cmdstanr") 
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
})

set.seed(09021980)
# setwd("C:/Users/goran/downloads") # Set custom, local output path for plots.
out_dir <- file.path(getwd(), "output")
dir.create(out_dir, showWarnings = FALSE)

## ---- 1. Simulate the experiment ---------------------------------------------
## True parameters on the logit scale. The point of simulating is to know these,
## so that section 4 can test whether the model recovers them.

truth <- c(
  intercept        = -0.40,  # baseline: 40% choose the higher-quality program
  cost_monthly     =  0.35,  # monthly price framing helps
  quality_text     =  0.25,  # adding text to the star rating helps a little
  comparison_side  =  0.60,  # side-by-side comparison helps the most
  cost_x_compare   = -0.30,  # ...but monthly framing adds less when side-by-side
  sd_site          =  0.35,  # sites differ in their baseline rate
  sd_site_compare  =  0.25   # ...and in how much side-by-side helps them
)

n_sites        <- 24
n_per_site     <- 160                     # 3,840 parents in total

sites <- tibble(
  site      = factor(seq_len(n_sites)),
  u_site    = rnorm(n_sites, 0, truth["sd_site"]),
  u_compare = rnorm(n_sites, 0, truth["sd_site_compare"])
)

dat <- expand_grid(site = sites$site, rep = seq_len(n_per_site)) |>
  # Independent randomization of each factor, as in the field design.
  mutate(
    cost_format = sample(c("annual", "monthly"), n(), replace = TRUE),
    quality_cue = sample(c("stars", "stars_text"), n(), replace = TRUE),
    comparison  = sample(c("single", "side_by_side"), n(), replace = TRUE)
  ) |>
  left_join(sites, by = "site") |>
  mutate(
    cost_m = as.integer(cost_format == "monthly"),
    qual_t = as.integer(quality_cue == "stars_text"),
    comp_s = as.integer(comparison  == "side_by_side"),
    eta = truth["intercept"] + u_site +
      truth["cost_monthly"]    * cost_m +
      truth["quality_text"]    * qual_t +
      (truth["comparison_side"] + u_compare) * comp_s +
      truth["cost_x_compare"]  * cost_m * comp_s,
    chose_quality = rbinom(n(), 1, plogis(eta))
  ) |>
  select(site, cost_format, quality_cue, comparison, chose_quality) |>
  # Fix reference levels explicitly (the status quo format) so that each
  # coefficient reads as "effect of switching away from the current practice".
  mutate(
    cost_format = factor(cost_format, levels = c("annual", "monthly")),
    quality_cue = factor(quality_cue, levels = c("stars", "stars_text")),
    comparison  = factor(comparison,  levels = c("single", "side_by_side"))
  )

cat("Simulated", nrow(dat), "parents across", n_sites, "sites.\n")
cat("Raw choice rate by format:\n")
print(dat |> count(cost_format, quality_cue, comparison, wt = chose_quality,
                   name = "n_quality") |>
        left_join(count(dat, cost_format, quality_cue, comparison),
                  by = c("cost_format", "quality_cue", "comparison")) |>
        mutate(rate = round(n_quality / n, 3)) |> select(-n_quality))

## ---- 2. Priors ---------------------------------------------------------------
## Weakly informative priors on the logit scale. A normal(0, 1) prior on a main
## effect says an odds ratio outside roughly [0.14, 7.4] is implausible for an
## information nudge; this rules out nothing a program office would believe.
## The prior on the varying-effect SDs keeps site differences moderate.

priors <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 1),   class = "b"),
  prior(exponential(2), class = "sd")
)

## A quick prior predictive check: what choice rates do these priors imply?
## If the implied rates were piled up near 0 or 1, the priors would be too wide.

model_formula <- bf(
  chose_quality ~ (cost_format + quality_cue + comparison)^2 +
    (1 + comparison | site)
)

fit_prior <- brm(
  model_formula, data = dat, family = bernoulli(), prior = priors,
  sample_prior = "only", chains = 2, iter = 1000, seed = 1, refresh = 0,
  silent = 2
)

prior_rates <- posterior_epred(fit_prior, ndraws = 200)
cat(sprintf(
  "\nPrior predictive choice rates: median %.2f, 90%% interval [%.2f, %.2f]\n",
  median(prior_rates), quantile(prior_rates, .05), quantile(prior_rates, .95)
))

## ---- 3. Fit the model --------------------------------------------------------

fit <- brm(
  model_formula, data = dat, family = bernoulli(), prior = priors,
  chains = 4, iter = 2000, warmup = 1000, seed = 2, refresh = 0, silent = 2,
  control = list(adapt_delta = 0.9)
)

## Convergence: every Rhat should be ~1.00 and bulk ESS comfortably > 400.
diag <- summarise_draws(as_draws_df(fit), "rhat", "ess_bulk") |>
  filter(!grepl("^r_|^lp__|^z_|^L_", variable))
stopifnot(all(diag$rhat < 1.01), all(diag$ess_bulk > 400))
cat("\nConvergence OK: max Rhat =", round(max(diag$rhat), 3),
    "; min bulk ESS =", round(min(diag$ess_bulk)), "\n")

## ---- 4. Parameter recovery ---------------------------------------------------
## The model's population-level coefficients should bracket the values the data
## were generated from. This is the test that separates a model that works from
## one that merely runs.

recovery <- fixef(fit, probs = c(.05, .95)) |>
  as.data.frame() |>
  tibble::rownames_to_column("term") |>
  mutate(truth = c(
    truth["intercept"], truth["cost_monthly"], truth["quality_text"],
    truth["comparison_side"],
    0,                          # cost x quality: no true interaction
    truth["cost_x_compare"],
    0                           # quality x comparison: no true interaction
  )) |>
  mutate(recovered = truth >= Q5 & truth <= Q95)

cat("\nParameter recovery (90% intervals should contain the truth):\n")
print(recovery |> mutate(across(where(is.numeric), \(x) round(x, 2))))
cat(sprintf("%d of %d population-level parameters recovered; with 90%% intervals, about %.1f misses are expected by chance.\n",
            sum(recovery$recovered), nrow(recovery), 0.1 * nrow(recovery)))

sd_recovery <- VarCorr(fit)$site$sd |> as.data.frame() |>
  tibble::rownames_to_column("term") |>
  mutate(truth = c(truth["sd_site"], truth["sd_site_compare"]))
cat("\nSite-level SDs:\n")
print(sd_recovery |> mutate(across(where(is.numeric), \(x) round(x, 2))))

## ---- 5. Posterior predictive check -----------------------------------------
## Does the fitted model reproduce the observed choice rate in each of the eight
## cells? Plotted as observed vs. predicted with 90% intervals.

cells <- dat |> distinct(cost_format, quality_cue, comparison)
ppc <- posterior_epred(fit, newdata = cells, re_formula = NA) # population level
cells_pred <- cells |>
  mutate(
    pred_mean = colMeans(ppc),
    pred_lo   = apply(ppc, 2, quantile, .05),
    pred_hi   = apply(ppc, 2, quantile, .95)
  ) |>
  left_join(
    dat |> group_by(cost_format, quality_cue, comparison) |>
      summarise(observed = mean(chose_quality), .groups = "drop"),
    by = c("cost_format", "quality_cue", "comparison")
  ) |>
  mutate(format = paste(cost_format, quality_cue, comparison, sep = " / "))

p_ppc <- ggplot(cells_pred, aes(x = reorder(format, pred_mean))) +
  geom_pointrange(aes(y = pred_mean, ymin = pred_lo, ymax = pred_hi),
                  colour = "grey40") +
  geom_point(aes(y = observed), colour = "firebrick", size = 2.5) +
  coord_flip() +
  labs(x = NULL, y = "Share choosing the higher-quality program",
       title = "Posterior predictive check by information format",
       subtitle = "Grey: model (90% interval, population level). Red: observed.") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "01_ppc_by_format.png"), p_ppc, width = 8, height = 4.5,
       dpi = 150)

## ---- 6. What the program office needs ---------------------------------------
## (a) Probability that each factor increases quality choices, at the average
##     level of the other factors. Averaging over the other factors matters when
##     interactions are present; a main-effect coefficient alone does not answer
##     the policy question.

epred_cells <- as_draws_matrix(ppc)      # draws x 8 cells, population level
colnames(epred_cells) <- cells_pred$format

effect_of <- function(factor_col, high, low) {
  hi <- rowMeans(epred_cells[, cells[[factor_col]] == high, drop = FALSE])
  lo <- rowMeans(epred_cells[, cells[[factor_col]] == low,  drop = FALSE])
  hi - lo
}

effects <- tibble(
  factor = c("Monthly price framing", "Star rating + text", "Side-by-side comparison"),
  draws  = list(
    effect_of("cost_format", "monthly", "annual"),
    effect_of("quality_cue", "stars_text", "stars"),
    effect_of("comparison",  "side_by_side", "single")
  )
) |>
  mutate(
    mean_pp   = sapply(draws, mean),
    lo_90     = sapply(draws, quantile, .05),
    hi_90     = sapply(draws, quantile, .95),
    pr_helps  = sapply(draws, \(d) mean(d > 0)),
    pr_gt_3pp = sapply(draws, \(d) mean(d > 0.03))   # a practically meaningful gain
  ) |>
  select(-draws)

cat("\nAverage effect of each factor (percentage points, averaged over the others):\n")
print(effects |> mutate(across(where(is.numeric), \(x) round(x, 3))))

## (b) Which of the eight formats to deploy, and how confident we can be.
best_idx <- apply(epred_cells, 1, which.max)
deploy <- cells_pred |>
  mutate(pr_best = as.numeric(table(factor(best_idx, levels = 1:8))) / nrow(epred_cells)) |>
  arrange(desc(pred_mean)) |>
  select(format, pred_mean, pred_lo, pred_hi, pr_best)

cat("\nFormats ranked by expected choice rate (pr_best = posterior probability of being the best):\n")
print(deploy |> mutate(across(where(is.numeric), \(x) round(x, 3))))

## (c) The interaction that matters in practice: monthly framing helps a lot
##     for single-program displays and much less once programs are side by side.

gain_monthly_single <- epred_cells[, "monthly / stars / single"] -
                       epred_cells[, "annual / stars / single"]
gain_monthly_side   <- epred_cells[, "monthly / stars / side_by_side"] -
                       epred_cells[, "annual / stars / side_by_side"]

cat(sprintf(
  "\nGain from monthly framing: %.1f pp with single display, %.1f pp side by side (P(difference > 0) = %.2f)\n",
  100 * mean(gain_monthly_single), 100 * mean(gain_monthly_side),
  mean(gain_monthly_single > gain_monthly_side)
))

## (d) Site heterogeneity: how much does the best format's benefit vary by site?
site_eff <- ranef(fit)$site[, , "comparisonside_by_side"] |> as.data.frame() |>
  tibble::rownames_to_column("site")
p_sites <- ggplot(site_eff, aes(x = reorder(site, Estimate), y = Estimate,
                                ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(colour = "grey30") +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  labs(x = "Site", y = "Site deviation in the side-by-side effect (logit)",
       title = "Does side-by-side comparison help every site?",
       subtitle = "Partial-pooling estimates with 95% intervals") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "01_site_heterogeneity.png"), p_sites, width = 7,
       height = 5, dpi = 150)

## ---- 7. Plain-language summary ------------------------------------------------
## Generated from the posterior, not typed by hand, so it cannot contradict the
## numbers above when the data change.

describe <- function(e) {
  size <- if (e$mean_pp >= 0.05) "a clear" else if (e$mean_pp >= 0.02) "a modest" else "a negligible"
  conf <- if (e$pr_helps >= 0.95) "" else if (e$pr_helps >= 0.80) ", though with some uncertainty" else ", but the evidence is weak"
  sprintf("%s gives %s gain of %.0f points (P(helps) = %.2f%s)", e$factor, size,
          100 * e$mean_pp, e$pr_helps, conf)
}

ranked <- effects |> arrange(desc(mean_pp))
top <- deploy[1, ]
status_quo <- cells_pred$pred_mean[cells_pred$format == "annual / stars / single"]

cat("\n--- Summary for the program office ---\n")
cat(sprintf(
  paste0("%s; %s; %s. Monthly framing and side-by-side display overlap: monthly ",
         "framing adds %.0f points when programs are shown one at a time but %.0f ",
         "points once they are side by side (if negative, it slightly hurts). The recommended format is '%s' ",
         "(posterior probability of being best: %.2f): an expected %.0f%% of parents ",
         "choosing the higher-quality program (90%% interval %.0f-%.0f%%), against ",
         "%.0f%% under the current annual / stars / single format. Sites differ in ",
         "how much they gain, so a staged rollout with site-level monitoring is advisable.\n"),
  describe(ranked[1, ]), describe(ranked[2, ]), describe(ranked[3, ]),
  100 * mean(gain_monthly_single), 100 * mean(gain_monthly_side),
  top$format, top$pr_best, 100 * top$pred_mean, 100 * top$pred_lo, 100 * top$pred_hi,
  100 * status_quo
))

sessionInfo()
