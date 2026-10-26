## =============================================================================
##  Identifying Top-Performing Practices when the Measure is Noisy
##  -----------------------------------------------------------------------------
##  Question. A payer wants to reward the best-performing primary care practices
##  on a quality measure (say, the share of diabetic patients with controlled
##  blood sugar). Practices differ enormously in size: a 20-patient practice can
##  post 95% by luck; a 900-patient practice cannot. Which practices are really
##  in the top decile, and how sure can we be?
##
##  Approach. A hierarchical (partial-pooling) binomial model with a practice-level
##  random effect and a covariate for patient mix. Small practices are shrunk
##  toward what the covariate predicts for them; large ones barely move. Ranking
##  is then done on the posterior, which answers the question directly:
##  P(practice is in the top decile), not "whose point estimate is highest".
##
##  The script compares three rankings against the simulated truth:
##    (1) naive, by observed rate;
##    (2) naive with a minimum-size rule, a common fix;
##    (3) hierarchical, by posterior probability of being top-decile.
##  Because the truth is known, the comparison is exact.
##
##  This is the "performance measure stabilization" problem in health care
##  evaluation. The model is deliberately simple; the point is the decision
##  framework and the honest accounting of uncertainty.
##
##  Reproducible: set.seed() is called once; no files are read.
##  Requires: brms (>= 2.20), rstan or cmdstanr, dplyr, tidyr, ggplot2, posterior.
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
setwd("C:/Users/goran/downloads") # Set custom, local output path for plots.
out_dir <- file.path(getwd(), "output")
dir.create(out_dir, showWarnings = FALSE)

## ---- 1. Simulate practices ---------------------------------------------------
n_practices <- 150
truth <- list(
  mu          = qlogis(0.58),  # typical practice: 58% of patients controlled
  beta_risk   = -0.45,         # sicker patient mix lowers the rate (per SD of risk)
  sd_practice =  0.50          # true between-practice spread in quality (logit)
)

practices <- tibble(
  practice   = factor(seq_len(n_practices)),
  # Heavily skewed sizes, as in real panels: many small, a few very large.
  n_patients = pmax(10L, as.integer(round(rlnorm(n_practices, log(70), 1.2)))),
  risk_score = as.numeric(scale(rnorm(n_practices))),  # standardized patient-mix risk
  u_quality  = rnorm(n_practices, 0, truth$sd_practice)  # true quality effect
) |>
  mutate(
    p_true     = plogis(truth$mu + truth$beta_risk * risk_score + u_quality),
    y          = rbinom(n(), n_patients, p_true),
    rate_obs   = y / n_patients,
    # The quantity we want to rank on: quality net of patient mix, i.e. the
    # practice's own contribution, u_quality. "Top decile" is defined on it.
    top_true   = u_quality >= quantile(u_quality, 0.9)
  )

cat(sprintf("Simulated %d practices; sizes from %d to %d patients (median %d).\n",
            n_practices, min(practices$n_patients), max(practices$n_patients),
            as.integer(median(practices$n_patients))))

## ---- 2. Two naive rankings --------------------------------------------------
## (1) Rank on the observed rate, ignoring size and patient mix.
## (2) The same, but exclude practices below 50 patients, which is how many
##     real programs handle the small-denominator problem.
k_top <- ceiling(0.10 * n_practices)
naive <- practices |>
  mutate(
    top_naive  = rank(-rate_obs) <= k_top,
    top_minsz  = n_patients >= 50 & rank(-ifelse(n_patients >= 50, rate_obs, -Inf)) <= k_top
  )

## ---- 3. Hierarchical model ---------------------------------------------------
## Binomial likelihood, logit link, practice-level random intercept, risk
## adjustment. Weakly informative priors: the practice SD prior puts most mass
## below 1 on the logit scale (odds ratios up to ~7 between a typical and an
## extreme practice), which is generous for a quality measure.
priors <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 1),   class = "b"),
  prior(exponential(2), class = "sd")
)

fit <- brm(
  y | trials(n_patients) ~ risk_score + (1 | practice),
  data = practices, family = binomial(), prior = priors,
  chains = 4, iter = 2000, warmup = 1000, seed = 3, refresh = 0, silent = 2
)

diag <- summarise_draws(as_draws_df(fit), "rhat", "ess_bulk") |>
  filter(!grepl("^r_|^lp__|^z_", variable))
stopifnot(all(diag$rhat < 1.01), all(diag$ess_bulk > 400))
cat("Convergence OK: max Rhat =", round(max(diag$rhat), 3),
    "; min bulk ESS =", round(min(diag$ess_bulk)), "\n")

## ---- 4. Parameter recovery ---------------------------------------------------
pop <- fixef(fit, probs = c(.05, .95)) |> as.data.frame() |>
  tibble::rownames_to_column("term") |>
  mutate(truth = c(truth$mu, truth$beta_risk),
         recovered = truth >= Q5 & truth <= Q95)
sdp <- VarCorr(fit)$practice$sd |> as.data.frame() |>
  mutate(term = "sd(practice)", truth = truth$sd_practice,
         recovered = truth >= Q2.5 & truth <= Q97.5)
cat("\nParameter recovery (intervals should contain the truth; about 1 miss in 10 is expected by chance):\n")
print(bind_rows(pop |> select(term, Estimate, lo = Q5, hi = Q95, truth, recovered),
                sdp |> select(term, Estimate, lo = Q2.5, hi = Q97.5, truth, recovered)) |>
        mutate(across(where(is.numeric), \(x) round(x, 2))))

## ---- 5. Posterior ranking ----------------------------------------------------
## Draws of each practice's quality effect (its random intercept). For each
## posterior draw, find which practices are in the top decile; the share of
## draws in which a practice makes it is P(top decile).
re_draws <- as_draws_matrix(fit, variable = "^r_practice\\[", regex = TRUE)
colnames(re_draws) <- levels(practices$practice)   # brms keeps level order
cut_draws <- apply(re_draws, 1, quantile, 0.9)     # top-decile cutoff per draw
in_top    <- sweep(re_draws, 1, cut_draws, ">=")

ranked <- practices |>
  mutate(
    quality_post  = colMeans(re_draws),
    quality_lo    = apply(re_draws, 2, quantile, .05),
    quality_hi    = apply(re_draws, 2, quantile, .95),
    pr_top        = colMeans(in_top),
    rank_post     = rank(-pr_top, ties.method = "first"),
    top_hier      = rank_post <= k_top
  ) |>
  left_join(naive |> select(practice, top_naive, top_minsz), by = "practice")

## ---- 6. Scorecard: which method finds the real top decile? --------------------
score <- function(flag, truth) {
  tp <- sum(flag & truth); fp <- sum(flag & !truth); fn <- sum(!flag & truth)
  tibble(selected = sum(flag), correct = tp, wrong = fp,
         precision = tp / (tp + fp), recall = tp / (tp + fn))
}
scorecard <- bind_rows(
  "Naive: observed rate"              = score(ranked$top_naive, ranked$top_true),
  "Naive with 50-patient minimum"     = score(ranked$top_minsz, ranked$top_true),
  "Hierarchical: P(top decile)"       = score(ranked$top_hier,  ranked$top_true),
  .id = "method"
)
cat("\nWho gets rewarded, and how many of them deserve it (truth known by simulation):\n")
print(scorecard |> mutate(across(c(precision, recall), \(x) round(x, 2))))

## The size bias in the naive ranking, stated plainly.
cat(sprintf(
  "\nMedian size of practices selected: naive %d, hierarchical %d, all practices %d.\n",
  as.integer(median(ranked$n_patients[ranked$top_naive])),
  as.integer(median(ranked$n_patients[ranked$top_hier])),
  as.integer(median(ranked$n_patients))
))

## ---- 7. Figures ----------------------------------------------------------------
## (a) Shrinkage: observed rate vs. the model's fitted rate for the same
##     practice (its own patient mix), point size = practice size. Large
##     practices sit on the diagonal; small ones are pulled toward what their
##     patient mix predicts.
ranked <- ranked |>
  mutate(rate_fit = fitted(fit)[, "Estimate"] / n_patients)
p_shrink <- ggplot(ranked, aes(x = rate_obs, y = rate_fit, size = n_patients)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(alpha = 0.6, colour = "steelblue4") +
  scale_size_area(max_size = 8, name = "Patients") +
  labs(x = "Observed rate", y = "Model estimate for the same practice",
       title = "Partial pooling shrinks small practices toward their expected rate",
       subtitle = "Large practices sit on the diagonal; small ones are pulled in") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "02_shrinkage.png"), p_shrink, width = 7, height = 5, dpi = 150)

## (b) Decision view: P(top decile) with the truth revealed by colour.
p_rank <- ranked |>
  arrange(pr_top) |>
  mutate(order = row_number()) |>
  ggplot(aes(x = order, y = pr_top, colour = top_true)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_point(aes(size = n_patients), alpha = 0.8) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"),
                      name = "Truly top decile") +
  scale_size_area(max_size = 6, name = "Patients") +
  labs(x = "Practices, ordered by P(top decile)", y = "Posterior P(top decile)",
       title = "Probability of being a top-decile practice",
       subtitle = "A payer can set the reward threshold on this probability, not on a point estimate") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "02_pr_top_decile.png"), p_rank, width = 8, height = 5, dpi = 150)

## (c) Caterpillar plot of the 25 highest-probability practices.
p_cat <- ranked |> slice_min(rank_post, n = 25) |>
  ggplot(aes(x = reorder(practice, pr_top), y = quality_post,
             ymin = quality_lo, ymax = quality_hi, colour = top_true)) +
  geom_pointrange() +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"),
                      name = "Truly top decile") +
  coord_flip() +
  labs(x = "Practice", y = "Quality effect (logit), 90% interval",
       title = "Top 25 practices by P(top decile)") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "02_top25_caterpillar.png"), p_cat, width = 7, height = 6, dpi = 150)

## ---- 8. Plain-language summary ------------------------------------------------
## Generated from the results, so it cannot contradict the tables above.
h <- scorecard[3, ]; nv <- scorecard[1, ]; mn <- scorecard[2, ]
sure   <- sum(ranked$pr_top >= 0.8)
unsure <- sum(ranked$pr_top > 0.2 & ranked$pr_top < 0.8)
small_wrong <- ranked |> filter(top_naive & !top_true) |> pull(n_patients)
minsize_note <- if (mn$correct <= nv$correct) {
  "A minimum-size rule does not fix this; it drops small practices that are genuinely excellent."
} else {
  "A minimum-size rule helps a little but still rewards on a noisy point estimate."
}
cat("\n--- Summary for the payer ---\n")
cat(sprintf(
  paste0("Ranking practices on their raw rate rewards %d practices, of which %d truly ",
         "belong to the top decile; the %d mistaken awards go to practices with a ",
         "median of %d patients, where a lucky year is enough to top the table. %s ",
         "The hierarchical ranking rewards %d practices, of which %d truly belong. ",
         "Beyond picking better, it says how sure we are: %d practices have at least ",
         "an 80%% chance of being top-decile, and %d sit in a grey zone (20-80%%) where ",
         "another year of data is worth more than a decision now. Recommended rule: ",
         "reward practices whose P(top decile) exceeds a threshold the payer chooses, ",
         "and report that probability with the award.\n"),
  nv$selected, nv$correct, nv$wrong, as.integer(median(small_wrong)), minsize_note,
  h$selected, h$correct, sure, unsure
))


sessionInfo()
