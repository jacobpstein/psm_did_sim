###############################
# This code implements a simulation showing 
# the consequences of conditioning on a treatment
# affected outcome in a propensity score matching model
# Session Info:
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.5
################################

# load packages
library(tidyverse) # for data cleaning
library(MatchIt) # matching
library(broom) # making nice output

# set seed
set.seed(42)

# create a function to implement and iterate over our simulation
simulate_once <- function(n = 10000) {
  # Simulate baseline data
  id_data <- tibble(
    id = 1:n
    , x = rnorm(n)
    , trend_slope = runif(n, -0.3, 0.3)
  )|>
    mutate(
      treat = rbinom(n, 1, plogis(0.5 * x + 1.5 * trend_slope))
    )
  
  # Create long data for 4 time points
  long_data <- id_data|>
    crossing(time = 0:3)|>
    mutate(
       post = ifelse(time == 3, 1, 0)
      , time_trend = trend_slope * time
      , y = 2 + 0.5 * x + time_trend + 2 * treat * post + rnorm(n * 4)
    )
  
  #  Create post-treatment covariate at time 2
  # This simulates a post-treatment outcome-derived covariate
  y_time2 <- long_data|>
    filter(time == 2) |>
    select(id, y_time2 = y)
  
  long_data <- long_data|>
    left_join(y_time2, by = "id")
  
  
  # Big no no match on x and post-treatment y_time2
  psm_bad <- matchit(treat ~ x + y_time2, data = pre_treatment, method = "nearest")
  matched_ids_bad <- match.data(psm_bad)$id
  
  did_psm_bad <- long_data |>
    filter(id %in% matched_ids_bad) |>
    lm(y ~ treat * post, data = .) |>
    tidy() |>
    filter(term == "treat:post") |>
    mutate(method = "Improper PSM (post-treatment) + DiD")
  
  #  Implement DR DiD using inverse probability weights
  long_data <- long_data|>
    mutate(
      ps = plogis(0.5 * x + 1.5 * trend_slope)
      , weight = ifelse(treat == 1, 1 / ps, 1 / (1 - ps))
    )
  
  dr_did <- lm(y ~ treat * post, data = long_data, weights = weight) |>
    tidy(). |>
    filter(term == "treat:post"). |>
    mutate(method = "DR DiD")
  
  # Combine all estimates
  bind_rows(did_psm_bad, dr_did)
}

# Run it!
results <- map_dfr(1:10000, ~simulate_once())

# Summarize results
summary_stats <- results |>
  group_by(method) |>
  summarise(
    mean_estimate = mean(estimate)
    , se = sd(estimate)
    , .groups = "drop"
  )


# Plot results
p1 <- results|>
  filter(method != "Proper PSM + DiD") |> 
  ggplot(aes(x = estimate, fill = method)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  labs(title = "10,000 runs of PSM + DiD with the treatment affected outcome included\nin the matching model compared to DR DiD without the bias of treatment affected variables"
       , subtitle = "Dotted line shows the true effect of 2.0"
       , x = "Treatment Effect Estimate", y = "") +
  geom_text(data = data.frame(x = c(1.97182622474081, 2.11703245041384 ),
                              y = c(1.84100819808392, 1.84100819808392),
                              label = c("Doubly Robust", "PSM + DID")),
            mapping = aes(x = x, y = y, label = label),
            colour = "white", family = "Gill Sans", inherit.aes = FALSE, size = 6) +
  usaid_plot() +
  theme(axis.text.y = element_blank()
        )


ggsave('did_comps.png', p1, height = 8, width = 12, dpi = 600)
