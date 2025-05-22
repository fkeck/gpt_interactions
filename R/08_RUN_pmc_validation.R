
library(tidyverse)

n_boot = 1000L
# TEL Validation
val_tel <- read_csv("data/validation_tel.csv") |> 
  mutate(validation = ifelse(validation == "XGen", "Ambiguous", validation))

val_tel |> 
  count(validation) |> 
  filter(!validation %in% c("OK", "SG", "SF", "WRONG")) |> 
  mutate(prop = n / sum(n) * 100) |> 
  arrange(desc(prop))

val_tel |> 
  count(validation) |> 
  filter(validation %in% c("OK", "SG", "SF", "WRONG")) |> 
  mutate(prop = n / sum(n) * 100)

1 - sum(val_tel$validation == "WRONG", na.rm = TRUE) / sum(val_tel$validation %in% c("OK", "SG", "SF", "WRONG"), na.rm = TRUE)


# LLM Validation stratified

compute_metrics <- function(tp, fp, fn) {
  res <- list()
  res$Precision <- tp / (tp + fp)
  res$Recall <- tp / (tp + fn)
  res$F1 <- 2 * ((res$Precision * res$Recall) / (res$Precision + res$Recall))
  return(unlist(res))
}

compute_metrics_multivars <- function(x, vars = c("valid_all",
                                                  "valid_interaction",
                                                  "valid_label",
                                                  "valid_category")) {
  found <- sum(!is.na(x$from))
  missed <- sum(x$missed, na.rm = TRUE)
  metrics <- vector("list", length(vars))
  names(metrics) <- vars
  
  for(i in vars) {
    wrong_interaction <- sum(x[[i]] == "FALSE", na.rm = TRUE)
    correct_interaction <- found - wrong_interaction
    
    metrics[[i]] <- compute_metrics(tp = correct_interaction, fp = wrong_interaction, fn = missed)
  }
  return(metrics |> bind_rows(.id = "variable"))
}

stratified_dat <- read_csv("data/validation_stratified.csv") |> 
  mutate(valid_interaction = valid_from & valid_to,
         valid_all = valid_interaction & valid_label & valid_category)


stratified_values <- compute_metrics_multivars(stratified_dat) |> 
  pivot_longer(Precision:F1, names_to = "metric")

stratified_dat_bs <- stratified_dat |> 
  rsample::group_bootstraps(group = ID, times = n_boot) |> 
  mutate(bs_val = map(splits, function(x) x |>
                        rsample::analysis() |>
                        compute_metrics_multivars()))

stratified_conf <- unnest(stratified_dat_bs, bs_val) |> 
  select(-splits) |> 
  pivot_longer(Precision:F1, names_to = "metric") |> 
  group_by(variable, metric) |> 
  summarise(ci_lower = quantile(value, probs = 0.025),
            ci_upper = quantile(value, probs = 0.975))

stratified_res <- left_join(stratified_values,
                            stratified_conf,
                            by = c("variable", "metric")) |> 
  mutate(metric = fct_relevel(metric, "Precision", "Recall", "F1"),
         variable = case_when(variable == "valid_all" ~ "Interaction (all)",
                              variable == "valid_interaction" ~ "Species",
                              variable == "valid_label" ~ "Label",
                              variable == "valid_category" ~ "Category") |> 
           fct_relevel("Interaction (all)", "Species", "Label", "Category"))

# Non stratified
nonstratified_dat <- read_csv("data/validation_nonstratified.csv") |> 
  mutate(valid_interaction = valid_from & valid_to,
         valid_all = valid_interaction & valid_label & valid_category)

nonstratified_values <- compute_metrics_multivars(nonstratified_dat) |> 
  pivot_longer(Precision:F1, names_to = "metric")

nonstratified_dat_bs <- nonstratified_dat |> 
  rsample::group_bootstraps(group = ID, times = n_boot) |> 
  mutate(bs_val = map(splits, function(x) x |>
                        rsample::analysis() |>
                        compute_metrics_multivars()))

nonstratified_conf <- unnest(nonstratified_dat_bs, bs_val) |> 
  select(-splits) |> 
  pivot_longer(Precision:F1, names_to = "metric") |> 
  group_by(variable, metric) |> 
  summarise(ci_lower = quantile(value, probs = 0.025),
            ci_upper = quantile(value, probs = 0.975))

nonstratified_res <- left_join(nonstratified_values,
                            nonstratified_conf,
                            by = c("variable", "metric")) |> 
  mutate(metric = fct_relevel(metric, "Precision", "Recall", "F1"),
         variable = case_when(variable == "valid_all" ~ "Interaction (all)",
                              variable == "valid_interaction" ~ "Species",
                              variable == "valid_label" ~ "Label",
                              variable == "valid_category" ~ "Category") |> 
           fct_relevel("Interaction (all)", "Species", "Label", "Category"))


# Plots
nonstratified_plot <-
  nonstratified_res |> 
  ggplot() +
  geom_col(aes(metric, value), fill = "grey60") +
  geom_linerange(aes(x = metric, ymin = ci_lower, ymax = ci_upper), alpha = 0.9) +
  facet_wrap(~variable, nrow = 1) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Original data")

stratified_plot <-
  stratified_res |> 
  ggplot() +
  geom_col(aes(metric, value), fill = "grey60") +
  geom_linerange(aes(x = metric, ymin = ci_lower, ymax = ci_upper), alpha = 0.9) +
  facet_wrap(~variable, nrow = 1) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Stratified data")

library(patchwork)
nonstratified_plot / stratified_plot + plot_annotation(tag_levels = "a")

