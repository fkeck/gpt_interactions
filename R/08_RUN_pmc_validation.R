
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


cat_missed <- nonstratified_dat |> 
  group_by(corrected_label) |> 
  summarise(n_missed = sum(missed, na.rm = TRUE)) |> 
  rename(category = corrected_label)

cat_error <- nonstratified_dat |> 
  filter(!valid_all) |>
  group_by(corrected_label) |> 
  count(name = "n_error") |> 
  rename(category = corrected_label)


cat_test <- nonstratified_dat |> 
  filter(valid_all, !is.na(from)) |> 
  group_by(category) |> 
  count() |> 
  left_join(cat_missed) |> 
  left_join(cat_error) |> 
  replace_na(replace = list(n_missed = 0, n_error = 0))

cat_test |> 
  mutate(category = str_to_title(category)) |> 
  setNames(c("Category", "Correct", "Missed", "Incorrect"))

rbind(cond1 = cat_test$n, cond2 = cat_test$n_missed) |> 
  chisq.test(simulate.p.value = TRUE, B = 200000)

rbind(cond1 = cat_test$n, cond2 = cat_test$n_error) |> 
  chisq.test(simulate.p.value = TRUE, B = 200000)

tab <- rbind(cond1 = cat_test$n, cond2 = cat_test$n_missed)
colnames(tab) <- cat_test$category




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


##### Assess LLM consistency #####

# Make batch
nonstratified_dat <- read_csv("data/validation_nonstratified.csv") |> 
  mutate(valid_interaction = valid_from & valid_to,
         valid_all = valid_interaction & valid_label & valid_category)

target <- nonstratified_dat |> 
  filter(ID %in% sample(unique(nonstratified_dat$ID), 100, replace = FALSE)) |> 
  distinct(ID, text)
instruct <- jsonlite::read_json("data/instructions_gpt35turbo2.json")

model = "gpt-4o-2024-05-13"
temperature = 0
sys_instruct = "system_instruction_v3"
user_instruct = "user_ex1_v3"
asst_instruct = "assistant_ex1_v3"
ff <- paste0("PMC/pmc_ecology/batch_consistency_input.jsonl")
con <- file(ff, open = "wb")

for(s in 1:100) {
  for(i in seq_len(nrow(target))) {
    jsonlite::toJSON(
      list(
        custom_id = paste0(target$ID[i], "_S", sprintf("%03d", s)),
        method = "POST",
        url = "/v1/chat/completions",
        body = list(
          model = model,
          temperature = temperature,
          messages = list(
            list(
              role = "system",
              type = "text",
              content = instruct[[sys_instruct]]
            ),
            list(
              role = "user",
              type = "text",
              content = instruct[[user_instruct]]
            ),
            list(
              role = "assistant",
              type = "text",
              content = instruct[[asst_instruct]]
            ),
            list(
              role = "user",
              type = "text",
              content = target$text[i]
            )
          )
        )
      ), auto_unbox = TRUE) |> 
      write_lines(con, append = TRUE)
    cat('\r', i)
    flush.console()
  }
}

close(con)

# Upload and processing via the web interface
# TEL
source("R/FUN_taxonomic_entity_linking.R")

dat <- list("PMC/pmc_ecology/batch_consistency_output.jsonl") |> 
  map(function(x) {
    con <- file(x, open = "r")
    dat <- jsonlite::stream_in(con)
    close(con)
    dat <- as_tibble(dat)
    return(dat)
  })

dat <- bind_rows(dat)

dat <- dat$response$body$choices |>
  bind_rows() |>
  as_tibble() |>
  mutate(ID = dat$custom_id, .before = 1)

dat$finish_reason |> table()
dat <- filter(dat, finish_reason == "stop")


dat$message$content[!str_detect(dat$message$content, "^from\\|to\\|label\\|category")] <- 
  dat$message$content[!str_detect(dat$message$content, "^from\\|to\\|label\\|category")] |>
  str_remove("^.*\\n")

dat$message$content <- str_remove(dat$message$content, "\\n```$")

dat$message$content[dat$message$content == "```"] <- ""

dat <- filter(dat, message$content != "")

res <- dat$message$content %>%
  map(function(ct) {
    ct <- ifelse(ct == "" | !str_detect(ct, "^from\\|to\\|label\\|category"), "from|to|label|category", ct)
    tryCatch({
      read.delim(text = ct, sep = "|")
    }, error = function(e) {
      read.delim(text = "from|to|label|category", sep = "|")
    })%>%
      as_tibble() %>% 
      mutate(across(where(function(x) !is.character(x)), as.character))
  })

res <- res |> 
  bind_rows(.id = "ID") |> 
  mutate(ID = dat$ID[as.numeric(ID)])

res <- res |> filter(!(from == "hawks" & to == "lemmings"),
                     !(from == "bats" & to == "Carica papaya"),
                     !(from == "capra ibex" & to == "festuca"))


tax <- unique(c(res$from, res$to))


tax_links <- tibble(old = tax,
                    new = tax |>
                      clean_special_chr() |> 
                      clean_vernacular_names() |>
                      clean_scientific_name() |> 
                      clean_phage_names()) |> 
  mutate(old == new) |>
  filter(nchar(new) > 2L) |> 
  distinct(new, .keep_all = TRUE)

write_rds(dat, "PMC/pmc_ecology/save_R_dat.rds")
write_rds(res, "PMC/pmc_ecology/save_R_res.rds")
write_rds(tax_links, "PMC/pmc_ecology/save_R_tax_links.rds")

res <- read_rds("PMC/pmc_ecology/save_R_res.rds")
tax_links <- read_rds("PMC/pmc_ecology/save_R_tax_links.rds")


##### Perform TEL #####
itis <- read_csv("data/vernacular_itis.csv")
eol <- read_csv("data/vernacular_eol.csv")
mix_itis_eol <- inner_join(itis, eol) |>
  distinct()

uids <- vector("list", length = length(tax_links$new))

for (i in 51:length(tax_links$new)) {
  uids[[i]] <- resolve_taxo_retry(tax_links$new[i],
                                  vernacular_itis_ref = mix_itis_eol,
                                  vernacular_eol_ref = mix_itis_eol)
  cli::cli_text(i)
  
  if(i %% 500L == 0L) {
    write_rds(uids, "PMC/pmc_ecology/save_R_uids_temp.rds")
  }
}

# Second pass
for (i in which(map_int(uids, ~attributes(.x) |> length()) == 0)) {
  uids[[i]] <- resolve_taxo_retry(tax_links$new[i],
                                  vernacular_itis_ref = itis,
                                  vernacular_eol_ref = eol)
  cli::cli_text(i)
}


write_rds(uids, "PMC/pmc_ecology/save_R_uids_2.rds")
uids <- read_rds("PMC/pmc_ecology/save_R_uids_2.rds")

uids <-
  map(uids, function(x) {
    tibble(uid = as.character(x),
           match = attr(x, "match"),
           multiple_matches = attr(x, "multiple_matches"),
           match_confidence = attr(x,"taxonomic_match_confidence"))}) |>
  list_rbind()

tax_links <- tax_links |> 
  mutate(uid = uids$uid)

uids <- uids |> 
  filter(!is.na(uid)) |>
  distinct(uid, .keep_all = TRUE)



##### Collect NCBI classification #####
uids_split <- split(uids, seq(nrow(uids)) %/% 100L)
uids_split_classif_raw <- vector("list", length = length(uids_split))

for(i in seq_along(uids_split)) {
  time_out = 30L
  delay_retry = 10L
  n_retry = 3L
  res <- "error"
  while (identical(res, "error") & n_retry > 0L) {
    res <- tryCatch({
      R.utils::withTimeout({
        as.uid(uids_split[[i]]$uid, check = FALSE) |>
          classification(db = "ncbi")
      }, timeout = time_out)
    },
    TimeoutException = function(cond) {
      message("Timeout. Skipping.")
      return("error")
    },
    error = function(cond) {
      message("\nSomething went wrong:")
      message(cond)
      message("\n")
      for (i in delay_retry:0) {
        cat("\rRetrying in", i, "s.  ")
        Sys.sleep(1)
      }
      cat("\n")
      return("error")
    })
    n_retry <- n_retry - 1L
  }
  
  if(identical(res, "error")) {
    message(i, " --> ERROR")
    uids_split_classif_raw[[i]] <- NULL
  } else {
    message(i, " --> OK")
    uids_split_classif_raw[[i]] <- res
  }
}

write_rds(uids_split_classif_raw, "PMC/pmc_ecology/save_R_uids_split_classif_raw.rds")
uids_split_classif_raw <- read_rds("PMC/pmc_ecology/save_R_uids_split_classif_raw.rds")


uids <- uids_split_classif_raw |>
  map(unclass) |>
  map(function(x) x[map_lgl(x, is.data.frame)]) |> 
  map(function(x) map(x, function(x) {
    lowest <- tibble(name = x$name[nrow(x)],
                     rank = "lowest_rank",
                     id = x$id[nrow(x)])
    x <- bind_rows(x, lowest)
    return(x)
  })) |>
  map(bind_rows, .id = "uid") |>
  bind_rows() |> 
  arrange(uid) |> 
  as_tibble()


##### Link TEL ids to classification #####
tax_links <- tax_links |>
  left_join(uids, by = "uid")

tax_links_long <- tax_links |>
  select(old, uid, name, rank, id) |> 
  dplyr::filter(rank %in% c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "lowest_rank"))

g_data <- res |> 
  left_join(tax_links_long |>
              rename(from_uid = id, from_name = name), by = c("from" = "old"), relationship = "many-to-many") |> 
  left_join(tax_links_long |>
              rename(to_uid = id, to_name = name), by = c("to" = "old", "rank" = "rank")) |> 
  filter(!is.na(from_uid), !is.na(to_uid))

write_csv(g_data, "PMC/pmc_ecology/save_rep_gdata.csv")
g_data <- read_csv("PMC/pmc_ecology/save_rep_gdata.csv")

g_data |> 
  mutate(ID_para = str_sub(ID, 1, -6), .after = ID) |> 
  filter(rank == "lowest_rank")

g_data |> 
  mutate(ID_para = str_sub(ID, 1, -6), .after = ID) |> 
  filter(rank == "species") |> 
  distinct(ID, from_name, to_name, label, category, .keep_all = TRUE) |> 
  mutate(comb = paste(from_name, to_name, label, category, sep = "|")) |> 
  group_by(ID_para) |> 
  count(comb) |> 
  ungroup() |> 
  mutate(p_hat = n/100) |> 
  mutate(pw_agreement = p_hat^2 + (1 - p_hat)^2) |> 
  summarise(mean(pw_agreement))


g_data |> 
  mutate(ID_para = str_sub(ID, 1, -6), .after = ID) |> 
  filter(rank == "lowest_rank") |> 
  distinct(ID, from_name, to_name, category, .keep_all = TRUE) |> 
  mutate(comb = paste(from_name, to_name, category, sep = "|")) |> 
  group_by(ID_para) |> 
  count(comb) |> 
  ungroup() |> 
  mutate(p_hat = n/100) |> 
  ggplot() +
  geom_histogram(aes(p_hat), bins = 10)
