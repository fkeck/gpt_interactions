library(tidyverse)
source("R/FUN_taxonomic_entity_linking.R")

dat <- list.files("PMC/pmc_ecology/batch_openai/output/", full.names = TRUE) |> 
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

for (i in 1:length(tax_links$new)) {
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

write_csv(g_data, "data/save_R_gdata_2.csv")

