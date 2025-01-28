
library(tidyverse)

pmc_txt <- read_csv("PMC/PMC_taxoNERD/pmc_ecology.csv")
pmc_ner <- read_csv("PMC/PMC_taxoNERD/pmc_ecology_txn_silver.csv")

derep_sp <- function(x) {
  vec <- unique(x)
  if(length(vec) < 2) return(1L)
  d <- as.matrix(stringdist::stringdistmatrix(vec))
  d[upper.tri(d, diag = TRUE)] <- 9999L
  idx <- which(apply(d, 1, min) > 1)
  idx <- length(idx)
  idx <- ifelse(any(str_detect(vec, "[Aa]nimals?")), idx - 1, idx)
  idx <- ifelse(any(str_detect(vec, "[Pp]lants?")), idx - 1, idx)
  idx
}

n_ner <- pmc_ner %>%
  group_by(ID) %>%
  summarise(n = derep_sp(Text))

pmc_txt <- left_join(pmc_txt, n_ner, by = "ID")

pmc_txt <- filter(pmc_txt, !is.na(n))
pmc_txt <- filter(pmc_txt, n > 1)

pat_interaction <- read_lines("data/patterns_interaction.txt") |>
  paste0(collapse = "|")
match_pat_interaction <- str_detect(pmc_txt$text, pat_interaction)
pmc_txt <- mutate(pmc_txt, match_pat_interaction)

pmc_txt <- filter(pmc_txt, match_pat_interaction)

write_csv(pmc_txt, file = "PMC/PMC_taxoNERD/pmc_ecology_sel4GPT.csv")

pmc_ner <- pmc_ner |> filter(ID %in% pmc_txt$ID)
write_csv(pmc_ner, file = "PMC/PMC_taxoNERD/pmc_ecology_txn_silver_sel4GPT.csv")


pmc_txt <- read_csv(file = "PMC/PMC_taxoNERD/pmc_ecology_sel4GPT.csv")

nc <- nchar(pmc_txt$text)

sort(nc, decreasing = TRUE)[1:50]

