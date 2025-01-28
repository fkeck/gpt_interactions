library(tidyverse)

source("R/FUN_proc_pmc.R")

files <- list.files("PMC/extracted/",
                    full.names = TRUE, recursive = TRUE)


con_out <- file("PMC/PMC_taxoNERD/pmc_ecology.jsonl", open = "wb")

for(i in seq_along(files)) {
  xml <- read_xml(files[i])
  res <- process_pmc(xml)
  jsonlite::stream_out(res, con_out, pagesize = 1000L)
  print(i)
}

close(con_out)

####### RUN TAXONERD FROM PYTHON SCRIPT ########

source("R/FUN_nerviz.R")

spacy_json_to_csv(json = "PMC/PMC_taxoNERD/pmc_ecology_silver_ann.json",
                  csv = "PMC/PMC_taxoNERD/pmc_ecology_silver_ann.csv")

px <- read_csv("PMC/PMC_taxoNERD/pmc_ecology_silver_ann.csv")

px <- px |> select(-Label)

write_csv(px, "PMC/PMC_taxoNERD/pmc_ecology_txn_silver.csv")

qx <- jsonlite::stream_in(file("PMC/PMC_taxoNERD/pmc_ecology.jsonl", open = "rb"))
qx <- as_tibble(qx)
write_csv(qx, "PMC/PMC_taxoNERD/pmc_ecology.csv")

