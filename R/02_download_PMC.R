
###### Match PMC OA and openalex ecology  and download from PMC #####

library(tidyverse)

openalex <- read_csv("PMC/openalex_ecology_IDs.csv")
pmc <- read_csv("https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv")

openalex_PMC <- openalex |>
  filter(!is.na(pmcid)) |> 
  select(pmcid) |> 
  deframe() |> 
  str_remove("https://www.ncbi.nlm.nih.gov/pmc/articles/")
  
openalex_PMC <- str_c("PMC", openalex_PMC)

pmc_filtered <- pmc |> filter(`Accession ID` %in% openalex_PMC)

pmc_files <- str_c("https://ftp.ncbi.nlm.nih.gov/pub/pmc/", pmc_filtered$File)
dest_files <- str_c("PMC/PMC_ecology/", basename(pmc_filtered$File))


for (i in 1:length(pmc_files)) {
  curl::curl_download(pmc_files[i], destfile = dest_files[i])
  cat(i, "\n")
}



###### Extract xml files from PMC archives #####

ff <- list.files("PMC/PMC_ecology", pattern = ".tar.gz$", full.names = TRUE)

setwd("PMC/extracted/")

for (i in ff) {
  cmd <- paste0("tar -xf ", i, " --wildcards --no-anchored '*.nxml'")
  system(command = cmd)
  cat(i, "\n")
}

list.files("PMC/extracted/", pattern = ".nxml$", recursive = TRUE) |> length()
file.remove(ff)
