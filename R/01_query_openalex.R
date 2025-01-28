
# Use the openalex API to download all the relevant references

library(tidyverse)
library(openalexR)
library(httr2)

query_url <-
  oa_query(entity = "works",
           options = list(select = c("type", "title", "publication_year", "ids")),
           topics.subfield.id = "subfields/2303") |>
  paste0("&mailto=XXXX",
         "&per-page=200")

query_url <- paste0(query_url, "&cursor=")
res_path <- "PMC/openalex/"
cursor <- "*"
iter <- 1L
block <- 1L
dat <- vector("list", 50L)
cols <- c(title = NA_character_,
          publication_year = NA_integer_,
          openalex = NA_character_,
          doi = NA_character_,
          mag = NA_character_,
          pmid = NA_character_,
          pmcid = NA_character_,
          type = NA_character_)

while(!is.null(cursor)) {

  cat("Block:", block, "........ Iteration ", iter, "/50\n", sep = "")
  
  query_url_cursor <- paste0(query_url, cursor)
  req <- request(query_url_cursor)
  resp <- req_perform(req)
  
  if(resp$status_code != 200) {
    stop("Server error")
  }
  
  json <- resp_body_json(resp)
  
  if(length(json$results) == 0) {
    dat <- dplyr::bind_rows(dat)
    readr::write_csv(dat, file.path(res_path, sprintf("%08d.csv", block)))
    break
  }
  
  dat[[iter]] <- json$results |>
    works2df() %>% 
    tidyr::unnest_wider(ids) %>% 
    add_column(!!!cols[setdiff(names(cols), names(.))]) %>% 
    dplyr::select(title, publication_year, type,
                  openalex, doi, mag, pmid, pmcid)
  
  if(iter == 50L) {
    dat <- dplyr::bind_rows(dat)
    readr::write_csv(dat, file.path(res_path, sprintf("%08d.csv", block)))
    iter <- 0L
    block <- block + 1L
    dat <- vector("list", 50L)
  }
  
  cursor = json$meta$next_cursor
  iter = iter + 1L
}


# Merging files
dat <- list.files(res_path, full.names = TRUE) %>% 
  map(read_csv, col_types = "cicccccc")

dat <- bind_rows(dat)

write_csv(dat, "PMC/openalex_ecology_IDs.csv")
