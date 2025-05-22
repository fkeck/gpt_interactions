
library(tidyverse)
library(taxize)
library(pluralize)
library(cli)


clean_scientific_name <- function(x) {
  
  x <- stringr::str_remove_all(x, # Words containing numbers if preceded by another word
                                 "\\w*[0-9]+\\w*\\s*")
  
  x <- stringr::str_remove_all(x, # Words between parentheses
                               "\\(([^\\)]+)\\)")
  
  x <- stringr::str_remove_all(x, # Words after opening parenthesis
                               "\\(.*")
  x <- stringr::str_remove_all(x, # Remove remaining dots
                               "\\.{2,}$| \\.+||^\\.+")
  x <- stringr::str_remove_all(x, # Remove words after sp. but not if the first is ending with a dot
                               "(?<= sp\\.) [a-z]+($| ).*")
  x <- stringr::str_remove_all(x, # Remove word after sp. nov.
                               "(?<=[A-Z][a-z]{2,30} [a-z]{2,30} sp\\. nov\\.).*")
  
  x <- stringr::str_remove_all(x, # Remove all punctuation at the end of the string
                               "[:punct:]+$")
  x <- stringr::str_remove_all(x, # Remove single character left at the end of the strings
                               "[:space:].$")

  .REGEX_UNCERTAIN <- "(?i)\\b(aff|cf|prox|nr|sp\\.? inc|indet|ca)\\b\\.?"
  x <- str_replace_all(x, pattern = .REGEX_UNCERTAIN, replacement = " ")
  
  x <- str_replace_all(x, pattern = "(?i)\\b(fam|gen|sp|spp)\\b[^[:alpha:]]*", replacement = "")
  x <- str_remove_all(x, pattern = "(?i)-like\\b")
  
  
  x <- stringr::str_squish(x)
  x <- stringr::str_trim(x)
  
  return(x)
}

clean_vernacular_names <- function(x) {
  x <- stringr::str_remove_all(x, #List of qualifiers
                               "(?i)\\b(groups?|clades?|of|species|types?|strains?|organisms?|females?|males?|offsprings?|larva[le]?|juveniles?|adults?|chicks?|natives?|seeds?|leaves|leaf|flowers?|aliens?|non.?natives?|exotic|introduced|invasives?|indigenous|residents?|omnivorous|fungivorous|bacterivorous|bactivorous|herbivorous|carnivorous|insectivorous|piscivorous|benthivorous|granivorous|planktivorous|parasitic|parasites?|(ecto|endo)parasites?|hosts?|other|pelagic|benthic|planktonic|sessile|terrestrial|freshwater|marine|deep-sea|soil|populations?|assemblages?|community|communities)\\b")
  x <- stringr::str_remove_all(x, "\\b(OTUs?|MOTUs?)\\b")
    x <- stringr::str_squish(x)
  x <- stringr::str_trim(x)
  return(x)
}

clean_phage_names <- function(x) {
  phage_match <- stringr::str_detect(x, "(?i)(\\bphages?|(bacterio|coli|cyano|pelagi|vibrio|anky|podo).?phages?)\\b")
  x[phage_match] <- "unclassified bacterial viruses"
  return(x)
}

clean_special_chr <- function(x) {
  x <- str_replace_all(x, pattern = "\\p{Pd}|−", replacement = "-")
  x <- str_replace_all(x, pattern = "'|’|‘|ʼ|ʹ|ʻ|ʼ", replacement = "'")
  x <- str_remove_all(x, pattern = "“|”")
  return(x)
}


wikipedia_classif <- function(x) {
  
  wk_p <- wikitaxa::wt_wiki_url_build(wiki = "en", type = "wikipedia", page = x, utf8 = TRUE) |>
    wikitaxa::wt_wiki_page()
  
  if(!is.null(wk_p$response_headers_all[[1]]$`mediawiki-api-error`)) {
    return(NA)
  }
  
  wk_html <- wk_p$content |>
    rawToChar() |> 
    xml2::read_html()
  
  wk_table <- wk_html |> 
    xml2::xml_find_first("//table[contains(@class, 'infobox')]")
  
  xml2::xml_find_all(wk_table, '//small') |> 
    xml2::xml_remove()
  
  
  if(is.na(wk_table)) {
    return(NA)
  }
  
  wk_table <- wk_table |> 
    rvest::html_table(header = FALSE, convert = FALSE)
  
  if(ncol(wk_table) < 2L) {
    return(NA)
  }
  
  if(!is.character(wk_table$X1) | !is.character(wk_table$X2)) {
    return(NA)
  }
  
  wk_table <- wk_table |> 
    mutate(across(where(is.character), function(x) {
      str_replace_all(x, "\\\\n", " ") |> 
        str_remove_all(":") |> 
        str_remove_all("(?<=[a-z])[A-Z].*") |> 
        str_trim() |> 
        str_squish()
    })) |> 
    filter(str_detect(X1, "^Domain|^Kingdom|^Phylum|^Class|^Order|^Family|^Genus|^Species")) |> 
    filter(!is.na(X2))
  
  if(nrow(wk_table) == 0L) {
    return(NA)
  }
  
  wk_table <- wk_table |> 
    mutate(X2 = str_remove_all(X2, "†"))
  
  if(str_starts(wk_table$X1[nrow(wk_table)], "Species")) {
    if(str_detect(wk_table$X2[nrow(wk_table)], "[A-Z]\\. ")) {
      gen <- wk_table$X2[wk_table$X1 == "Genus"]
      wk_table$X2[nrow(wk_table)] <- str_replace(wk_table$X2[nrow(wk_table)], "[A-Z]\\.", gen)
    }
  }
  return(wk_table)
}


#' Resolve names to NCBI taxonomy UID
#'
#' @param x the name of the taxon (scientific or vernacular).
#' @param verbose logical. If TRUE prints information to the console.
#' @param use_itis_links logical. If TRUE uses Wikipedia links to ITIS. Require ITIS API.
#' @param vernacular_itis_ref a dataframe with two columns \code{vernacular} and \code {scientific}.
#' @param vernacular_eol_ref a dataframe with two columns \code{vernacular} and \code {scientific}.
#'
#' @return an NCBI taxonomic identifier as returned by the \code{get_uid} function from \pkg{taxize}.
#' 
resolve_taxo <- function(x, verbose = TRUE,
                         use_itis_links = FALSE,
                         use_wikipedia_box = TRUE,
                         vernacular_itis_ref = NULL,
                         vernacular_eol_ref = NULL) {
  
  if(verbose) {
    cli_div(theme = list(span.emph = list(color = "blue")))
    cli_h3("Processing {.emph {x}}.")
    cli_end()
  }
    
  
  conf <- NA

  match_eol <- c(vernacular_eol_ref$scientific[which(vernacular_eol_ref$vernacular == pluralize::singularize(str_to_lower(x)))],
                  vernacular_eol_ref$scientific[which(vernacular_eol_ref$vernacular == pluralize::pluralize(str_to_lower(x)))])
  
  match_itis <- c(vernacular_itis_ref$scientific[which(vernacular_itis_ref$vernacular == pluralize::singularize(str_to_lower(x)))],
                  vernacular_itis_ref$scientific[which(vernacular_itis_ref$vernacular == pluralize::pluralize(str_to_lower(x)))])
  
  match_vern <- c(match_eol, match_itis) |>
    table() |>
    sort(decreasing = TRUE) |>
    names()

  if(length(match_vern) > 0) {
    if(length(match_vern) == 1) {
      x <- match_vern
      if(verbose) cli_alert_info("Matched in itis/eol vernacular names: translated to {x}.")
    } else {
      match_vern_uid <- get_uid(match_vern, ask = FALSE, messages = FALSE)
      match_vern <- match_vern[!is.na(match_vern_uid)][1]
      if(!is.na(match_vern)) {
        x <- match_vern
        if(verbose) cli_alert_info("Matched in itis/eol vernacular names: translated to {x}.")
      }
    }
  }
  if(!is.na(taxize::gna_verifier(x, data_sources = 1:12)$dataSourceId[1])) {
    gnr <- taxize::gna_verifier(x, 1:12, all_matches = TRUE) |>
      janitor::clean_names() |>
      distinct(data_source_id, .keep_all = TRUE)
  } else {
    gnr <- tibble()
  }

  if(nrow(gnr) > 0) {
    if(sum(gnr$submitted_name == gnr$matched_name) > 3) {
      gnr_sci <- TRUE
      if(verbose) cli_alert_info("GNR (>3 matches).")
    } else {
      if("NCBI" %in% gnr$data_source_title_short) {
        gnr <- gnr |>
          dplyr::filter(data_source_title_short == "NCBI") |> 
          dplyr::mutate(exact = submitted_name == matched_name)
        
        if(max(gnr$sort_score) > 9.4) {
          x <- gnr$matched_name[which.max(gnr$sort_score)]
          gnr_sci <- TRUE
          if(verbose) cli_alert_info("GNR (Score NCBI >0.98).")
        } else {
          if(any(gnr$exact)) {
            x <- gnr$matched_name[gnr$exact][1]
            gnr_sci <- TRUE
            if(verbose) cli_alert_info("GNR (NCBI exact match).")
          } else {
            gnr_sci <- FALSE
          }
        }
      } else {
        gnr_sci <- FALSE
      }
    }
  } else {
    gnr_sci <- FALSE
  }
  
  # Scientific name
  if(gnr_sci) {
    uid <- get_uid(x, ask = FALSE, messages = FALSE)
    
    if(is.na(uid)) {
      uid_syn <- get_uid(x, ask = FALSE, messages = FALSE, modifier = "Synonym")
      if(!is.na(uid_syn)) {
        uid <- uid_syn
        if(verbose) cli_alert_success("Solved with NCBI scientific name (synonym).")
        conf <- 20
      }
    } else {
      if(verbose) cli_alert_success("Solved with NCBI scientific name.")
      conf <- 10
    }
    
    if(is.na(uid)) {
      if(attr(uid, "match") == "NA due to not found") {
        if(verbose) cli_alert_danger("Scientific name not found.")
      }
      if(attr(uid, "match") == "NA due to ask=FALSE & > 1 result") {
        
        uid_df <- get_uid_(x, messages = FALSE)[[1]]
        if(length(unique(uid_df$division)) == 1L) {
          rank_order <- str_split(taxize::rank_ref$ranks, ",") |> unlist()
          highest_ranked_row <- which.min(match(uid_df$rank, rank_order))
          lowest_ranked_row <- which.max(match(uid_df$rank, rank_order))
          lowest_ranked_classif <- get_uid(x, rows = lowest_ranked_row, ask = FALSE, messages = FALSE) |> classification()
          if(all(uid_df$scientificname %in% lowest_ranked_classif[[1]]$name)) {
            uid <- get_uid(x, rows = highest_ranked_row, ask = FALSE, messages = FALSE)
          }
        }
        if(!is.na(uid)) {
          conf <- 30
          if(verbose) cli_alert_success("Solved after resolving name conflicts.")
        } else {
          if(verbose) cli_alert_danger("Scientific name conflicts")
        }
      }
    }
    
    if(is.na(uid)) {
      wk_links <- wikitaxa::wt_wikipedia(x, wiki = "en")$externallinks
      wk_links_ncbi <- wk_links[str_detect(wk_links, "www.ncbi.nlm.nih.gov/Taxonomy/")]
      if(length(wk_links_ncbi) > 0L) {
        uid_vec <- str_extract(wk_links_ncbi, "(?<=&id=)[0-9]+")
        if(length(uid_vec) == 1L) {
          uid <- uid_vec
        }
        if(!is.na(uid)) {
          conf <- 40
          if(verbose) cli_alert_success("Solved scientific name using Wikipedia links to NCBI.")
        }
      }
    }
    
    if(use_itis_links) {
      if(is.na(uid)) {
        wk_links <- wikitaxa::wt_wikipedia(x, wiki = "en")$externallinks
        wk_links_itis <- wk_links[str_detect(wk_links, "https://www.itis.gov/")]
        if(length(wk_links_itis) > 0L) {
          itis_tsn <- str_extract(wk_links_itis, "(?<=&search_value=)[0-9]+")
          if(length(itis_tsn) == 1L) {
            itis_sci_name <- itis_hierarchy(itis_tsn, "up")$taxonname
            uid <- get_uid(itis_sci_name, ask = FALSE, messages = FALSE)
          }
          if(!is.na(uid)) {
            conf <- 75
            if(verbose) cli_alert_success("Solved scientific name using Wikipedia links to ITIS.")
          }
        }
      }
    }
    
    # Common name
  } else {
    uid <- get_uid(sci_com = x, modifier = "Common Name", messages = FALSE, ask = FALSE)
    
    if(!is.na(uid)) {
      conf <- 50
      if(verbose) cli_alert_success("Solved with direct NCBI vernacular name match.")
    }
    
    # Deal with multiple matches
    if(is.na(uid) & attr(uid, "multiple_matches")) {
      all_matches <- get_uid_(sci_com = x, modifier = "Common Name", messages = verbose)[[1]]
      string_dist <- stringdist::stringdist(all_matches$commonname, x)
      best_row <- which(string_dist == min(string_dist))
      if(length(best_row) == 1L) {
        uid <- get_uid(sci_com = x, modifier = "Common Name", messages = FALSE, ask = FALSE, rows = best_row)
      }
      if(!is.na(uid)) {
        conf <- 55
        if(verbose) cli_alert_success("Solved with fuzzy NCBI vernacular name match.")
      }
    }
    
    # Try plural/singular
    if(is.na(uid)) {
      if(pluralize::is_plural(x)) {
        uid <- get_uid(sci_com = pluralize::singularize(x), modifier = "Common Name", messages = FALSE, ask = FALSE)
      } else {
        uid <- get_uid(sci_com = pluralize::pluralize(x), modifier = "Common Name", messages = FALSE, ask = FALSE)
      }
      if(!is.na(uid)) {
        conf <- 60
        if(verbose) cli_alert_success("Solved with fuzzy NCBI vernacular name match.")
      }
    }
    
    # Try solving through Wikipedia page links to NCBI
    if(is.na(uid)) {
      wk_links <- wikitaxa::wt_wikipedia(pluralize::singularize(x), wiki = "en")$externallinks
      wk_links_ncbi <- wk_links[str_detect(wk_links, "www.ncbi.nlm.nih.gov/Taxonomy/")]
      if(length(wk_links_ncbi) > 0L) {
        uid_vec <- str_extract(wk_links_ncbi, "(?<=&id=)[0-9]+")
        if(length(uid_vec) == 1L) {
          uid <- uid_vec
        }
        if(!is.na(uid)) {
          conf <- 70
          if(verbose) cli_alert_success("Solved using Wikipedia links to NCBI.")
        }
      }
    }
    if(is.na(uid)) {
      wk_links <- wikitaxa::wt_wikipedia(pluralize::pluralize(x), wiki = "en")$externallinks
      wk_links_ncbi <- wk_links[str_detect(wk_links, "www.ncbi.nlm.nih.gov/Taxonomy/")]
      if(length(wk_links_ncbi) > 0L) {
        uid_vec <- str_extract(wk_links_ncbi, "(?<=&id=)[0-9]+")
        if(length(uid_vec) == 1L) {
          uid <- uid_vec
        }
        if(!is.na(uid)) {
          conf <- 70
          if(verbose) cli_alert_success("Solved using Wikipedia links to NCBI.")
        }
      }
    }
    
    # Try solving through Wikipedia page links to ITIS
    if(use_itis_links) {
      if(is.na(uid)) {
        wk_links <- wikitaxa::wt_wikipedia(pluralize::singularize(x), wiki = "en")$externallinks
        wk_links_itis <- wk_links[str_detect(wk_links, "https://www.itis.gov/")]
        if(length(wk_links_itis) > 0L) {
          itis_tsn <- str_extract(wk_links_itis, "(?<=&search_value=)[0-9]+")
          if(length(itis_tsn) == 1L) {
            itis_sci_name <- itis_hierarchy(itis_tsn, "up")$taxonname
            uid <- get_uid(itis_sci_name, ask = FALSE, messages = FALSE)
          }
          if(!is.na(uid)) {
            conf <- 75
            if(verbose) cli_alert_success("Solved using Wikipedia links to ITIS.")
          }
        }
      }
      if(is.na(uid)) {
        wk_links <- wikitaxa::wt_wikipedia(pluralize::pluralize(x), wiki = "en")$externallinks
        wk_links_itis <- wk_links[str_detect(wk_links, "https://www.itis.gov/")]
        if(length(wk_links_itis) > 0L) {
          itis_tsn <- str_extract(wk_links_itis, "(?<=&search_value=)[0-9]+")
          if(length(itis_tsn) == 1L) {
            itis_sci_name <- itis_hierarchy(itis_tsn, "up")$taxonname
            uid <- get_uid(itis_sci_name, ask = FALSE, messages = FALSE)
          }
          if(!is.na(uid)) {
            conf <- 75
            if(verbose) cli_alert_success("Solved using Wikipedia links to ITIS.")
          }
        }
      }
    }
  }
  
  if(is.na(uid)) {
    if(str_detect(x, "(virus|viruses)$")) {
      x <- "Viruses"
      uid <- get_uid(x, messages = FALSE, ask = FALSE)
      if(!is.na(uid)) {
        conf <- 90
        if(verbose) cli_alert_success("Solved with exception 'Viruses'.")
      }
    }
  }
  
  # Recurrent passes
  if(is.na(uid)) {
    if(str_detect(x, "\\b[A-Z][a-z]+ [a-z]{2,}") && (str_count(x, "[:space:]") > 1L)) {
      gensp <- str_extract_all(x, "\\b[A-Z][a-z]+ [a-z]{2,}")[[1]]
      if(length(gensp) == 1L) {
        uid <- resolve_taxo(gensp, use_wikipedia_box = use_wikipedia_box)
        if(!is.na(uid)) {
          conf <- 100
          if(verbose) cli_alert_success("Solved by extracting genus and species.")
        }
      }
    }
  }
  
  if(is.na(uid)) {
    if(str_detect(x, "\\b[A-Z][a-z]+\\b") && (str_count(x, "[:space:]") > 0L)) {
      gen <- str_extract_all(x, "\\b[A-Z][a-z]+\\b")[[1]]
      if(length(gen) == 1L) {
        uid <- resolve_taxo(gen, use_wikipedia_box = use_wikipedia_box)
        if(!is.na(uid)) {
          conf <- 110
          if(verbose) cli_alert_success("Solved by extracting genus.")
        }
      }
    }
  }
  
  if(is.na(uid)) {
    if(str_detect(x, "^[a-z-]+ ")) {
      shrank <- str_remove(x, "^[a-z-]+ ")[[1]]
      uid <- resolve_taxo(shrank,
                          use_wikipedia_box = use_wikipedia_box,
                          vernacular_itis_ref = vernacular_itis_ref,
                          vernacular_eol_ref = vernacular_eol_ref)
      if(!is.na(uid)) {
        conf <- 120
        if(verbose) cli_alert_success("Solved by shrining left words.")
      }
    }
  }
  
  if(is.na(uid)) {
    if(use_wikipedia_box) {
      x_wiki <- wikipedia_classif(x)
      if(!is.na(x_wiki)) {
        x_wiki <- x_wiki$X2[nrow(x_wiki)]
        if(verbose) cli_alert_info("Matched in Wikipedia box: translated to {x_wiki}.")
        uid <- resolve_taxo(x_wiki, use_wikipedia_box = FALSE)
        if(!is.na(uid)) {
          conf <- 130
        }
      }
    }
  }

  if(is.na(uid)) {
    conf <- -999
    if(verbose) cli_alert_danger("Fail to solve.")
  }

  attr(uid, "taxonomic_match_confidence") <- conf
  return(uid)
}

resolve_taxo_retry <- function(x, ..., time_out = 15, delay_retry = 5, n_retry = 3L) {
  
  res <- "error"
  
  while (identical(res, "error") & n_retry > 0L) {
    
    res <- tryCatch({
      R.utils::withTimeout({
        resolve_taxo(x, ...)
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
    return(NA)
  } else {
    return(res)
  }
}




