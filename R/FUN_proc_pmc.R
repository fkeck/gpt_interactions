
library(tidyverse)
library(xml2)


#' Process PMC XML files
#' 
#' Parser for PMC XML format. Returned dataframe is not directly
#' compatible with other editor process_* functions.
#'
#' @param x a PMC XML document (xml2 xml_document object)
#'
#' @return a tibble with paragraphs classified into sections
#' 
process_pmc <- function(x, remove_references = TRUE,
                        expand_genus_names = TRUE,
                        remove_formulas = TRUE,
                        remove_tables = TRUE,
                        remove_figures = TRUE) {
  
  # Extract PMCID
  pmc_id <-
    x |>
    xml_find_first("//article-id[@pub-id-type='pmc']") |>
    xml_text()
  
  # Extract Main Text
  p_nodes <-
    x |>
    xml_find_all("//body//p")
  
  parents <-
    p_nodes |>
    map(xml_parents)
  
  sec_parent <-
    parents |>
    map_lgl(~.x[[1]] |>
              xml_name() |>
              magrittr::is_in(c("sec", "body")))
  
  p_nodes <- p_nodes[sec_parent]
  parents <- parents[sec_parent]
  
  # Extract sections
  sections <-
    parents |>
    map_df(function(x) {
      tag <- xml_name(x)
      
      if(tag[1] == "body") {
        section_title <- NA_character_
        section_type <- NA_character_
        subsection_title <- NA_character_
        subsection_type <- NA_character_
      } else {
        n_tag <- length(tag)
        n_sec <- sum(tag == "sec")
        is_sec <- which(tag == "sec")
        
        section <- x[[max(is_sec)]]
        section_title <- xml_find_first(section, "title") |>
          xml_text() |> 
          str_remove("^[0-9\\.]* *") |> 
          str_remove("\\.$") |> 
          str_replace_all(" & ", " and ") |> 
          str_to_title()
        section_type <- xml_attr(section, "sec-type")
        
        if(n_sec > 1L) {
          subsection <- x[[sort(is_sec, TRUE)[2]]]
          subsection_title <- xml_find_first(subsection, "title") |>
            xml_text() |> 
            str_remove("^[0-9\\.]* *") |> 
            str_remove("\\.$") |> 
            str_replace_all(" & ", " and ") |> 
            str_to_title()
          subsection_type <- xml_attr(subsection, "sec-type")
        } else {
          subsection_title <- NA_character_
          subsection_type <- NA_character_
        }
        
      }
      res <- tibble(section_title, section_type, subsection_title, subsection_type)
      
      return(res)
    })
  
  if(remove_formulas) {
    formulas <- xml_find_all(x, "//inline-formula|//disp-formula-group|//disp-formula")
    xml_remove(formulas)
  }
  
  if(remove_tables) {
    tables <- xml_find_all(x, "//table-wrap-group|//table-wrap|//table-wrap-foot|//table")
    xml_remove(tables)
  }
  
  if(remove_figures) {
    figures <- xml_find_all(x, "//fig-group|//fig")
    xml_remove(figures)
  }
  
  
  #' Remove references from paragraphs
  #'
  #' @param x a paragraph xml node.
  #' @param mode method to identify and remove references (xref or regex).
  #'
  #' @return 
  #' A cleaned character string.
  remove_ref <- function(x, mode) {
    
    if(mode == "xref") {
      sups <- xml_find_all(x, "//sup[//xref]")
      xml_remove(sups)
      
      xrefs <- xml_find_all(x, "//xref")
      xml_remove(xrefs)
      
      txt <- xml_text(x)
      txt <- str_remove_all(txt, " *\\[[;,\\-– ]*\\] *")
      txt <- str_remove_all(txt, " *\\([;,\\-– ]*\\) *")
    }
    
    if(mode == "regex") {
      txt <- xml_text(x) %>% 
        str_remove_all("[:space:]\\([^\\)]*(1[7-9][0-9]{2}|20[0-2][0-9])[^\\)]*\\)") %>%
        str_remove_all("[:space:][\\[\\(][^[:alpha:]]*[\\]\\)]")
    }
    return(txt)
  }
  
  if(remove_references) {
    p_text <- p_nodes |> map_chr(remove_ref, mode = "regex")
    
  } else {
    p_text <- p_nodes |> xml_text()
  }
  
  if(expand_genus_names) {
    
    fulltext <- paste0(p_text, collapse = " ")
    
    gens <- p_nodes %>%
      xml_find_all('//italic') %>%
      xml_text() %>%
      unique() %>%
      enframe(name = NULL, value = "genus") %>%
      filter(str_detect(genus, "^[A-Z][a-z]*$")) %>% 
      filter(!str_detect(genus, "ae$")) %>% 
      filter(nchar(genus) > 1) %>% 
      deframe()
    
    gens_sps <- p_nodes %>%
      xml_find_all('//italic') %>%
      xml_text() %>%
      str_split(",[:blank:]|;[:blank:]") %>% 
      unlist() |> 
      unique() %>%
      enframe(name = NULL, value = "genus") %>%
      filter(str_detect(genus, "^[A-Z][a-z]+[:blank:][a-z]{2,}$")) %>% 
      deframe()
    
    p_text <- map_chr(p_text, function(x){
      sp_short <- str_extract_all(x, "[A-Z][a-z]?\\.[:blank:][a-z]{2,}") %>%
        unlist() %>%
        unique()
      for(i in sp_short) {
        sp_short_gen <- str_split(i, "\\.[:blank:]")[[1]][1]
        sp_short_epithet <- str_split(i, "\\.[:blank:]")[[1]][2]
        pat_long <- paste0(sp_short_gen, "[a-z]+[:blank:]", sp_short_epithet)
        sp_long <- str_extract_all(fulltext, pat_long)[[1]] %>%
          unique()
        
        if(length(sp_long) != 1) {
          sp_long <- gens_sps[str_detect(gens_sps, pat_long)]
        }

        if(length(sp_long) != 1) {
          gen_long <- gens[str_detect(gens, paste0("^", sp_short_gen))]
          if(length(gen_long) < 1) { next }
          if(length(gen_long) > 1) {
            gen_long <- gen_long[str_detect(fulltext, paste0("\\bgen[a-z\\.]*\\b[:blank:]", gen_long))]
          }
          if(length(gen_long) != 1) { next }
          sp_long <- paste(gen_long, sp_short_epithet)
        }
        x <- str_replace_all(x, i, sp_long)
      }
      return(x)
    })
  }
  
  res <- mutate(sections, text = p_text) |> 
    mutate(ID = paste0("PMC", pmc_id, "_P", str_pad(row_number(), 3, "left", 0)), .before = 1L)
  
  return(res)
}
