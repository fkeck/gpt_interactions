nerviz <- function(texts, ents,
                   html_file = tempfile(fileext = "html"),
                   browse = TRUE,
                   ent_style = "none") {
  
  html <- vector("list", nrow(texts))
  html_titles <- vector("list", nrow(texts))
  html_p <- vector("list", nrow(texts))
  
  ent_col <- 
    c("LIVB" = "#f89880")
  
  for(i in seq_len(nrow(texts))) {
    ents_i <- filter(ents, ID == texts$ID[i])

    # Text
    sst <- str_split(texts$text[i], "")[[1]]
    
    min0 <- function(x) {if (length(x)>0) min(x) else Inf}
    
    if(min0(ents_i$Start_char) == 0) {
      splits_1 <- c(ents_i$Start_char +1 , ents_i$End_char + 1) %>% unique() %>% sort()
      splits_2 <- c(ents_i$Start_char[-1], ents_i$End_char, length(sst)) %>% unique() %>% sort()
    } else {
      splits_1 <- c(0, ents_i$Start_char +1 , ents_i$End_char + 1) %>% unique() %>% sort()
      splits_2 <- c(ents_i$Start_char, ents_i$End_char, length(sst)) %>% unique() %>% sort()
    }
    
    if(length(splits_1) > length(splits_2)) splits_1 <- splits_1[-length(splits_1)]
    sst <- map2_chr(splits_1, splits_2, ~ paste0(sst[.x:.y], collapse = ""))
    
    for(it in seq_len(nrow(ents_i))) {
      idx <- which(c(0, splits_2) == ents_i$Start_char[it])
      if(ent_style == "normal") {
        sst[idx] <- paste(sst[idx], ents_i$Label[it])
      }
      if(ent_style == "superscript") {
        sst[idx] <- paste(sst[idx], "<sup>", ents_i$Label[it], "</sup>")
      }
      sst[idx] <- paste0('<mark title= ', ents_i$Label[it], ' style="background: ', ent_col[ents_i$Label[it]], ';padding: 0.3em 0.25em; margin: 0 0.1em; line-height: 1; border-radius: 0.35em;">', sst[idx], "</mark>")
    }
    
    
    html_titles[[i]] <- paste0(URLdecode(texts$file[i]) %>% str_remove("\\.[a-z]+$"),
                              "<br>",
                              texts$section[i],
                              ifelse(is.na(texts$subsection[i]), "", paste0(" / ", texts$subsection[i]))) %>% 
       paste0('<h1 style="font-size: 1em;line-height: 2; direction: ltr">', ., "</h1>")
    html_p[[i]] <- paste0(sst, collapse = "") %>% 
      paste0('<p style="font-size: 0.9em;line-height: 1.8; direction: ltr; text-align: justify; text-justify: inter-word;">', ., "</p1>")
    
    html[[i]] <- paste('<div style="background: #f1f2f4; margin-bottom: 4rem; margin-left: 1.5rem; margin-right: 4rem; border-radius: 12px;
box-shadow: var(--ds-shadow-raised, 0px 1px 1px #091e4240, 0px 0px 1px #091e424f); padding: 12px; padding-top: 2px; padding-bottom: 2px">',
                       html_titles[[i]],
                       html_p[[i]],
                       "</div>") %>% 
      htmltools::HTML()
  }
  
  html %>% htmltools::save_html(file = html_file)
  
  if(browse) {
    browseURL(html_file)
  }
  
  cat(html_file)
  
}



doi_to_url <- function(x) {
  paste0("http://dx.doi.org/", x)
}

url_to_doi <- function(x) {
  str_remove(x, "^https://doi\\.org/")
}

doi_to_file <- function(x) {
  URLencode(x, reserved = TRUE)
}

file_to_doi <- function(x) {
  URLdecode(x)
}


insert_sp_tag <- function(texts, ents) {

    html_p <- vector("character", nrow(texts))
  
  for(i in seq_len(nrow(texts))) {
    ents_i <- filter(ents, ID == texts$ID[i])
    
    # Text
    sst <- str_split(texts$text[i], "")[[1]]
    
    min0 <- function(x) {if (length(x)>0) min(x) else Inf}
    
    if(min0(ents_i$Start_char) == 0) {
      splits_1 <- c(ents_i$Start_char +1 , ents_i$End_char + 1) %>% unique() %>% sort()
      splits_2 <- c(ents_i$Start_char[-1], ents_i$End_char, length(sst)) %>% unique() %>% sort()
    } else {
      splits_1 <- c(0, ents_i$Start_char +1 , ents_i$End_char + 1) %>% unique() %>% sort()
      splits_2 <- c(ents_i$Start_char, ents_i$End_char, length(sst)) %>% unique() %>% sort()
    }
    
    if(length(splits_1) > length(splits_2)) splits_1 <- splits_1[-length(splits_1)]
    sst <- map2_chr(splits_1, splits_2, ~ paste0(sst[.x:.y], collapse = ""))
    
    for(it in seq_len(nrow(ents_i))) {
      idx <- which(c(0, splits_2) == ents_i$Start_char[it])
      sst[idx] <- paste0('<sp>', sst[idx], "</sp>")
    }
    
    html_p[[i]] <- paste0(sst, collapse = "")
    
  }
    return(html_p)
}



spacy_json_to_csv <- function(json, csv) {
  paste("jq -r '
  [[\"ID\",\"Text\",\"Label\",\"Start_char\",\"End_char\"]], 
   [ .[] 
     | [ .ID ] + 
       (.ents[] | [.text,.label,.start_char,.end_char]) ] 
  | .[] 
  | @csv'", json, ">", csv) %>% 
    system()
}


