
library(tidyverse)

txt <- read_csv("PMC/PMC_taxoNERD/pmc_ecology_sel4GPT.csv")
txt <- txt |> distinct_all()

chunk <- 50000L
n <- nrow(txt)
r  <- rep(1:ceiling(n/chunk),each = chunk)[1:n]
txt_spl <- split(txt, r)

for(g in seq(1, max(r))) {
  target <- txt_spl[[g]]
  instruct <- jsonlite::read_json("data/instructions_gpt35turbo2.json")
  
  model = "gpt-4o"
  temperature = 0
  sys_instruct = "system_instruction_v3"
  user_instruct = "user_ex1_v3"
  asst_instruct = "assistant_ex1_v3"
  ff <- paste0("PMC/pmc_ecology/batch_gpt4o_50K_", str_pad(g, 2, pad = "0"), "_input.jsonl")
  con <- file(ff, open = "wb")
  
  for(i in seq_len(nrow(target))) {
    jsonlite::toJSON(
      list(
        custom_id = target$ID[i],
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
  
  close(con)

}
