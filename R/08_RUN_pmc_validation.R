
library(tidyverse)


dat <- read_csv("data/validation_data.csv")

n_all <- nrow(dat)
n_para <- dat$ID |> unique() |> length()

inter_found <- sum(!is.na(dat$from))
inter_missed <- sum(dat$n_missed, na.rm = TRUE)
inter_wrong <- sum(dat$valid == "FALSE") + sum(!is.na(dat$corrected_label))
inter_correct <- inter_found - inter_wrong

inter_real <- inter_correct + inter_missed

true_pos <- inter_correct
false_pos <- inter_wrong
false_neg <- inter_missed


(precision <- true_pos / (true_pos + false_pos))
(recall <- true_pos / (true_pos + false_neg))


