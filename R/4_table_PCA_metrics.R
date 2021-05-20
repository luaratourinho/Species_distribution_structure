library(dplyr)
library(tidyr)
library(purrr)
library(readr)


# Organizing the table ----------------------------------------------------


arquivos <- list.files("./Results_2_landscapemetrics", full.names = TRUE)

read_2 <- function(x) {
  sp.name <- gsub("./Results_2_landscapemetrics/", "", x)
  sp.name <- gsub(".csv", "", sp.name)
  
  read_csv2(x) %>%
    select(function_name, value) %>%
    spread(function_name, value) %>%
    mutate(Species = sp.name)
}

data <- arquivos %>%
    purrr::map(read_2) %>%  
    reduce(rbind) %>% 
    separate(Species, sep = "_", 
             c("Genus", "Species", "Data", "Time", "Projection"))


data2 <- data %>% 
  filter(Time == "cu") %>% #select time
  mutate(Species_complete = paste(Genus, Species)) # Juntar nomes
  


# PCA ---------------------------------------------------------------------


