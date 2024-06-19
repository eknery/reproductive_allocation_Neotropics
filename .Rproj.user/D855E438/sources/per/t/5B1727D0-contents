### libraries
if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("Hmisc")) install.packages("Hmisc"); library("Hmisc")
if (!require("PupillometryR")) install.packages("PupillometryR"); library("PupillometryR")
if (!require("plyr")) install.packages("plyr"); library("plyr")
if (!require("ape")) install.packages("ape"); library("ape")

### read data
habitat_mtx = read.table("0_data/habitat_matrix.csv", 
                         h=T, sep=",", na.strings = "na")

################################ CLASSIFICATION ##############################

### count habtiats
sort(apply(habitat_mtx[,-1], FUN = sum, MARGIN = 2))

### habitat names
habitat_names = habitat_mtx %>% 
  select(! species) %>% 
  colnames()

forest = c(
  "Floresta_Ombrofila", 
  "Floresta_Ombrofila_Mista",
  "Floresta_de_Terra_Firme",
  "Floresta_Ciliar"
  )

open = c( 
   "Afloramento_rochoso",
   "Campinarana", 
   "Carrasco",
   "Cerrado",
   "Campo_de_Varzea",
   "Campo_Rupestre", 
   "Restinga",
   "Savana_Amazonica"
   )

other = habitat_names[!habitat_names %in% c(forest, open)]

### presence per habitat
habitat_pres = habitat_mtx %>% 
  pivot_longer(cols = any_of(habitat_names), 
               names_to = "habitat",
               values_to = "presence")

### assigning habitat types
habitat_type = habitat_pres %>% 
  mutate(type = case_when(
    (habitat %in% forest)     & presence == 1  ~ "forest",
    (habitat %in% open)       & presence == 1  ~ "open",
    (habitat %in% other)      & presence == 1  ~ "other"
  )
  ) %>% 
  filter( !is.na(type) )

### counting habitat types
habitat_count = habitat_type %>% 
  pivot_wider(names_from = type, 
              values_from = presence,
              values_fill = 0) %>% 
  group_by(species) %>% 
  reframe(
    n_forest = sum(forest), 
    n_open = sum(open),
    n_other = sum(other)
  )

### defining habtiat range
habitat_range = habitat_count %>% 
  mutate(
    range = case_when(
      n_forest >= 1  & n_open == 0  & n_other == 0  ~ "forest_specialist",
      n_forest == 0  & n_open >= 1  & n_other == 0  ~ "open_specialist",
      TRUE ~ "generalist"
    )
  ) %>% 
  mutate(range = factor(range, levels = c("forest_specialist",
                                          "generalist",
                                          "open_specialist"
  )
  )
  )

table(habitat_range$range)

### with state probibility
habitat_range_prob = habitat_range %>% 
  mutate(prob = 1) %>% 
  pivot_wider(names_from = range,
              names_expand = T,
              values_from = prob,
              values_fill = 0
  )

### exporting 
saveRDS(habitat_range, "1_habitat_results/habitat_range.RDS")
