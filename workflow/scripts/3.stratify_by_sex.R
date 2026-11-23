#stratify_by_sex.R
#This script stratifies the count matix by sex and diagnosis

#Libraries --- --- 

pacman::p_load("tidyverse",
               "vroom")

#args --- ---
args <- commandArgs(trailingOnly = TRUE)
counts <- args[1]
metadata <- args[2]

# Read RDS counts
counts <- readRDS(counts)
counts <- as.data.frame(counts)

# Read metadata
metadata <- vroom(metadata)

# Stratify metadata --- ---

#AD women and men

metadata_msex0_cerad_1 <- metadata %>% filter(sex == "0") %>% filter(ceradsc == "1")
metadata_msex1_cerad_1 <- metadata %>% filter(sex == "1")  %>% filter(ceradsc == "1")

#Healthy women and men 

metadata_msex0_cerad_4 <- metadata %>% filter(sex == "0") %>% filter(ceradsc == "4")
metadata_msex1_cerad_4 <- metadata %>% filter(sex == "1")  %>% filter(ceradsc == "4")

# Filter counts --- ---

counts_msex0_cerad_1 <- counts %>% dplyr::select(all_of(metadata_msex0_cerad_1$specimenID))
counts_msex1_cerad_1 <- counts %>% dplyr::select(all_of(metadata_msex1_cerad_1$specimenID))

counts_msex0_cerad_4 <- counts %>% dplyr::select(all_of(metadata_msex0_cerad_4$specimenID))
counts_msex1_cerad_4 <- counts %>% dplyr::select(all_of(metadata_msex1_cerad_4$specimenID))

#Save data --- ---

arg_counts_msex0_cerad_1 <- args[3]
arg_counts_msex1_cerad_1 <-  args[4]

arg_counts_msex0_cerad_4 <-  args[5]
arg_counts_msex1_cerad_4 <-  args[6]

saveRDS(counts_msex0_cerad_1, arg_counts_msex0_cerad_1)
saveRDS(counts_msex1_cerad_1, arg_counts_msex1_cerad_1)
saveRDS(counts_msex0_cerad_4, arg_counts_msex0_cerad_4)
saveRDS(counts_msex1_cerad_4, arg_counts_msex1_cerad_4)

#Save metadata --- ----

arg_metadata_msex0_cerad_1 <- args[7]
arg_metadata_msex1_cerad_1 <-  args[8]

arg_metadata_msex0_cerad_4 <-  args[9]
arg_metadata_msex1_cerad_4 <-  args[10]

vroom::vroom_write(metadata_msex0_cerad_1, arg_metadata_msex0_cerad_1)
vroom::vroom_write(metadata_msex1_cerad_1, arg_metadata_msex1_cerad_1)
vroom::vroom_write(metadata_msex0_cerad_4, arg_metadata_msex0_cerad_4)
vroom::vroom_write(metadata_msex1_cerad_4, arg_metadata_msex1_cerad_4)

#END