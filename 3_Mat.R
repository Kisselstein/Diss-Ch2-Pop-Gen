## Determine Sample Mating Type & Clonality 
## after cleaning, making database, and filtering
## before phasing or plotting allele frequencies

#Load packages
library(tidyverse)
library(readxl)
library(lubridate)

#Load filtered dataset
# Note: read counts>10,
# allele frequency n>1 across all samples for each locus, and
# off-target alleles removed
## NOTE: Did not remove samples and loci that didn't work well overall
## LATER: Look for more off-target alleles in new loci for output_z_g
dat_filtered2 <- read_rds("r_data/filtered_data.rds")

unique(dat_filtered2$year)

# Only want to look at my own samples and data
dat_filtered3 <- dat_filtered2
#dat_filtered3 <- dat_filtered2 %>% 
#  filter(year %in% c('2015','2016','2017','2018','2019','2020','Control'))


## Determine Sample Mating Type
unique(dat_filtered2$loci)
unique(dat_filtered3$loci)
mat1_1_loci <- c("EnMAT1_1_3_A", "EnMAT1_1_1_B", "EnMAT1_1_1_D", 
                 "EnMAT1_1_3_C", "EnMAT1_1_1_A", "EnMAT1_1_3_B",
                 "HPM175mat1_11", "PodaphMAT1_1")
# LATER: What happened to last 2 and Enalpha????

mat1_2_loci <- c("EnHMG", "EnHMG2", "PodaphMAT1_2", "HPM200MAT1_21")
# LATER: What happened to the last 2?


dat1_1 <- dat_filtered3 %>% 
  filter(loci %in% mat1_1_loci) %>% 
  mutate(mat1 = "1") %>% 
  select(sample_id,mat1) %>% 
  distinct()
head(dat1_1)  

dat1_2 <- dat_filtered3 %>% 
  filter(loci %in% mat1_2_loci) %>% 
  mutate(mat2 = "2") %>% 
  select(sample_id,mat2) %>% 
  distinct()
head(dat1_2)

dat_filtered4 <- full_join(dat1_1,dat1_2)
head(dat_filtered4)
nrow(dat_filtered4) == dat_filtered3 %>% filter(loci %in% mat1_1_loci | loci %in% mat1_2_loci) %>% select(-loci,-all,-reads) %>% distinct() %>% nrow()

dat_filtered5 <- dat_filtered3 %>% 
  full_join(dat_filtered4) %>% 
  mutate(mat = case_when(
    mat1 == "1" & mat2 == "2" ~ "3",
    mat1 == "1" & is.na(mat2) ~ "1",
    is.na(mat1) & mat2 == "2" ~ "2",
    is.na(mat1) & is.na(mat2) ~ "0",
    TRUE ~ "oops"
  )) %>% 
  select(-mat1,-mat2) #keep r1_raw for now
head(dat_filtered5)
nrow(dat_filtered2) == nrow(dat_filtered5) #should have added 3 columns and that's it
unique(dat_filtered5$mat)


write_rds(dat_filtered5,"r_data/filtered_data_mat.rds")
## LATER: Determine mat using a 1/20 ratio instead? (Normalize reads first?? Naw.)


