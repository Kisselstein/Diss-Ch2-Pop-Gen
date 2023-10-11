## Data Filtering
## December 4, 2022

#Load packages
library(tidyverse)
library(readxl)
library(lubridate)

setwd("Dissertation_Final_Draft/")

#Load data or just continue from 1B_Database script
new_dat7b <- read_rds("output/database_clean.rds")
str(new_dat7b)

#Change allele and read columns to numeric
new_dat7b$all1 <- as.numeric(new_dat7b$all1)
new_dat7b$all2 <- as.numeric(new_dat7b$all2)
new_dat7b$reads1 <- as.numeric(new_dat7b$reads1)
new_dat7b$reads2 <- as.numeric(new_dat7b$reads2)
str(new_dat7b)

# Filter data based on read count using tidyverse
# Only keep rows where there are at least 10 reads for either allele
dat10 <- new_dat7b %>% 
  filter(reads1 >= 10 | reads2 >= 10)
# Sometimes even though the first allele has enough reads, the second does not; replace those instances with NAs
dat10 <- dat10 %>% 
  mutate(reads2 = case_when(
    reads2 < 10 ~ NA_real_,
    TRUE ~ reads2
  )) %>% 
  mutate(all2 = case_when(
    is.na(reads2) ~ NA_real_,
    TRUE ~ all2
  ))

# Let's change the format of the data frame to include one column for allele and one for reads only:
# First step: split the data in two parts
dat10_pt1 <- dat10 %>% 
  mutate(all = all1,
         reads = reads1) %>% 
  select(-all1,-reads1,-all2,-reads2)
dat10_pt2 <- dat10 %>% 
  mutate(all = all2,
         reads = reads2) %>% 
  select(-all1,-reads1,-all2,-reads2)
# Bind the two parts back together
dat10 <- rbind.data.frame(dat10_pt1, dat10_pt2)

#If the filtering worked, this will return zero rows
dat10 %>% filter(reads < 10)

# Filter data based on allele frequencies
# An allele must be present in at least two samples to be kept. 
alleles_keep <- dat10 %>% 
  group_by(loci, all) %>% 
  tally() %>% 
  filter(n > 1 & !is.na(all))

# Keep only alleles that make it past the cutoff:
dat_filtered <- dat10 %>% 
  filter(loci %in% alleles_keep$loci & all %in% alleles_keep$all)

#view(dat_filtered)
head(dat_filtered)
unique(dat_filtered$site)



## Remove off-target alleles

dat_filtered %>% pull(loci) %>% unique()

## NEW METHOD - USED IN CH. 1 AS WELL 
loci_off_target <- c("EnHMG","EnMS2","SDHb_h267y","EnCDnew_24","EnCDnew_37","EnMS11",
                     "EnMS4","KHJ33757Pr9","KHJ34540PR5EX3","KHJ35605EX1pr9","EnMAT1_1_1_E",
                     "EnSLA2_F","EnHMG2","Walt_G143A") #updated

dat_filtered2a <- dat_filtered %>%
  filter(loci == "EnHMG" & all != "3"|
           loci == "EnMS2" & all != "8"|
           loci == "SDHb_h267y" & !all %in% c("5","11","9","8","7","10","15","19","14","13","16","17","26")|
           loci == "EnCDnew_24" & !all %in% c("11","1","10","13","12","18","17","16")|
           #loci == "EnCDnew_37" & all != ""|
           loci == "EnMS11" & all != "1"|
           loci == "EnMS4" & !all %in% c("10")|
           loci == "KHJ33757Pr9" & !all %in% c("4","8","24","23")|
           #loci == "KHJ34540PR5EX3" & all != ""| 
           loci == "KHJ35605EX1pr9" & !all %in% c("6","2","3","4","5","10","8","7","11")|
           #loci == "EnMAT1_1_1_E" & all != "3"| #alleles 1 and 2 both match 1-1-1 gene but very diff in length...
           loci == "EnSLA2_F" & all != "2"|
           loci == "EnHMG2" & all != "2"|
           loci == "Walt_G143A" & !all %in% c("5","15","10","17")
  )
dat_filtered2a %>% select(loci, all) %>% unique() %>% arrange(loci,all) %>% view() #test that it worked correctly

## LATER: Check EnCSEPs for off-target alleles


# All loci remaining
dat_filtered2b <- dat_filtered %>%
  filter(!loci %in% loci_off_target)

#Combine
dat_filtered2 <- rbind.data.frame(dat_filtered2a,dat_filtered2b)
nrow(dat_filtered)-nrow(dat_filtered2) #14,007 rows removed

#Save filtered data: 
# Note: read counts>10,
# allele frequency n>1 across all samples for each locus, and
# off-target alleles removed
saveRDS(dat_filtered2, "r_data/filtered_data.rds")


head(dat_filtered2)
## NOTE: Did not remove samples and loci that didn't work well overall
## LATER: Look for more off-target alleles in new loci for old output_z_g and new output_h and EnCSEPs



