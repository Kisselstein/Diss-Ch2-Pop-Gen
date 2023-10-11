## Add meta data

#Load packages
library(tidyverse)
library(readxl)
library(lubridate)

#setwd("..")

# Save cleaned data
new_dat4 <- readRDS("Dissertation_Final_Draft/r_data/cleaned_data.rds")

## Add meta data for samples: date collected, disease inc and sev, etc
# Add date, disease incidence, and disease severity data to our ampseq data

# Import disease data for 2019
disease19 <- read_excel("../../Data/Tutoring/Ch2/input/2019_Disease_summary_fixed2.xlsx", sheet = 1) %>%
  mutate(year = 2019) %>% 
  mutate(year = as.character(year)) %>% 
  mutate(date = as.character(date))
head(disease19)

disease2019 <- disease19 %>%
  mutate(avg_inc = as.character(avg_inc)) %>% 
  mutate(avg_sev = as.character(avg_sev)) %>% 
  mutate(variety = case_when(
    site == "AR" ~ "CH",
    site == "BH5" ~ "CN5",
    site == "BHA" ~ "AU",
    site == "FRC" ~ "CH",
    site == "LLA" ~ "AU",
    site == "LLC" ~ "CH"
  )) %>%
  mutate(site = case_when(
    str_detect(site, "CR[0-9]*") ~ str_replace(site, "CR","CR_"),
    TRUE ~ site
  )) %>%
  mutate(trt = case_when(
    str_detect(site, "CR_[0-9]*") ~ word(site, 2, 2, "_"),
    TRUE ~ NA_character_
  )) %>%
  mutate(site = case_when(
    str_detect(site, "CR_[0-9]*") ~ "CR",
    site == "AR" ~ "AR",
    site == "BH5" ~ "BH",
    site == "BHA" ~ "BH",
    site == "FRC" ~ "FR",
    site == "LLA" ~ "LL",
    site == "LLC" ~ "LL"
  )) %>%
  mutate(variety = case_when(
    site == "CR" ~ "CH",
    TRUE ~ variety
  )) %>%
  mutate(timepoint = case_when(
    site == "CR" & timepoint == "1" ~ "T1",
    site == "CR" & timepoint == "2" ~ "T2",
    site == "CR" & timepoint == "3" ~ "T3",
    site == "CR" & timepoint == "4" ~ "T4",
    TRUE ~ timepoint
  )) %>% 
  mutate(avg_inc = case_when(
    avg_inc == "0" ~ NA_character_,
    TRUE ~ avg_inc
  )) %>% 
  mutate(avg_sev = case_when(
    avg_sev == "0" ~ NA_character_,
    TRUE ~ avg_sev
  ))
head(disease2019)
unique(disease2019$site)
unique(disease2019$variety)
unique(disease2019$trt)
unique(disease2019$timepoint)


# Import disease data for 2018
disease18 <- read_excel("../../Data/Tutoring/Combined/input/Samples_Date_Disease_Summary.xlsx", sheet = "2018") %>%
  select(site = Site, date = Date, avg_inc = `Avg Incidence (%)`, avg_sev = `Avg Severity (%)`) %>% 
  mutate(year = 2018) %>% 
  mutate(year = as.character(year))
head(disease18)

# Split timepoint info in new column
disease18 <- disease18 %>% 
  
  mutate(timepoint = word(string = site, start = 2, end = 2, sep = "_")) %>% 
  mutate(site = word(string = site, start = 1, end = 1, sep = "_"))
unique(disease18$site)

disease2018 <- disease18 %>%
  mutate(variety = case_when(
    site == "FR" ~ "CH",
    site == "LL" ~ "CH",
    site == "AR" ~ "CH",
    site == "BH" ~ "AU",
    site == "BHD" ~ "DE",
    site == "BH5" ~ "CN5",
  )) %>%
  mutate(site = case_when(
    str_detect(site, "CR[0-9]*") ~ str_replace(site, "CR","CR_"),
    site == "BHD" ~ "BH",
    TRUE ~ site
  )) %>%
  mutate(trt = case_when(
    str_detect(site, "CR_[0-9]*") ~ word(site, 2, 2, "_"),
    TRUE ~ NA_character_
  )) %>%
  mutate(site = case_when(
    str_detect(site, "CR[0-9]*") ~ "CR",
    site == "BH5" ~ "BH",
    TRUE ~ site
  )) %>% 
  mutate(timepoint = case_when(
    site == "LLO" ~ NA_character_,
    TRUE ~ timepoint
  )) %>% 
  mutate(variety = case_when(
    is.na(variety) & site == "CR" | site == "AR" | site == "FR" | site == "LL" ~ "CH",
    TRUE ~ variety
  )) %>% 
  mutate(avg_inc = case_when(
    avg_inc == "NA" ~ NA_character_,
    TRUE ~ avg_inc
  )) %>% 
  mutate(avg_sev = case_when(
    avg_sev == "NA" ~ NA_character_,
    TRUE ~ avg_sev
  )) %>% 
  mutate(date = case_when(
    date == "NA" ~ NA_character_,
    TRUE ~ date
  )) 
unique(disease2018$year)
unique(disease2018$site)
unique(disease2018$variety)
disease2018 %>% filter(is.na(variety))
unique(disease2018$trt)
unique(disease2018$timepoint)
disease2018 %>% filter(is.na(timepoint))
unique(disease2018$date)


# Import disease data for 2017

# Make this data match our ampseq data
# CR trts:
# S-10 Luna Experience (fluopyram and low tebuconazole) 
# Q-13 (azoxystrobin and flutriafol)
# D-25 (tetraconazole)
# E-12 (flutriafol)
# Control-33

disease17 <- read_excel("../../Data/Tutoring/Combined/input/Samples_Date_Disease_Summary.xlsx", sheet = "2017") %>%
  select(site = Site_New, variety = Variety_New, date = Date, avg_inc = `Avg Incidence (%)`, avg_sev = `Avg Severity (%)`) %>% 
  mutate(year = "2017")
head(disease17)


disease2017 <- disease17 %>%
  mutate(trt = case_when(
    site == "CR" ~ variety,
    TRUE ~ NA_character_
  )) %>% 
  mutate(variety = case_when(
    site == "CR" ~ "CH",
    variety == "NA" ~ NA_character_,
    TRUE ~ variety
  )) %>% 
  mutate(avg_inc = case_when(
    avg_inc == "NA" ~ NA_character_,
    TRUE ~ avg_inc
  )) %>% 
  mutate(avg_sev = case_when(
    avg_sev == "NA" ~ NA_character_,
    TRUE ~ avg_sev
  )) %>% 
  mutate(date = case_when(
    date == "NA" ~ NA_character_,
    TRUE ~ date
  )) 
head(disease2017,10)
unique(disease2017$site)
unique(disease2017$variety)
unique(disease2017$trt)


# Import disease data for 2016
new_dat4 %>% filter(year=="2016") %>% select(site,variety,trt) %>% distinct()

disease16 <- read_excel("../../Data/Tutoring/Combined/input/Samples_Date_Disease_Summary.xlsx", sheet = "2016") %>%
  select(site = Site_New, variety = Variety_New, date = Date, avg_inc = `Avg Incidence (%)`, avg_sev = `Avg Severity (%)`) %>% 
  mutate(year = "2016")
head(disease16)


# Make this data match our ampseq data
disease2016 <- disease16 %>%
  mutate(date = as.character(date)) %>% 
  mutate(trt = case_when(
    site == "CR" ~ variety,
    TRUE ~ NA_character_
  )) %>%  
  mutate(variety = case_when(
    site == "CR" ~ "CH",
    TRUE ~ variety
  )) %>% 
  mutate(avg_inc = case_when(
    avg_inc == "NA" ~ NA_character_,
    TRUE ~ avg_inc
  )) %>% 
  mutate(avg_sev = case_when(
    avg_sev == "NA" ~ NA_character_,
    TRUE ~ avg_sev
  )) %>% 
  mutate(date = case_when(
    date == "NA" ~ NA_character_,
    TRUE ~ date
  )) 
head(disease2016)
unique(disease2016$site)
unique(disease2016$variety)
unique(disease2016$trt)


# Import disease data for 2015
disease15 <- read_excel("../../Data/Tutoring/Combined/input/Samples_Date_Disease_Summary.xlsx", sheet = "2015") %>%
  dplyr::select(site = Site_New, variety = Variety_New, date = Date, avg_inc = `Avg Incidence (%)`, avg_sev = `Avg Severity (%)`) %>% 
  mutate(year = "2015") %>% 
  mutate(date = as.character(date))
head(disease15)

# Make this data match our ampseq data
disease2015 <- disease15 %>%
  mutate(trt = case_when(
    site == "CR" ~ variety,
    TRUE ~ NA_character_
  )) %>%  
  mutate(variety = case_when(
    site == "CR" ~ "CH",
    TRUE ~ variety
  )) %>% 
  mutate(avg_inc = case_when(
    avg_inc == "NA" ~ NA_character_,
    TRUE ~ avg_inc
  )) %>% 
  mutate(avg_sev = case_when(
    avg_sev == "NA" ~ NA_character_,
    TRUE ~ avg_sev
  )) %>% 
  mutate(date = case_when(
    date == "NA" ~ NA_character_,
    TRUE ~ date
  )) %>% 
  filter(!is.na(site))
head(disease2015)
unique(disease2015$site)
unique(disease2015$variety)
unique(disease2015$trt)



# Add date and disease data to the ampseq data
disease2019 <- disease2019 %>% mutate(date = ymd(date)) %>% 
  select(year,site,variety,trt,timepoint,avg_inc,avg_sev,date)
disease2018 <- disease2018 %>% mutate(date = ymd(date)) %>% 
  select(year,site,variety,trt,timepoint,avg_inc,avg_sev,date)
disease2017 <- disease2017 %>% mutate(date = ymd(date)) %>% 
  mutate(timepoint = NA_character_) %>% 
  select(year,site,variety,trt,timepoint,avg_inc,avg_sev,date)
disease2016 <- disease2016 %>% mutate(date = ymd(date)) %>% 
  mutate(timepoint = NA_character_) %>% 
  select(year,site,variety,trt,timepoint,avg_inc,avg_sev,date)
disease2015 <- disease2015 %>% mutate(date = ymd(date)) %>% 
  mutate(timepoint = NA_character_) %>% 
  select(year,site,variety,trt,timepoint,avg_inc,avg_sev,date)

fix_na_var <- data.frame (year = c("2018", "2018", "2018", "2018", "2018", "2018", "2018"),
                          site = c("W", "W", "W", "W", "LLO", "LLO", "LLO"),
                          variety = c("VR", "VL", "VR", "VL", "AU", "CH", NA_character_),
                          timepoint = c("T1", "T1", "T2", "T2", NA_character_, NA_character_, NA_character_),
                          date = c("20180905", "20180905", "20181012","20181012", "20180809", "20180809", "20180809")) %>% 
  mutate(trt=NA_character_) %>% 
  mutate(avg_inc = case_when(
    site == "LLO" ~ "40",
    TRUE ~ NA_character_
  )) %>% 
  mutate(avg_sev = case_when(
    site == "LLO" ~ "20",
    TRUE ~ NA_character_
  )) %>% 
  mutate(date = ymd(date)) %>% 
  select(year,site,variety,trt,timepoint,avg_inc,avg_sev,date)
fix_na_var
str(fix_na_var)

disease_data <- rbind.data.frame(disease2019,disease2018,disease2017,disease2016,disease2015,fix_na_var) %>% 
  distinct()
str(disease_data)
str(disease2019)

head(disease_data)
unique(disease_data$year)
unique(disease_data$site)
unique(disease_data$variety)
disease_data %>% filter(is.na(variety)) #2018 W + LLO already fixed, 2017 P is still TBD!!
unique(disease_data$trt)
unique(disease_data$timepoint)
unique(disease_data$avg_inc)
unique(disease_data$avg_sev)
unique(disease_data$date)

## NOTE: Need to add 2020 data
## NOTE: No variety for 2017 P !!

new_dat5 <- left_join(new_dat4,disease_data, by = c("year", "site", "variety", "trt", "timepoint"))
nrow(new_dat4)==nrow(new_dat5)
head(new_dat5)


# see if all date data ended up in new_dat5 dataset
disease_data %>% filter(!is.na(date)) %>% select(year,site,variety,trt,timepoint) %>%
  distinct() %>% nrow() #111
new_dat5 %>% filter(!is.na(date)) %>% select(year,site,variety,trt,timepoint) %>%
  distinct() %>% nrow() #90

sort(unique(disease_data$date))
sort(unique(new_dat5$date))

## Dates missing from new_dat5:
# "2019-08-16" "2019-08-20" "2019-09-18"
disease_data %>% filter(date == "2019-08-16") # 19 LL AU T1 had no disease
disease_data %>% filter(date == "2019-08-20") # 19 AR CH T2 had no disease
disease_data %>% filter(date == "2019-09-18") # 19 AR CH T3 had no disease

new_dat5 %>% filter(year=="2019" & site=="AR") %>% select(year,site,variety,trt,timepoint) %>% distinct()


# see if all disease data ended up in the new_dat5 dataset
disease_data %>% filter(!is.na(avg_inc)) %>% select(year,site,variety,trt,timepoint) %>%
  distinct() %>% nrow() #65
new_dat5 %>% filter(!is.na(avg_inc)) %>% select(year,site,variety,trt,timepoint) %>%
  distinct() %>% nrow() #64

sort(unique(disease_data$avg_inc))
sort(unique(new_dat5$avg_inc))

disease_data %>% filter(!is.na(avg_inc)) %>% select(year,site,variety,trt,timepoint) %>%
  distinct() %>% arrange(year,site,variety,trt,timepoint) #65
new_dat5 %>% filter(!is.na(avg_inc)) %>% select(year,site,variety,trt,timepoint) %>%
  distinct() %>% arrange(year,site,variety,trt,timepoint) #64, missing LL CH T1


## LATER: Add additional meta data (dates,who collected,etc)


new_dat6 <- new_dat5 %>% 
  mutate(trt_desc = case_when(
    year == "2019" & trt == "3" ~ "Rhyme2.08EC_and_Induce",
    year == "2019" & trt == "10" ~ "Luna_Exp._and_Induce",
    year == "2019" & trt == "11" ~ "Fervent_and_Induce",
    year == "2019" & trt == "13" ~ "Kenja_and_Induce",
    year == "2019" & trt == "18" ~ "UTC",
    year == "2018" & trt == "1" ~ "Rhyme2.08EC_and_Induce",
    year == "2018" & trt == "10" ~ "Luna_Exp._and_Induce",
    year == "2018" & trt == "11" ~ "Luna_Sensation_and_Induce",
    year == "2018" & trt == "15" ~ "Mettle125ME_and_Induce",
    year == "2018" & trt == "29" ~ "UTC_JMS_Stylet_Oil",
    year == "2017" & trt == "S" ~ "trt10_SDHI_fluopyram_and_low_tebuconazole",
    year == "2017" & trt == "D" ~ "trt25_DMI_tetraconazole",
    year == "2017" & trt == "E" ~ "trt12_DMI_flutriafol",
    year == "2017" & trt == "Q" ~ "trt13_QoI_azoxystrobin_and_flutriafol",
    year == "2017" & trt == "33" & str_detect(sample,"[1,2][S]")==TRUE ~ "trt33_UTC_single_colony",
    year == "2017" & trt == "33" & str_detect(sample,"[1,2][C]")==TRUE ~ "trt33_UTC_coaleascing_colony",
    year == "2016" & trt == "2" ~ "UBI-4319_4SC_and_Induce_weak_DMI",
    year == "2016" & trt == "29" ~ "29_UTC",
    year == "2016" & trt == "5" ~ "Revus_Top_and_Induce_moderate_DMI",
    year == "2015" & trt == "2Q" ~ "Quinoxyfen",
    year == "2015" & trt == "1D" ~ "DMI",
    year == "2015" & trt == "11S" ~ "SDHI_and_DMI",
    year == "2015" & trt == "16Q" ~ "QoI",
    year == "2015" & trt == "UTC" ~ "rows_2_5_9",
    TRUE ~ NA_character_
  ))

# Now we want to replace each missing value (coded as ./.:0 in the raw data) with a format that R understands
# and recognizes as missing data: NA
new_dat6[new_dat6$reads=="./.:0",]$reads <- NA_character_
new_dat6 %>% filter(reads =="./.:0")
head(new_dat6)

# Now we want to break down the reads column into four columns: 
# allele 1, allele 2, reads for allele 1, reads for allele 2
new_dat6$alleles <- word(new_dat6$reads, 1, 1, sep = ":")
new_dat6$all_reads <- word(new_dat6$reads, 2, 2, sep = ":")
new_dat6$all1 <- word(new_dat6$alleles, 1, 1, sep = "/")
new_dat6$reads1 <- word(new_dat6$all_reads, 1, 1, sep = ",")
new_dat6$all2 <- word(new_dat6$alleles, 2, 2, sep = "/")
new_dat6$reads2 <- word(new_dat6$all_reads, 2, 2, sep = ",")

# Sometimes when there is only one allele the second id is not an NA as it should be.
# Replace the second allele with NA every time the number of reads for the second allele is missing.
new_dat6[is.na(new_dat6$reads2),]$all2 <- NA

# Now we can drop the columns we no longer need and make reads1+2 numeric 
new_dat7 <- new_dat6 %>% 
  select(-reads, -alleles, -all_reads) %>% 
  mutate(sample_id_new = paste(run,year,site,variety,trt,timepoint,sample,plate,well,sep = "_")) %>% 
  mutate(reads1 = as.numeric(reads1)) %>% 
  mutate(reads2 = as.numeric(reads2))
head(new_dat7) #check
length(unique(new_dat7$sample_id))==length(unique(new_dat7$sample_id_new)) #differences but should have same number
unique(new_dat7$all1)
#new_dat7 %>% filter(is.na(all1)) %>% view()
#new_dat7 %>% filter(site=="LLO") %>% pull(sample_id) %>% unique() 
#new_dat7 %>% filter(site=="LLO") %>% pull(sample_id_new) %>% unique() 

# Check data types: keep everything a character
str(new_dat7)

# Save AmpSeq database as RDS and excel
write_rds(new_dat7,"Dissertation_Final_Draft/r_data/database.rds")
write_rds(new_dat7,"Dissertation_Final_Draft/output/database/database.rds")
write_csv(new_dat7,"Dissertation_Final_Draft/output/database/database.csv")

# Take a look at the data!
str(new_dat7)
unique(new_dat7$loci) #83 loci
unique(new_dat7$run) #9 runs
unique(new_dat7$year)
unique(new_dat7$site)
unique(new_dat7$variety)
unique(new_dat7$trt)
unique(new_dat7$timepoint)
unique(new_dat7$sample)
unique(new_dat7$plate)
unique(new_dat7$well)
unique(new_dat7$date)
unique(new_dat7$avg_inc)
unique(new_dat7$avg_sev)
unique(new_dat7$trt_desc)
unique(new_dat7$all1)
#max(new_dat7$reads1)
unique(new_dat7$all2)
#max(new_dat7$reads2)
length(unique(new_dat7$sample_id_new)) #5082 (NOTE:some were sequenced more than once!!!)


# Save AmpSeq database as RDS and excel
# This time without rows of missing loci/alleles 
new_dat7b <- new_dat7 %>% filter(!is.na(all1))
write_rds(new_dat7b,"Dissertation_Final_Draft/r_data/database_clean.rds")
write_rds(new_dat7b,"Dissertation_Final_Draft/output/database_clean.rds")
write_csv(new_dat7b,"Dissertation_Final_Draft/output/database_clean.csv")



## LATER (DATABASE):
## Add the GPS coordinates, esp for each of the wild samples

# Add files to folder
## CR trt maps, BH map
## Sampling summaries from me and Hema and team
## Sample files and primer pools





## Ch. 2 Appendix Table: Date of Collection for 2017-2020

setwd("../../../Desktop/Data/Tutoring")
new_dat7b <- readRDS("Dissertation_Final_Draft/r_data/database_clean.rds")

head(new_dat7b)

new_dat7b %>% filter(!is.na(date)) %>% pull(year) %>% unique()

unique(new_dat7b$run)
unique(new_dat7b$year)
unique(new_dat7b$site)
unique(new_dat7b$variety)

# only samples used in Ch. 2 analyses
sample_date1 <- new_dat7b %>% 
  filter(year %in% c(2015:2020) &
           !run %in% c('201512','201611-AV9') &
           !site %in% c('LLO','Control') &
           variety != "DE" &
           !is.na(variety)) %>% 
  select(year,site,variety,timepoint,date) %>% #no trt needed: all cr trts collected same day
  distinct() %>% 
  arrange(year,site,variety,timepoint)
sample_date1 #need to add dates for 2019 W, 2020 all

# add dates for 2019 W, 2020 all
sample_date2 <- sample_date1 %>% 
  mutate(date_new = as.character(date)) %>% 
  mutate(date_new = case_when(
    year == "2017" & site == "CR" ~ NA_character_,#"2017-01-01", #LATER: find
    year == "2017" & site == "BH" & variety == "CN19" ~ "2017-09-21",
    year == "2019" & site == "W" & timepoint == "T1" ~ "2019-08-24", #Aug. 23-25
    year == "2019" & site == "W" & timepoint == "T2" ~ "2019-09-30",
    year == "2020" & site == "BH" & variety == "AU" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "BH" & variety == "CN5" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "FR" & variety == "CH" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "LL" & variety == "AU" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "LL" & variety == "CH" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "W" & timepoint == "T1" ~ "2020-09-04-17",
    year == "2020" & site == "W" & timepoint == "T2" ~ "2020-09-24", #20200923-25
    TRUE ~ date_new
  ))
sample_date2 %>% filter(is.na(date_new))
head(sample_date2)

## LATER: FIX 2020 W TP CUTOFF:
## SAMPLES 1:41 - T1
## SAMPLES 42:117

## LATER: ADD ALL DATES TO APPENDIX TABLE

#rename site and variety
sample_date3 <- sample_date2 %>%  
  mutate(site_new=case_when(
    site == "AR" ~ "C-A",
    site == "BH" & variety == "AU" ~ "C-B1",
    site == "BH" & variety == "CN5"  ~ "C-B2",
    site == "BH" & variety == "CN19"  ~ "C-B3",
    site == "BH" & variety == "CN35"  ~ "C-B4",
    site == "FR"  & variety == "CH" ~ "C-C1",
    site == "FR"  & variety == "LM" ~ "C-C2",
    site == "LL" & variety == "AU" ~ "C-D1",
    site == "LL" & variety == "CH" ~ "C-D2",
    site == "LL" & variety == "PB" ~ "C-D3",
    site == "LL" & variety == "PN" ~ "C-D4",
    site == "CR" ~ "T",
    site == "W" ~ "W",
    TRUE ~ site
  )) %>% 
  mutate(variety_new=case_when(
    variety == "AU" ~ "Aurore",
    variety == "CH" ~ "Chardonnay",
    variety %in% c('CN19','CN35','CN5') ~ "Chancellor",
    variety == "LM" ~ "Lemberger",
    variety == "PB" ~ "Pinot blanc",
    variety == "PN" ~ "Pinot noir",
    variety == "VA" ~ "V. aestivalis",
    variety == "VL" ~ "V. labrusca",
    variety == "VR" ~ "V. riparia",
  )) %>% 
  select(year,site_new,variety_new,timepoint,date_new) %>% 
  distinct() 
head(sample_date3)

sample_date4 <- sample_date3 %>% 
  filter(date_new != "2015-10-01") %>% 
  mutate_at("timepoint",replace_na,"T1") %>% 
  pivot_wider(id_cols = c('year','site_new','variety_new'), names_from = "timepoint", values_from = "date_new") %>% 
  arrange(year,site_new) %>% 
  select(-variety_new)
view(sample_date4) #need to add trt or remove:
# 15 T 2015-10-01, 17 T ?

write_csv(sample_no5,"output/sample_dates_by_yr_trt.csv")





## Ch. 2 Appendix Table: Disease Inc and Sev for 2018 and 2019

# only samples used in Ch. 2 analyses
sample_dis1 <- new_dat7b %>% 
  filter(year %in% c(2015:2020) &
           !run %in% c('201512','201611-AV9') &
           !site %in% c('LLO','Control') &
           variety != "DE" &
           !is.na(variety)) %>% 
  select(year,site,variety,trt,timepoint,avg_inc,avg_sev) %>% #no trt needed: all cr trts collected same day
  distinct() %>% 
  arrange(year,site,variety,trt,timepoint)
sample_dis1 #need to add dates for 2019 W, 2020 all

# add dates for 2019 W, 2020 all
sample_date2 <- sample_date1 %>% 
  mutate(date_new = as.character(date)) %>% 
  mutate(date_new = case_when(
    year == "2017" & site == "CR" ~ NA_character_,#"2017-01-01", #LATER: find
    year == "2017" & site == "BH" & variety == "CN19" ~ "2017-09-21",
    year == "2019" & site == "W" & timepoint == "T1" ~ "2019-08-24", #Aug. 23-25
    year == "2019" & site == "W" & timepoint == "T2" ~ "2019-09-30",
    year == "2020" & site == "BH" & variety == "AU" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "BH" & variety == "CN5" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "FR" & variety == "CH" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "LL" & variety == "AU" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "LL" & variety == "CH" & timepoint == "T1" ~ NA_character_,#"2020-01-01",#LATER: find
    year == "2020" & site == "W" & timepoint == "T1" ~ "2020-09-04-17",
    year == "2020" & site == "W" & timepoint == "T2" ~ "2020-09-24", #20200923-25
    TRUE ~ date_new
  ))
sample_date2 %>% filter(is.na(date_new))
head(sample_date2)

## LATER: FIX 2020 W TP CUTOFF:
## SAMPLES 1:41 - T1
## SAMPLES 42:117

## LATER: ADD ALL DATES TO APPENDIX TABLE

#rename site and variety
sample_date3 <- sample_date2 %>%  
  mutate(site_new=case_when(
    site == "AR" ~ "C-A",
    site == "BH" & variety == "AU" ~ "C-B1",
    site == "BH" & variety == "CN5"  ~ "C-B2",
    site == "BH" & variety == "CN19"  ~ "C-B3",
    site == "BH" & variety == "CN35"  ~ "C-B4",
    site == "FR"  & variety == "CH" ~ "C-C1",
    site == "FR"  & variety == "LM" ~ "C-C2",
    site == "LL" & variety == "AU" ~ "C-D1",
    site == "LL" & variety == "CH" ~ "C-D2",
    site == "LL" & variety == "PB" ~ "C-D3",
    site == "LL" & variety == "PN" ~ "C-D4",
    site == "CR" ~ "T",
    site == "W" ~ "W",
    TRUE ~ site
  )) %>% 
  mutate(variety_new=case_when(
    variety == "AU" ~ "Aurore",
    variety == "CH" ~ "Chardonnay",
    variety %in% c('CN19','CN35','CN5') ~ "Chancellor",
    variety == "LM" ~ "Lemberger",
    variety == "PB" ~ "Pinot blanc",
    variety == "PN" ~ "Pinot noir",
    variety == "VA" ~ "V. aestivalis",
    variety == "VL" ~ "V. labrusca",
    variety == "VR" ~ "V. riparia",
  )) %>% 
  select(year,site_new,variety_new,timepoint,date_new) %>% 
  distinct() 
head(sample_date3)

sample_date4 <- sample_date3 %>% 
  filter(date_new != "2015-10-01") %>% 
  mutate_at("timepoint",replace_na,"T1") %>% 
  pivot_wider(id_cols = c('year','site_new','variety_new'), names_from = "timepoint", values_from = "date_new") %>% 
  arrange(year,site_new) %>% 
  select(-variety_new)
view(sample_date4) #need to add trt or remove:
# 15 T 2015-10-01, 17 T ?

write_csv(sample_no5,"output/sample_dates_by_yr_trt.csv")






