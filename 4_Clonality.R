
## Determine Sample Clonality

#Load packages
library(tidyverse)
library(readxl)
library(lubridate)


#Filtered data with column added for mating type (0-3) and clonality (0/1)
full_data <- readRDS("r_data/filtered_data_mat.rds")
str(full_data)
unique(full_data$run)
unique(full_data$year)
unique(full_data$site)



## Make tables with perc het for loci in different groups of samples



## Determine heterozygosity for each sample-locus pair
# Isolate each sample and assign het = TRUE if it has two alleles per locus, het = FALSE if it has one.
str(full_data)
count_het <- full_data %>%
  group_by(loci,sample_id) %>% 
  #sample_id is all we need here to account for all columns except for but all & reads
  tally() %>% 
  ungroup() %>%
  mutate(het = case_when(
    n == 2 ~ TRUE,
    n == 1 ~ FALSE
  ))
unique(count_het$n)

dat_joined <- full_data %>%
  full_join(count_het) %>% #left_join/full_join doesn't matter
  select(-n) %>% 
  mutate(days=yday(date))
head(dat_joined)
unique(dat_joined$het)
nrow(dat_joined)==nrow(full_data)
dat_joined %>% select(sample_id,run,year,site,variety,trt,timepoint,sample,plate,well,date,avg_inc,avg_sev,mat) %>% 
  distinct() %>% nrow()==length(unique(dat_joined$sample_id))


## Only use samples of interest for Ch. 2
unique(dat_joined$year)
unique(dat_joined$site)
dat_joined2 <- dat_joined %>% 
  filter(year %in% c(2015:2020,'Control')) %>% 
  filter(site %in% c('BH','LL','FR','W','CR','AR','Control'))
unique(dat_joined2$year)
unique(dat_joined2$site)

# All Ch. 2 samples
# Count true/het per locus
perc_het_true <- dat_joined2 %>%
  select(loci, sample_id, het) %>%
  distinct() %>%
  group_by(loci) %>%
  tally(het == "TRUE") %>%
  mutate(het = n) %>%
  ungroup() %>%
  select(loci, het)

#count false/hom per locus
perc_het_false <- dat_joined2 %>%
  select(loci, sample_id, het) %>%
  distinct() %>%
  group_by(loci) %>%
  tally(het == "FALSE") %>%
  mutate(hom = n) %>%
  ungroup() %>%
  select(loci, hom)

#join together  
perc_het <- perc_het_true %>%
  left_join(perc_het_false)

#then calculate T/T+F
perc_het_t <- perc_het %>%
  mutate(total = het + hom) %>%
  mutate(perc_het = (het/total)*100) %>%
  arrange(perc_het)
#view(perc_het_t)

#export csv
write.csv(perc_het_t, "output/perc_het_by_all_loci.csv")


## Same but remove controls
dat_joined3 <- dat_joined2 %>%
  filter(year != "Control" & site != "Control")
perc_het_true2 <- dat_joined3 %>%
  select(loci, sample_id, het) %>%
  distinct() %>% #can also use unique()
  group_by(loci) %>%
  tally(het == "TRUE") %>%
  mutate(het = n) %>%
  ungroup() %>%
  select(loci, het)
#count false/hom per locus
perc_het_false2 <- dat_joined3 %>%
  select(loci, sample_id, het) %>%
  distinct() %>%
  group_by(loci) %>%
  tally(het == "FALSE") %>%
  mutate(hom = n) %>%
  ungroup() %>%
  select(loci, hom)
#join together  
perc_het2 <- perc_het_true2 %>%
  left_join(perc_het_false2)
#then calculate T/T+F
perc_het_t2 <- perc_het2 %>%
  mutate(total = het + hom) %>%
  mutate(perc_het = (het/total)*100) %>%
  arrange(perc_het)
#view(perc_het_t2)
#export csv
write.csv(perc_het_t2, "output/perc_het_wout_controls.csv")


## only controls
dat_joined4 <- dat_joined2 %>%
  filter(year == "Control" | site == "Control")
unique(dat_joined4$sample_id) #45 samples
perc_het_true3 <- dat_joined4 %>%
  select(loci, sample_id, het) %>%
  distinct() %>% #can also use unique()
  group_by(loci) %>%
  tally(het == "TRUE") %>%
  mutate(het = n) %>%
  ungroup() %>%
  select(loci, het)
#count false/hom per locus
perc_het_false3 <- dat_joined4 %>%
  select(loci, sample_id, het) %>%
  distinct() %>%
  group_by(loci) %>%
  tally(het == "FALSE") %>%
  mutate(hom = n) %>%
  ungroup() %>%
  select(loci, hom)
#join together  
perc_het3 <- perc_het_true3 %>%
  left_join(perc_het_false3)
#then calculate T/T+F
perc_het_t3 <- perc_het3 %>%
  mutate(total = het + hom) %>%
  mutate(perc_het = (het/total)*100) %>%
  arrange(perc_het)
#view(perc_het_t3)
#export csv
write.csv(perc_het_t3, "output/perc_het_controls.csv")


## Only controls and samples from <=30% inc <=5% sev
dat_joined5 <- dat_joined2 %>%
  mutate(avg_inc = as.numeric(avg_inc)) %>% 
  mutate(avg_sev = as.numeric(avg_sev)) %>% 
  filter((avg_sev <= 5.0 &
           avg_inc <= 30.0 ) |
           (year=="Control"|site=="Control"))
unique(dat_joined5$avg_inc)
unique(dat_joined5$avg_sev)

perc_het_true4 <- dat_joined5 %>%
  select(loci, sample_id, het) %>%
  distinct() %>% #can also use unique()
  group_by(loci) %>%
  tally(het == "TRUE") %>%
  mutate(het = n) %>%
  ungroup() %>%
  select(loci, het)
#count false/hom per locus
perc_het_false4 <- dat_joined5 %>%
  select(loci, sample_id, het) %>%
  distinct() %>%
  group_by(loci) %>%
  tally(het == "FALSE") %>%
  mutate(hom = n) %>%
  ungroup() %>%
  select(loci, hom)
#join together  
perc_het4 <- perc_het_true4 %>%
  left_join(perc_het_false4)
#then calculate T/T+F
perc_het_t4 <- perc_het4 %>%
  mutate(total = het + hom) %>%
  mutate(perc_het = (het/total)*100) %>%
  arrange(perc_het)
#view(perc_het_t4)
#export csv
write.csv(perc_het_t4, "output/perc_het_controls_and_less30inc_and_less5sev.csv")


## All samples, but separate by year (no controls) (only 2018 and 2019)
#2018
dat_joined6 <- dat_joined2 %>%
  filter(year == "2018") %>% 
  filter(year!="Control"&site!="Control")
perc_het_true5 <- dat_joined6 %>%
  select(loci, sample_id, het) %>%
  distinct() %>% #can also use unique()
  group_by(loci) %>%
  tally(het == "TRUE") %>%
  mutate(het = n) %>%
  ungroup() %>%
  select(loci, het)
#count false/hom per locus
perc_het_false5 <- dat_joined6 %>%
  select(loci, sample_id, het) %>%
  distinct() %>%
  group_by(loci) %>%
  tally(het == "FALSE") %>%
  mutate(hom = n) %>%
  ungroup() %>%
  select(loci, hom)
#join together  
perc_het5 <- perc_het_true5 %>%
  left_join(perc_het_false5)
#then calculate T/T+F
perc_het_t5 <- perc_het5 %>%
  mutate(total = het + hom) %>%
  mutate(perc_het = (het/total)*100) %>%
  arrange(perc_het)
#view(perc_het_t5)
#export csv
write.csv(perc_het_t5, "output/perc_het_18.csv")

#2019
dat_joined7 <- dat_joined2 %>%
  filter(year == "2019") %>% 
  filter(year!="Control"&site!="Control")
perc_het_true6 <- dat_joined7 %>%
  select(loci, run, year, site, variety, trt, timepoint, sample, mat, het) %>%
  distinct() %>% #can also use unique()
  group_by(loci) %>%
  tally(het == "TRUE") %>%
  mutate(het = n) %>%
  ungroup() %>%
  select(loci, het)
#count false/hom per locus
perc_het_false6 <- dat_joined7 %>%
  select(loci, run, year, site, variety, trt, timepoint, sample, mat, het) %>%
  distinct() %>%
  group_by(loci) %>%
  tally(het == "FALSE") %>%
  mutate(hom = n) %>%
  ungroup() %>%
  select(loci, hom)
#join together  
perc_het6 <- perc_het_true6 %>%
  left_join(perc_het_false6)
#then calculate T/T+F
perc_het_t6 <- perc_het6 %>%
  mutate(total = het + hom) %>%
  mutate(perc_het = (het/total)*100) %>%
  arrange(perc_het)
#view(perc_het_t4)
#export csv
write.csv(perc_het_t6, "output/perc_het_19.csv")





## % het markers for each control sample
controls <- dat_joined2 %>% 
  filter(year=="Control"|site=="Control")
controls_het <- controls %>% 
  filter(het=="TRUE") %>% 
  group_by(run,sample,sample_id) %>% 
  tally() %>% 
  ungroup() %>% 
  select(run,sample,sample_id,het=n)
controls_hom <- controls %>% 
  filter(het=="FALSE") %>% 
  group_by(run,sample,sample_id) %>% 
  tally() %>% 
  ungroup() %>% 
  select(run,sample,sample_id,hom=n)
controls1 <- full_join(controls_het,controls_hom) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%  
  mutate(het = het/2) %>% #NOTE: het have 2 loci per sample so divide by 2!
  mutate(perc_het = (het/(het+hom))) %>% 
  arrange(perc_het)
view(controls1)
range(controls1$perc_het)
write_csv(controls1,"output/perc_het_by_control_sample_loci_not_listed.csv")

test1 <- controls %>% filter(het==TRUE) %>% select(run,sample,het_loci=loci) 
head(test1)
test2 <- controls1 %>% full_join(test1) %>% arrange(perc_het) %>% distinct()
head(test2)
write_csv(test2,"output/perc_het_by_control_sample_loci_listed.csv")

test2 %>% filter(!is.na(het_loci)) %>% pull(sample_id) %>% unique() %>% length()
#30 samples

test3 <- test2 %>% group_by(het_loci) %>% tally() %>% arrange(n) %>% 
  filter(!is.na(het_loci)) %>% 
  mutate(perc=n/30)
test3
write_csv(test3,"output/control_het_loci_count.csv")




## Ok now choose markers for determining clonality

view(perc_het_t2) #samples, without controls
view(perc_het_t3) #controls (isolates)
view(perc_het_t4) #controls and samples from <=30%/5% inc/sev


perc_het_t3 <- read_csv("../../../Desktop/Data/Tutoring/Dissertation_Final_Draft/output/perc_het_controls.csv")
perc_het_t2 <- read_csv("../../../Desktop/Data/Tutoring/Dissertation_Final_Draft/output/perc_het_wout_controls.csv")


6/23 #kept: 0.2609
#Walt G143A not included
t3loci <- perc_het_t3 %>% filter(perc_het < 7 & total > 22) %>% pull(loci) %>% unique()
t3loci #21 loci

controls_hom_loci <- perc_het_t2 %>% filter(perc_het > 10 & loci %in% t3loci) %>% 
  pull(loci) %>% unique()
controls_hom_loci #18 loci

hom_loci <- perc_het_t2 %>% filter(perc_het < 12.8) %>% pull(loci) %>% unique()
hom_loci #Make Walt_G143A the cutoff

# Make a list of known het loci
unique(dat_joined2$loci)
het_loci <- c('SDHb_h267y', 'Walt_G143A', 'CYP51_a495t_y136f', 
              'CYP51_a1119c', 'SDHbFH277RB', 'CYTb')
het_loci



## Calculate percent clonality
## Using controls_hom_loci and hom_loci
str(dat_joined2)
unique(dat_joined2$year)


# A. use controls_hom_loci (this is the new method of determining clonality)
# count number of loci amplified for each sample
perc_clonal_a1 <- dat_joined2 %>% 
  select(-all,-reads) %>% 
  distinct() %>% 
  filter(loci %in% controls_hom_loci & 
           !loci %in% het_loci) %>% 
  group_by(year,sample_id) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(tot_loc = as.numeric(n)) %>% 
  select(-n)
head(perc_clonal_a1)
unique(perc_clonal_a1$tot_loc)
range(perc_clonal_a1$tot_loc) #1-16

# count number of het loci for each sample
perc_clonal_a2 <- dat_joined2 %>% 
  select(-all,-reads) %>% 
  distinct() %>% 
  filter(loci %in% controls_hom_loci & 
           !loci %in% het_loci &
           het == "TRUE") %>% 
  group_by(year,sample_id) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(het_loc = as.numeric(n)) %>% 
  select(-n)
head(perc_clonal_a2)
unique(perc_clonal_a2$het_loc)
range(perc_clonal_a2$het_loc) #1-15

perc_clonal_a <- full_join(perc_clonal_a1,perc_clonal_a2) %>% 
  mutate_at(c('het_loc'), ~replace_na(.,0)) %>% 
  mutate(hom_loc = tot_loc-het_loc) %>% 
  mutate(perc_clonal = hom_loc/tot_loc) %>% 
  mutate(perc_clonal_simple = case_when(
    perc_clonal < .11 ~ .05,
    perc_clonal > .11 & perc_clonal < .21 ~ .15,
    perc_clonal > .21 & perc_clonal < .31 ~ .25,
    perc_clonal > .31 & perc_clonal < .41 ~ .35,
    perc_clonal > .41 & perc_clonal < .51 ~ .45,
    perc_clonal > .51 & perc_clonal < .61 ~ .55,
    perc_clonal > .61 & perc_clonal < .71 ~ .65,
    perc_clonal > .71 & perc_clonal < .81 ~ .75,
    perc_clonal > .81 & perc_clonal < .91 ~ .85,
    perc_clonal > .91 ~ .95
  ))
head(perc_clonal_a)
unique(perc_clonal_a$perc_clonal_simple)
range(perc_clonal_a$perc_clonal_simple) #0.05-0.95
#perc_clonal_a %>% filter(is.na(perc_clonal_simple)) %>% pull(sample_id) %>% unique()
unique(perc_clonal_a$perc_clonal)
range(perc_clonal_a$perc_clonal) #0-1
range(perc_clonal_a$tot_loc) #1-16
perc_clonal_a %>% filter(year=="2018") %>% pull(tot_loc) %>% range() #1-8
perc_clonal_a %>% filter(year=="2019") %>% pull(tot_loc) %>% range() #1-16

test_clonal <- perc_clonal_a %>% 
  filter(tot_loc > 5) %>% 
  #group_by(year,perc_clonal_simple) %>% tally() %>% 
  group_by(perc_clonal_simple) %>% tally()
test_clonal #there shouldn't be any NA's
# control samples are 85 and 95



# Add percent clonality data back to the full dataset
head(perc_clonal_a)
dat_joined_clonal <- full_join(dat_joined2,perc_clonal_a, by = c("sample_id", "year"))
nrow(dat_joined2)==nrow(dat_joined_clonal)
head(dat_joined_clonal)
unique(dat_joined_clonal$perc_clonal_simple)
unique(dat_joined_clonal$tot_loc)

str(dat_joined_clonal)

dat_joined_clonal <- dat_joined_clonal %>% 
  mutate(high_conf_clonal = case_when(
    perc_clonal < 0.3 & tot_loc > 6 ~ "0", #at least 5/7 loci het
    perc_clonal > 0.7 & tot_loc > 6 ~ "1" #at least 5/7 loci hom
)) 
unique(dat_joined_clonal$high_conf_clonal)

unique(dat_joined_clonal$run)
unique(dat_joined_clonal$year)
unique(dat_joined_clonal$site)

str(dat_joined_clonal)
## Export 2015-2020 & Control dataset
# cleaned, filtered, mat, clonality
write_rds(dat_joined_clonal,"r_data/filtered_data_mat_clonal.rds")



