

## Now figure out how to format the SSR data for pop stats
# https://grunwaldlab.github.io/Population_Genetics_in_R/Data_Preparation.html

## LATER: KEEP MORE LOCI AND SAMPLES???
ssr_dat_combined <- readRDS("r_data/ssr_dat.rds")
str(ssr_dat_combined)

ssr_dat_fmt <- ssr_dat_combined %>% 
  mutate(run = word(sample_id,1,1,"_")) %>%
  mutate(year = word(sample_id,2,2,"_")) %>%
  mutate(site = word(sample_id,3,3,"_")) %>%
  mutate(variety = word(sample_id,4,4,"_")) %>%
  mutate(trt = word(sample_id,5,5,"_")) %>%
  mutate(timepoint = word(sample_id,6,6,"_")) %>%
  mutate(sample = word(sample_id,6,6,"_")) %>%
  mutate(variety_simple = case_when(
    variety %in% c('CN5','CN19','CN35') ~ "CN",
    TRUE ~ variety
  )) %>% 
  filter(year!="Control"&site!="Control") %>% #remove controls
  filter(variety!="NA") #remove 1 random W NA variety sample
  
ssr_dat_fmt %>% select(site,variety) %>% distinct() %>% arrange(site,variety)

ssr_dat_fmt2 <- ssr_dat_fmt %>% 
  mutate(site = case_when(
    site == "AR" ~ "A",
    site == "BH" & variety == "AU" ~ "B1",
    site == "BH" & variety == "CN5"  ~ "B2",
    site == "BH" & variety == "CN19"  ~ "B3",
    site == "BH" & variety == "CN35"  ~ "B4",
    site == "FR"  & variety == "CH" ~ "C1",
    site == "FR"  & variety == "LM" ~ "C2",
    site == "LL" & variety == "AU" ~ "D1",
    site == "LL" & variety == "CH" ~ "D2",
    site == "LL" & variety == "PB" ~ "D3",
    site == "LL" & variety == "PN" ~ "D4",
    site == "CR" ~ "T",
    site == "W" ~ "W",
    TRUE ~ site
  )) %>% 
  mutate(ctw = case_when(
    str_detect(site,"[A-D]") ~ "C",
    str_detect(site,"T") ~ "T",
    str_detect(site,"W") ~ "W"
  )) %>% 
  select(run,year,ctw,site,variety,variety_simple,trt,timepoint,sample_id,loci,all,len)

head(ssr_dat_fmt2)
unique(ssr_dat_fmt2$run)
unique(ssr_dat_fmt2$year)
unique(ssr_dat_fmt2$ctw)
unique(ssr_dat_fmt2$site)
unique(ssr_dat_fmt2$variety)
unique(ssr_dat_fmt2$variety_simple)
unique(ssr_dat_fmt2$trt)
unique(ssr_dat_fmt2$timepoint)


# Group pop by everything
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",trt,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(sample_id,pop,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3)
#view(ssr_dat_fmt3)
write_csv(ssr_dat_fmt3,"output/poppr/data/ssr_dat.csv")

# sample no. per pop in this dataset
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/poppr/data/ssr_dat_sample_n.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #4616 samples
length(unique(ssr_dat_fmt3$pop)) #104 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name





## remove pops with < 10 samples
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop) %>% 
  filter(n >= 10)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/poppr/data/ssr_dat_sample_n_10min.csv")

keep <- ssr_dat_fmt3_sample_n %>% pull(pop) %>% unique()
ssr_dat_fmt3a <- ssr_dat_fmt3 %>% filter(pop %in% keep)
head(ssr_dat_fmt3a,5)
write_csv(ssr_dat_fmt3a,"output/poppr/data/ssr_dat_10min.csv")

# Open both in excel and add: loci no, sample no, pop no, then paste popname/nums
head(ssr_dat_fmt3a)
length(colnames(ssr_dat_fmt3a))-2 #13 loci
length(unique(ssr_dat_fmt3a$sample_id)) #4562 samples
length(unique(ssr_dat_fmt3a$pop)) #91 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name


view(unique(ssr_dat_fmt3a$sample_id))


# Make subset test dataset - only 2018 + 2019 A, B1, B2, T, W, T1, T4
unique(ssr_dat_fmt2$site)
ssr_dat_fmt_subset <- ssr_dat_fmt2 %>% 
  filter(year %in% c('2018','2019') &
           site %in% c('A','B1','B2','T','W') &
           timepoint %in% c('T1','T4') &
           trt %in% c('29','18','NA')) %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",trt,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(sample_id,pop,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt_subset)

# sample no. per pop in this dataset
ssr_dat_fmt_subset_sample_n <- ssr_dat_fmt_subset %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  filter(n>9) %>% 
  select(n,pop)
ssr_dat_fmt_subset_sample_n
sum(ssr_dat_fmt_subset_sample_n$n) #519
range(ssr_dat_fmt_subset_sample_n$n)
write_csv(ssr_dat_fmt_subset_sample_n,"output/poppr/data/ssr_dat_subset_sample_n.csv")

keep <- ssr_dat_fmt_subset_sample_n %>% pull(pop)

ssr_dat_fmt_subset <- ssr_dat_fmt_subset %>% 
  filter(pop %in% keep)
head(ssr_dat_fmt_subset)
write_csv(ssr_dat_fmt_subset,"output/poppr/data/ssr_dat_subset.csv")

length(colnames(ssr_dat_fmt_subset))-2 #13 loci
length(unique(ssr_dat_fmt_subset$sample_id)) #519 samples
length(unique(ssr_dat_fmt_subset$pop)) #11 populations











## LOOK AT C vs T (include CR trts)
## LINE 408 - 600??
## Remove samples that don't work in at least 7/13
## Remove samples that were collected in Sept + Oct
ssr_dat_fmt4 <- ssr_dat_fmt2 %>% 
  filter(run != "201512" & 
           str_detect(year,"201[8,9]") &
           site != "TBD" & site != "P" & site != "LLO" &
           variety != "DE" & !is.na(variety) & variety != "NA" &
           site != "W"&
           str_detect(sample_id,"2018_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2018_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2018_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T2")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_CR_CH_[.]+_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T4")==FALSE) %>% 
  mutate(variety = case_when(
    variety == "LM" ~ "LB",
    variety == "CN5" ~ "CN",
    variety == "CNBH" ~ "CN",
    TRUE ~ variety
  ))  %>%
  mutate(site = case_when(
    site == "AR" ~ "A",
    site == "BH" & variety == "AU" ~ "B1",
    site == "BH" & variety == "CN"  ~ "B2",
    site == "BH" & variety == "CNB"  ~ "B3",
    site == "FR" & variety == "CH" ~ "C1",
    site == "FR" & variety == "LB" ~ "C2",
    site == "LL" & variety == "AU" ~ "D1",
    site == "LL" & variety == "CH" ~ "D2",
    site == "LL" & variety == "PB" ~ "D3",
    site == "LL" & variety == "PN" ~ "D4",
    site == "CR" ~ "T",
    site == "W" ~ "W",
    TRUE ~ site
  )) %>% 
  mutate(ctw = case_when(
    str_detect(site,"[A-D]") ~ "C",
    str_detect(site,"T") ~ "T",
  )) %>% 
  mutate(trt = case_when(
    year == "2019" & trt == "3" ~ "D",
    year == "2019" & trt == "10" ~ "DS",
    year == "2019" & trt == "11" ~ "DS",
    year == "2019" & trt == "13" ~ "S",
    year == "2019" & trt == "18" ~ "UTC",
    year == "2018" & trt == "1" ~ "D",
    year == "2018" & trt == "10" ~ "DS",
    year == "2018" & trt == "11" ~ "SQ",
    year == "2018" & trt == "15" ~ "D",
    year == "2018" & trt == "29" ~ "UTCj",
    year == "2017" & trt == "S" ~ "DS",
    year == "2017" & trt == "D" ~ "D",
    year == "2017" & trt == "E" ~ "D",
    year == "2017" & trt == "Q" ~ "DQ",
    year == "2017" & trt == "33" & str_detect(sample_id,"[1,2][S]") ~ "UTCs",
    year == "2017" & trt == "33" & str_detect(sample_id,"[1,2][C]") ~ "UTCc",
    year == "2016" & trt == "2" ~ "D",
    year == "2016" & trt == "29" ~ "UTC",
    year == "2016" & trt == "5" ~ "D",
    year == "2015" & trt == "2Q" ~ "Quin",
    year == "2015" & trt == "1D" ~ "D",
    year == "2015" & trt == "11S" ~ "DS",
    year == "2015" & trt == "16Q" ~ "Q",
    year == "2015" & trt == "UTC" ~ "UTC",
    TRUE ~ NA_character_
  )) %>% 
  select(year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  distinct() 

head(ssr_dat_fmt2)
#unique(ssr_dat_fmt2$run)
unique(ssr_dat_fmt2$year)
unique(ssr_dat_fmt2$ctw)
unique(ssr_dat_fmt2$site)
unique(ssr_dat_fmt2$variety)
unique(ssr_dat_fmt2$trt)
unique(ssr_dat_fmt2$timepoint)


ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct() 
head(ssr_dat_fmt3)

names(ssr_dat_fmt3)
sample_discard <- ssr_dat_fmt3 %>% 
  select(-pop) %>% 
  group_by(sample_id) %>% 
  summarise(is.na(EnMS2),
            is.na(EnCDnew_24),
            #is.na(EnMS11),
            is.na(EnMS6),
            is.na(EnCDnew_25),
            is.na(EnCDnew_29),
            is.na(EnCDnew_30),
            is.na(EnCDnew_32),
            is.na(EnCDnew_39),
            is.na(EnCDnew_41),
            is.na(EnMS8),
            is.na(EnCDnew_27),
            is.na(EnCDnew_35)
  ) %>% 
  ungroup() %>% 
  select(sample_id,L1='is.na(EnMS2)',
         L2='is.na(EnCDnew_24)',
         #L3='is.na(EnMS11)',
         L4='is.na(EnMS6)',
         L5='is.na(EnCDnew_25)',
         L6='is.na(EnCDnew_29)',
         L7='is.na(EnCDnew_30)',
         L8='is.na(EnCDnew_32)',
         L9='is.na(EnCDnew_39)',
         L10='is.na(EnCDnew_41)',
         L11='is.na(EnMS8)',
         L12='is.na(EnCDnew_27)',
         L13='is.na(EnCDnew_35)') %>% 
  mutate(L1=as.character(L1),
         L2=as.character(L2),
         #L3=as.character(L3),
         L4=as.character(L4),
         L5=as.character(L5),
         L6=as.character(L6),
         L7=as.character(L7),
         L8=as.character(L8),
         L9=as.character(L9),
         L10=as.character(L10),
         L11=as.character(L11),
         L12=as.character(L12),
         L13=as.character(L13)) %>% 
  mutate(L1=str_replace(L1,"TRUE","1")) %>% 
  mutate(L2=str_replace(L2,"TRUE","1")) %>% 
  #mutate(L3=str_replace(L3,"TRUE","1")) %>% 
  mutate(L4=str_replace(L4,"TRUE","1")) %>% 
  mutate(L5=str_replace(L5,"TRUE","1")) %>% 
  mutate(L6=str_replace(L6,"TRUE","1")) %>% 
  mutate(L7=str_replace(L7,"TRUE","1")) %>% 
  mutate(L8=str_replace(L8,"TRUE","1")) %>% 
  mutate(L9=str_replace(L9,"TRUE","1")) %>% 
  mutate(L10=str_replace(L10,"TRUE","1")) %>% 
  mutate(L11=str_replace(L11,"TRUE","1")) %>% 
  mutate(L12=str_replace(L12,"TRUE","1")) %>% 
  mutate(L13=str_replace(L13,"TRUE","1")) %>% 
  mutate(L1=str_replace(L1,"FALSE","0")) %>% 
  mutate(L2=str_replace(L2,"FALSE","0")) %>% 
  #mutate(L3=str_replace(L3,"FALSE","0")) %>% 
  mutate(L4=str_replace(L4,"FALSE","0")) %>% 
  mutate(L5=str_replace(L5,"FALSE","0")) %>% 
  mutate(L6=str_replace(L6,"FALSE","0")) %>% 
  mutate(L7=str_replace(L7,"FALSE","0")) %>% 
  mutate(L8=str_replace(L8,"FALSE","0")) %>% 
  mutate(L9=str_replace(L9,"FALSE","0")) %>% 
  mutate(L10=str_replace(L10,"FALSE","0")) %>% 
  mutate(L11=str_replace(L11,"FALSE","0")) %>% 
  mutate(L12=str_replace(L12,"FALSE","0")) %>% 
  mutate(L13=str_replace(L13,"FALSE","0")) %>% 
  mutate(L1=as.numeric(L1),
         L2=as.numeric(L2),
         #L3=as.numeric(L3),
         L4=as.numeric(L4),
         L5=as.numeric(L5),
         L6=as.numeric(L6),
         L7=as.numeric(L7),
         L8=as.numeric(L8),
         L9=as.numeric(L9),
         L10=as.numeric(L10),
         L11=as.numeric(L11),
         L12=as.numeric(L12),
         L13=as.numeric(L13)) %>% 
  #mutate(na_sum=L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13) %>% 
  mutate(na_sum=L1+L2+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13) %>% 
  filter(na_sum>7) %>% 
  pull(sample_id) %>% 
  unique()
sample_discard

# Group pop by everything
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  filter(!sample_id %in% sample_discard) %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",trt,"_",timepoint)) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3)
write_csv(ssr_dat_fmt3,"output/poppr/data/ssr_dat_final_tvc.csv")

# sample no. per pop in this dataset
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/poppr/data/ssr_dat_final_sample_n_tvc.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #1091 samples
length(unique(ssr_dat_fmt3$pop)) #33 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name


## repeat but remove pops with < 10 samples
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop) %>% 
  filter(n >= 10)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/poppr/data/ssr_dat_final_sample_n_10min_tvc.csv")
keep <- ssr_dat_fmt3_sample_n %>% pull(pop) %>% unique()
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  filter(!sample_id %in% sample_discard) %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",trt,"_",timepoint)) %>% 
  select(sample_id,pop,loci,len)  %>% 
  filter(pop %in% keep) %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3,5)
write_csv(ssr_dat_fmt3,"output/poppr/data/ssr_dat_final_10min_tvc.csv")
# Open both in excel and add: loci no, sample no, pop no, then paste popname/nums
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #1055 samples
length(unique(ssr_dat_fmt3$pop)) #26 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name















# Now group pop up to variety (no trt and tp)
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3)
write_csv(ssr_dat_fmt3,"output/poppr/data/ssr_dat_final_ycsv.csv")

# sample no. per pop in this dataset
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/poppr/data/ssr_dat_final_ycsv_sample_n.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #1587 samples
length(unique(ssr_dat_fmt3$pop)) #12 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name

## repeat but remove pops with < 10 samples
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop) %>% 
  filter(n >= 10)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/poppr/data/ssr_dat_final_ycsv_sample_n_10min.csv")
keep <- ssr_dat_fmt3_sample_n %>% pull(pop) %>% unique()
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  filter(pop %in% keep) %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3,5)
write_csv(ssr_dat_fmt3,"output/poppr/data/ssr_dat_final_ycsv_10min.csv")
# Open both in excel and add: loci no, sample no, pop no, then paste popname/nums
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #1577 samples
length(unique(ssr_dat_fmt3$pop)) #10 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name




# Now group pop up to site (no var, trt, tp)
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3)
write_csv(ssr_dat_fmt3,"output/poppr/data/ssr_dat_final_ycs.csv")

# sample no. per pop in this dataset
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/poppr/data/ssr_dat_final_ycs_sample_n.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #1587 samples
length(unique(ssr_dat_fmt3$pop)) #12 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name

## repeat but remove pops with < 10 samples
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop) %>% 
  filter(n >= 10)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/poppr/data/ssr_dat_final_ycs_sample_n_10min.csv")
keep <- ssr_dat_fmt3_sample_n %>% pull(pop) %>% unique()
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  filter(pop %in% keep) %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3,5)
write_csv(ssr_dat_fmt3,"output/poppr/data/ssr_dat_final_ycs_10min.csv")
# Open both in excel and add: loci no, sample no, pop no, then paste popname/nums
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #1577 samples
length(unique(ssr_dat_fmt3$pop)) #10 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name






# Now group pop up to ctw (no site, var, trt, tp)
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  mutate(pop = paste0(year,"_",ctw)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3)
write_csv(ssr_dat_fmt3,"output/ssr_dat_final_yc.csv")

# sample no. per pop in this dataset
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/ssr_dat_final_yc_sample_n.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #4555 samples
length(unique(ssr_dat_fmt3$pop)) #13 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name

## repeat but remove pops with < 10 samples
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop) %>% 
  filter(n >= 10)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/ssr_dat_final_yc_sample_n_10min.csv")
keep <- ssr_dat_fmt3_sample_n %>% pull(pop) %>% unique()
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  mutate(pop = paste0(year,"_",ctw)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  filter(pop %in% keep) %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3,5)
write_csv(ssr_dat_fmt3,"output/ssr_dat_final_yc_10min.csv")
# Open both in excel and add: loci no, sample no, pop no, then paste popname/nums
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #4554 samples
length(unique(ssr_dat_fmt3$pop)) #13 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name











## Let's just compare commercial vs wild 2018 + 2019
# Now group pop up to site (no var, trt, tp)
head(ssr_dat_fmt2)
ssr_dat_fmt2 %>% filter(str_detect(year,"201[8,9]")) %>% filter(site!="T") %>% 
  select(year,ctw,site,variety,trt,timepoint,sample_id) %>% distinct() %>% 
  group_by(year,ctw,site,variety,trt,timepoint) %>% tally() %>% ungroup %>% view()
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  filter(str_detect(year,"201[8,9]") & site!="T") %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3)
write_csv(ssr_dat_fmt3,"output/ssr_dat_final_1819ycsvtp.csv")

# sample no. per pop in this dataset
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/ssr_dat_final_1819ycsvtp_sample_n.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #1992 samples
length(unique(ssr_dat_fmt3$pop)) #42 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name

## repeat but remove pops with < 10 samples
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop) %>% 
  filter(n >= 10)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/ssr_dat_final_1819ycsvtp_sample_n_10min.csv")
keep <- ssr_dat_fmt3_sample_n %>% pull(pop) %>% unique()
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  filter(str_detect(year,"201[8,9]") & site!="T") %>%
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(sample_id,pop,loci,len)  %>% 
  filter(pop %in% keep) %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3,5)
write_csv(ssr_dat_fmt3,"output/ssr_dat_final_1819ycsvtp_10min.csv")
# Open both in excel and add: loci no, sample no, pop no, then paste popname/nums
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #1965 samples
length(unique(ssr_dat_fmt3$pop)) #36 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name







## Let's just compare commercial vs wild 2018 + 2019
## Remove late commercial timepoints bc fungicides prb stop being applied?

# Now group pop up to site (no var, trt, tp)
test <- ssr_dat_fmt2 %>% filter(str_detect(year,"201[,89]")&site!="T") %>% mutate(test=word(sample_id,2,6,"_"))
unique(test$test)

head(ssr_dat_fmt2)
ssr_dat_fmt2 %>% filter(str_detect(year,"201[8,9]") & site!="T") %>% 
  filter(str_detect(sample_id,"2018_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2018_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2018_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T2")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_CR_CH_[.]+_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T4")==FALSE
  ) %>% 
  select(year,ctw,site,variety,trt,timepoint,sample_id) %>% distinct() %>% 
  group_by(year,ctw,site,variety,trt,timepoint) %>% tally() %>% ungroup %>% view()
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  filter(str_detect(year,"201[8,9]") & site!="T") %>% 
  filter(str_detect(sample_id,"2018_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2018_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2018_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T2")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_CR_CH_[.]+_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T4")==FALSE
  ) %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3)
write_csv(ssr_dat_fmt3,"output/ssr_dat_final_1819ycsvtp_NOSEPT.csv")

# sample no. per pop in this dataset
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/ssr_dat_final_1819ycsvtp_NOSEPT_sample_n.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #795 samples
length(unique(ssr_dat_fmt3$pop)) #24 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name

## repeat but remove pops with < 10 samples
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop) %>% 
  filter(n >= 10)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/ssr_dat_final_1819ycsvtp_NOSEPT_sample_n_10min.csv")
keep <- ssr_dat_fmt3_sample_n %>% pull(pop) %>% unique()
ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  filter(str_detect(year,"201[8,9]") & site!="T") %>%
  filter(str_detect(sample_id,"2018_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2018_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2018_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T2")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_CR_CH_[.]+_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T4")==FALSE
  ) %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(sample_id,pop,loci,len)  %>% 
  filter(pop %in% keep) %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3,5)
write_csv(ssr_dat_fmt3,"output/ssr_dat_final_1819ycsvtp_NOSEPT_10min.csv")
# Open both in excel and add: loci no, sample no, pop no, then paste popname/nums
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #768 samples
length(unique(ssr_dat_fmt3$pop)) #18 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name












## Let's just compare commercial vs wild 2018 + 2019
## Remove late commercial timepoints bc fungicides prb stop being applied?
## Remove samples that did not amplify in at least 8/13 loci

head(ssr_dat_fmt2)
ssr_dat_fmt2 %>% filter(str_detect(year,"201[8,9]") & site!="T") %>% 
  filter(str_detect(sample_id,"2018_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2018_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2018_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T2")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_CR_CH_[.]+_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T4")==FALSE
  ) %>% 
  select(year,ctw,site,variety,trt,timepoint,sample_id) %>% distinct() %>% 
  group_by(year,ctw,site,variety,trt,timepoint) %>% tally() %>% ungroup %>% view()

ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  filter(str_detect(year,"201[8,9]") & site!="T") %>% 
  filter(str_detect(sample_id,"2018_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2018_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2018_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T2")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_CR_CH_[.]+_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T4")==FALSE&
           !sample_id %in% sample_discard
  ) %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct() 
head(ssr_dat_fmt3)

names(ssr_dat_fmt3)
sample_discard <- ssr_dat_fmt3 %>% 
  select(-pop) %>% 
  group_by(sample_id) %>% 
  summarise(is.na(EnMS2),
            is.na(EnCDnew_24),
            is.na(EnMS11),
            is.na(EnMS6),
            is.na(EnCDnew_25),
            is.na(EnCDnew_29),
            is.na(EnCDnew_30),
            is.na(EnCDnew_32),
            is.na(EnCDnew_39),
            is.na(EnCDnew_41),
            is.na(EnMS8),
            is.na(EnCDnew_27),
            is.na(EnCDnew_35)
  ) %>% 
  ungroup() %>% 
  select(sample_id,L1='is.na(EnMS2)',
         L2='is.na(EnCDnew_24)',
         L3='is.na(EnMS11)',
         L4='is.na(EnMS6)',
         L5='is.na(EnCDnew_25)',
         L6='is.na(EnCDnew_29)',
         L7='is.na(EnCDnew_30)',
         L8='is.na(EnCDnew_32)',
         L9='is.na(EnCDnew_39)',
         L10='is.na(EnCDnew_41)',
         L11='is.na(EnMS8)',
         L12='is.na(EnCDnew_27)',
         L13='is.na(EnCDnew_35)') %>% 
  mutate(L1=as.character(L1),
         L2=as.character(L2),
         L3=as.character(L3),
         L4=as.character(L4),
         L5=as.character(L5),
         L6=as.character(L6),
         L7=as.character(L7),
         L8=as.character(L8),
         L9=as.character(L9),
         L10=as.character(L10),
         L11=as.character(L11),
         L12=as.character(L12),
         L13=as.character(L13)) %>% 
  mutate(L1=str_replace(L1,"TRUE","1")) %>% 
  mutate(L2=str_replace(L2,"TRUE","1")) %>% 
  mutate(L3=str_replace(L3,"TRUE","1")) %>% 
  mutate(L4=str_replace(L4,"TRUE","1")) %>% 
  mutate(L5=str_replace(L5,"TRUE","1")) %>% 
  mutate(L6=str_replace(L6,"TRUE","1")) %>% 
  mutate(L7=str_replace(L7,"TRUE","1")) %>% 
  mutate(L8=str_replace(L8,"TRUE","1")) %>% 
  mutate(L9=str_replace(L9,"TRUE","1")) %>% 
  mutate(L10=str_replace(L10,"TRUE","1")) %>% 
  mutate(L11=str_replace(L11,"TRUE","1")) %>% 
  mutate(L12=str_replace(L12,"TRUE","1")) %>% 
  mutate(L13=str_replace(L13,"TRUE","1")) %>% 
  mutate(L1=str_replace(L1,"FALSE","0")) %>% 
  mutate(L2=str_replace(L2,"FALSE","0")) %>% 
  mutate(L3=str_replace(L3,"FALSE","0")) %>% 
  mutate(L4=str_replace(L4,"FALSE","0")) %>% 
  mutate(L5=str_replace(L5,"FALSE","0")) %>% 
  mutate(L6=str_replace(L6,"FALSE","0")) %>% 
  mutate(L7=str_replace(L7,"FALSE","0")) %>% 
  mutate(L8=str_replace(L8,"FALSE","0")) %>% 
  mutate(L9=str_replace(L9,"FALSE","0")) %>% 
  mutate(L10=str_replace(L10,"FALSE","0")) %>% 
  mutate(L11=str_replace(L11,"FALSE","0")) %>% 
  mutate(L12=str_replace(L12,"FALSE","0")) %>% 
  mutate(L13=str_replace(L13,"FALSE","0")) %>% 
  mutate(L1=as.numeric(L1),
         L2=as.numeric(L2),
         L3=as.numeric(L3),
         L4=as.numeric(L4),
         L5=as.numeric(L5),
         L6=as.numeric(L6),
         L7=as.numeric(L7),
         L8=as.numeric(L8),
         L9=as.numeric(L9),
         L10=as.numeric(L10),
         L11=as.numeric(L11),
         L12=as.numeric(L12),
         L13=as.numeric(L13)) %>% 
  mutate(na_sum=L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13) %>% 
  filter(na_sum>7) %>% 
  pull(sample_id) %>% 
  unique()
sample_discard


ssr_dat_fmt3 <- ssr_dat_fmt2 %>% 
  filter(str_detect(year,"201[8,9]") & site!="T") %>% 
  filter(str_detect(sample_id,"2018_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2018_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2018_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2018_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T2")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T2")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T3")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_CR_CH_[.]+_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T3")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T3")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_AR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN5_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_CN_NA_T4")==FALSE&
           str_detect(sample_id,"2019_BH_AU_NA_T4")==FALSE&
           str_detect(sample_id,"2019_FR_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_CH_NA_T4")==FALSE&
           str_detect(sample_id,"2019_LL_AU_NA_T4")==FALSE&
           !sample_id %in% sample_discard
  ) %>% 
  mutate(pop = paste0(year,"_",ctw,"_",site,"_",variety,"_",timepoint)) %>% 
  #select(pop,year,ctw,site,variety,trt,timepoint,sample_id,loci,all,len) %>% 
  select(pop,sample_id,loci,len)  %>% 
  pivot_wider(names_from = loci, values_from = len) %>% 
  distinct()
head(ssr_dat_fmt3)
write_csv(ssr_dat_fmt3,"output/ssr_dat_final_1819ycsvtp_7-13loci_NOSEPT.csv")

# sample no. per pop in this dataset
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/ssr_dat_final_1819ycsvtp_7-13loci_NOSEPT_sample_n.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #795 samples
length(unique(ssr_dat_fmt3$pop)) #24 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name

## repeat but remove pops with < 10 samples
ssr_dat_fmt3_sample_n <- ssr_dat_fmt3 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop) %>% 
  filter(n >= 10)
head(ssr_dat_fmt3_sample_n)
range(ssr_dat_fmt3_sample_n$n)
write_csv(ssr_dat_fmt3_sample_n,"output/ssr_dat_final_1819ycsvtp_7-13loci_NOSEPT_sample_n_10min.csv")
keep <- ssr_dat_fmt3_sample_n %>% pull(pop) %>% unique()
head(ssr_dat_fmt3)
ssr_dat_fmt3 <- ssr_dat_fmt3 %>% 
  filter(pop %in% keep)
head(ssr_dat_fmt3,5)
write_csv(ssr_dat_fmt3,"output/ssr_dat_final_1819ycsvtp_7-13loci_NOSEPT_10min.csv")
# Open both in excel and add: loci no, sample no, pop no, then paste popname/nums
head(ssr_dat_fmt3)
length(colnames(ssr_dat_fmt3))-2 #13 loci
length(unique(ssr_dat_fmt3$sample_id)) #768 samples
length(unique(ssr_dat_fmt3$pop)) #18 populations
# for samples nos per pop and corresponding pop names
# save with "_genalex" pasted at end of file name










## Export all data now??
ssr_dat_fmt_wcc <- read_csv("output/ssr_dat_fmt.csv")
ssr_dat_fmt_wcc %>% 
  mutate(year = word(pop,2,2,"_")) %>% 
  pull(year) %>% 
  unique() #"2019" "2015" "2017" "2016"

test <- ssr_dat_fmt_wcc %>% 
  mutate(year = word(pop,2,2,"_")) %>% 
  filter(year == "2019") %>% 
  mutate(site = word(pop,3,3,"_")) %>%
  mutate(wcc = case_when(
    site %in% c('BH','FR','LL','AR') ~ "CM",
    site == "CR" ~ "CR",
    site == "W" ~ "W"
  )) %>% 
  mutate(pop = paste0("2019_",wcc,"_",(word(pop,3,-1,"_")))) %>% 
  select(-year,-site,-wcc) %>% 
  mutate(mlg = paste0(EnMS2,"-",EnCDnew_24,"-",EnCDnew_37,"-",EnMS6,"-",EnCDnew_29,"-",EnCDnew_30,"-",
                      EnCDnew_32,"-",EnCDnew_33,"-",EnCDnew_39,"-",EnCDnew_41,"-",
                      EnMS8,"-",EnCDnew_27,"-",EnCDnew_25,"-",EnMS10))
test2 <- test %>% 
  mutate(pop_new = word(pop,1,2,"_")) %>% 
  select(id,pop_new,mlg) %>% 
  arrange(pop_new,mlg)
view(test2)

test3 <- test2 %>% 
  group_by(pop_new,mlg) %>% 
  tally() %>% 
  ungroup()
view(test3)

test4 <- test %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup()
view(test4)

#Maybe only compare W VL and FR CH to isolate host types?
#colnames(ssr_dat_fmt_wcc)




















## OLD
#read in data
ssr_dat_combined <- readRDS("r_data/ssr_dat.rds")
head(ssr_dat_combined)

ssr_dat_fmt <- ssr_dat_combined %>% 
  mutate(pop = word(sample_id,1,-3,"_")) %>% 
  select(id=sample_id,pop,loci,len) %>% 
  pivot_wider(id_cols = c(id,pop), names_from = loci, values_from = len)
head(ssr_dat_fmt)

# remove samples that don't amplify in at least 50% loci??
samples_keep1 <- ssr_dat_combined %>% 
  mutate(pop = word(sample_id,1,-3,"_")) %>% 
  select(sample_id,pop,loci,len) %>% 
  group_by(sample_id) %>% 
  tally() %>% 
  ungroup() 
range(samples_keep1$n) # 1-13
head(samples_keep1)

ssr_dat_combined %>% 
  left_join(samples_keep1, by = "sample_id") %>% 
  mutate(run = word(sample_id,1,1,"_")) %>% 
  group_by(run) %>% 
  summarise(max(n))

samples_keep <- ssr_dat_combined %>% 
  left_join(samples_keep1, by = "sample_id") %>% 
  mutate(run = word(sample_id,1,1,"_")) %>% 
  mutate(keep = case_when(
    run == "201512" & n >= 1 ~ "1",
    run == "201611" & n >= 7 ~ "1",
    run == "201803" & n >= 7 ~ "1",
    run == "201904" & n >= 3.5 ~ "1",
    run == "201911" & n >= 7 ~ "1",
    run == "202103" & n >= 7 ~ "1",
    TRUE ~ "0"
  )) %>% #keep samples that amplify in at least half loci
  pull(sample_id) %>% 
  unique() 
length(samples_keep) #keeping 6953 samples
length(unique(ssr_dat_combined$sample_id)) #out of 6953 total


## LATER: Add back 201512 and 201904 later??
## Wait why would we remove 201904??

# remove loci that don't amplify in at least 20% samples
ssr_loci_discard <- ssr_dat_combined %>% 
  mutate(pop = word(sample_id,1,-3,"_")) %>% 
  select(id=sample_id,pop,loci,len) %>% 
  group_by(loci) %>% 
  tally() %>% 
  ungroup() %>% # max(ssr_loci_discard$n)=6364
  filter(n<1272.8) %>% #<20%
  pull(loci) %>% 
  unique() 
length(ssr_loci_discard) #no loci being removed
length(unique(ssr_dat_combined$loci)) #out of 13 loci originally 

ssr_dat_fmt <- ssr_dat_combined %>% 
  mutate(pop = word(sample_id,1,6,"_")) %>% ##LATER: FIX!!!
  select(id=sample_id,pop,loci,len) %>% 
  filter(id %in% samples_keep & #keep samples work in 50% loci for their run (all)
           !loci %in% ssr_loci_discard & #discard loci that don't work in at least 20% samples (none)
           str_detect(pop,"201512_")==FALSE & #erased the line with 201904_
           str_detect(pop,"2018_LLO")==FALSE &
           str_detect(pop,"2015_TBD_")==FALSE & 
           str_detect(pop,"201803_2017_P_NA_NA_NA")==FALSE &
           str_detect(pop,"202103_2019_W_NA_NA_T1")==FALSE &
           str_detect(pop,"Control")==FALSE &
           str_detect(pop,"Empty")== FALSE
  ) %>% 
  pivot_wider(id_cols = c(id,pop), names_from = loci, values_from = len) %>% 
  mutate(year = word(pop,2,2,"_")) %>%
  mutate(site = word(pop,3,3,"_")) %>%
  mutate(wcc = case_when(
    site %in% c('BH','FR','LL','AR') ~ "CM",
    site == "CR" ~ "CR",
    site == "W" ~ "W"
  )) %>% 
  mutate(pop = paste0(year,"_",wcc,"_",(word(pop,3,-1,"_")))) %>% 
  filter(str_detect(pop,"2020_NA_")==FALSE) %>% 
  select(-year,-site,-wcc)
## LATER: Add back 201512 and 201904 later?? CONTROLS?
## LATER: FIX 201611_2016_ COMMERCIAL SAMPLES NAMING IN DATA CLEANING
#view(ssr_dat_fmt)
head(ssr_dat_fmt)
unique(ssr_dat_fmt$pop)
write_csv(ssr_dat_fmt,"output/all_ssr_dat_fmt.csv")

# Open both in excel and add:
# Top Row: #samples, #loci, #samples, #pops & sample no per pop
# Second Row: corresponding pop names sarting in 4th column
head(ssr_dat_fmt)
length(colnames(ssr_dat_fmt))-2 #13 loci
length(unique(ssr_dat_fmt$id)) #4753 samples
unique(ssr_dat_fmt$pop) #101 populations
ssr_dat_fmt_sample_n <- ssr_dat_fmt %>% 
  group_by(pop) %>% tally() %>% ungroup() %>% select(n,pop)
write_csv(ssr_dat_fmt_sample_n,"output/all_ssr_dat_fmt_sample_n.csv") 
# open in excel, copy, and transpose paste into top 2 rows 
# for samples nos per pop and corresponding pop names

# save with "_genalex" pasted at end of file name








## Remove CR
ssr_dat_fmtnocr <- ssr_dat_fmt %>% 
  mutate(cww = word(pop,2,2,"_")) %>% 
  filter(cww != "CR") %>%
  select(-cww)
head(ssr_dat_fmtnocr)
unique(ssr_dat_fmtnocr$pop) #63
write_csv(ssr_dat_fmtnocr,"output/ssr_dat_fmt_nocr.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmtnocr)
length(colnames(ssr_dat_fmtnocr))-2 #13 loci
length(unique(ssr_dat_fmtnocr$id)) #3227 samples
unique(ssr_dat_fmtnocr$pop) #63 populations
ssr_dat_fmtnocr_sample_n <- ssr_dat_fmtnocr %>% 
  group_by(pop) %>% tally() %>% ungroup() %>% 
  select(n,pop)
write_csv(ssr_dat_fmtnocr_sample_n,"output/ssr_dat_fmtnocr_sample_n.csv") #open in excel, copy, and transpose paste into top 2 rows 
# for samples nos per pop and corresponding pop names

# save with "_genalex" pasted at end of file name






## Only keep 2018-2020
ssr_dat_fmt_1820 <- ssr_dat_fmt %>% 
  mutate(year = word(pop,1,1,"_")) %>%
  filter(year %in% c('2018','2019','2020')) %>% 
  select(-year)
head(ssr_dat_fmt_1820)
#ssr_dat_fmt_1820 %>% group_by(pop) %>% tally() %>% view()
unique(ssr_dat_fmt_1820$pop) #76
write_csv(ssr_dat_fmt_1820,"output/ssr_dat_fmt_1820.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt_1820)
length(colnames(ssr_dat_fmt_1820))-2 #13 loci
length(unique(ssr_dat_fmt_1820$id)) #3625 samples
unique(ssr_dat_fmt_1820$pop) #76 populations
ssr_dat_fmt_1820_sample_n <- ssr_dat_fmt_1820 %>% 
  group_by(pop) %>% tally() %>% ungroup() %>% 
  select(n,pop)
write_csv(ssr_dat_fmt_1820_sample_n,"output/ssr_dat_fmt_1820_sample_n.csv") #open in excel, copy, and transpose paste into top 2 rows 
# for samples nos per pop and corresponding pop names

# save with "_genalex" pasted at end of file name






## Only keep 2018-2020 and remove CR
ssr_dat_fmt_1820nocr <- ssr_dat_fmt_1820 %>% 
  mutate(year = word(pop,1,1,"_")) %>%
  mutate(cww = word(pop,2,2,"_")) %>% 
  filter(year %in% c('2018','2019','2020')) %>%
  filter(cww != "CR") %>%
  select(-year,-cww)
head(ssr_dat_fmt_1820nocr)
#ssr_dat_fmt_1820nocr %>% group_by(pop) %>% tally() %>% view()
unique(ssr_dat_fmt_1820nocr$pop) #50
write_csv(ssr_dat_fmt_1820nocr,"output/ssr_dat_fmt_1820nocr.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt_1820nocr)
length(colnames(ssr_dat_fmt_1820nocr))-2 #13 loci
length(unique(ssr_dat_fmt_1820nocr$id)) #2648 samples
unique(ssr_dat_fmt_1820nocr$pop) #50 populations
ssr_dat_fmt_1820nocr_sample_n <- ssr_dat_fmt_1820nocr %>% 
  group_by(pop) %>% tally() %>% ungroup() %>% 
  select(n,pop)
write_csv(ssr_dat_fmt_1820nocr_sample_n,"output/ssr_dat_fmt_1820nocr_sample_n.csv") #open in excel, copy, and transpose paste into top 2 rows 
# for samples nos per pop and corresponding pop names

# save with "_genalex" pasted at end of file name





## CLONAL SAMPLES ONLY 
str(ssr_dat_combined) #all samples
str(ssr_dat_pt2_final) #clonal samples only
ssr_dat_fmt_clonal <- ssr_dat_pt2_final %>% 
  mutate(pop = word(sample_id,1,-3,"_")) %>% 
  select(id=sample_id,pop,loci,len) %>% 
  pivot_wider(id_cols = c(id,pop), names_from = loci, values_from = len)
head(ssr_dat_fmt_clonal)

# remove samples that don't amplify in at least 50% loci??
samples_keep1 <- ssr_dat_pt2_final %>% 
  mutate(pop = word(sample_id,1,-3,"_")) %>% 
  select(sample_id,pop,loci,len) %>% 
  group_by(sample_id) %>% 
  tally() %>% 
  ungroup() 
range(samples_keep1$n) # 1-5
head(samples_keep1)

ssr_dat_pt2_final %>% 
  left_join(samples_keep1, by = "sample_id") %>% 
  mutate(run = word(sample_id,1,1,"_")) %>% 
  group_by(run) %>% 
  summarise(max(n))

samples_keep <- ssr_dat_pt2_final %>% 
  left_join(samples_keep1, by = "sample_id") %>% 
  mutate(run = word(sample_id,1,1,"_")) %>% 
  mutate(keep = case_when(
    run == "201512" & n >= 0.5 ~ "1",
    run == "201611" & n >= 2 ~ "1",
    run == "201803" & n >= 2.5 ~ "1",
    run == "201904" & n >= 2 ~ "1",
    run == "201911" & n >= 2 ~ "1",
    run == "202103" & n >= 1.5 ~ "1",
    TRUE ~ "0"
  )) %>% #keep samples that amplify in at least half loci
  pull(sample_id) %>% 
  unique() 
length(samples_keep) #keeping 502 samples
length(unique(ssr_dat_pt2_final$sample_id)) #out of 502 total

# remove loci that don't amplify in at least 20% samples
ssr_loci_discard <- ssr_dat_pt2_final %>% 
  mutate(pop = word(sample_id,1,-3,"_")) %>% 
  select(id=sample_id,pop,loci,len) %>% 
  group_by(loci) %>% 
  tally() %>% 
  ungroup() %>% # max(ssr_loci_discard$n)=6364
  filter(n<1272.8) %>% #<20%
  pull(loci) %>% 
  unique() 
length(ssr_loci_discard) #no loci being removed
length(unique(ssr_dat_combined$loci)) #out of 13 loci originally 

ssr_dat_fmt <- ssr_dat_combined %>% 
  mutate(pop = word(sample_id,1,6,"_")) %>% ##LATER: FIX!!!
  select(id=sample_id,pop,loci,len) %>% 
  filter(id %in% samples_keep & #keep samples work in 50% loci for their run (all)
           !loci %in% ssr_loci_discard & #discard loci that don't work in at least 20% samples (none)
           str_detect(pop,"201512_")==FALSE & #erased the line with 201904_
           str_detect(pop,"2018_LLO")==FALSE &
           str_detect(pop,"2015_TBD_")==FALSE & 
           str_detect(pop,"201803_2017_P_NA_NA_NA")==FALSE &
           str_detect(pop,"202103_2019_W_NA_NA_T1")==FALSE &
           str_detect(pop,"Control")==FALSE &
           str_detect(pop,"Empty")== FALSE
  ) %>% 
  pivot_wider(id_cols = c(id,pop), names_from = loci, values_from = len) %>% 
  mutate(year = word(pop,2,2,"_")) %>%
  mutate(site = word(pop,3,3,"_")) %>%
  mutate(wcc = case_when(
    site %in% c('BH','FR','LL','AR') ~ "CM",
    site == "CR" ~ "CR",
    site == "W" ~ "W"
  )) %>% 
  mutate(pop = paste0(year,"_",wcc,"_",(word(pop,3,-1,"_")))) %>% 
  filter(str_detect(pop,"2020_NA_")==FALSE) %>% 
  select(-year,-site,-wcc)
## LATER: Add back 201512 and 201904 later?? CONTROLS?
## LATER: FIX 201611_2016_ COMMERCIAL SAMPLES NAMING IN DATA CLEANING
#view(ssr_dat_fmt)
head(ssr_dat_fmt)
unique(ssr_dat_fmt$pop)
write_csv(ssr_dat_fmt,"output/all_ssr_dat_fmt.csv")




ssr_dat_fmt_clonal2 <- ssr_dat_fmt_clonal %>% 
  mutate(year = word(pop,1,1,"_")) %>%
  mutate(cww = word(pop,2,2,"_")) %>% 
  #filter(year %in% c('2018','2019','2020')) %>%
  #filter(cww != "CR") %>%
  select(-year,-cww)
head(ssr_dat_fmt_clonal2,5)
#ssr_dat_fmt_1820nocr %>% group_by(pop) %>% tally() %>% view()
unique(ssr_dat_fmt_clonal$pop) #50
write_csv(ssr_dat_fmt_1820nocr,"output/ssr_dat_fmt_1820nocr.csv")

# Open both in excel and add first row containing:
# loci no, sample no, pop no,
# then sample no per pop, then second row with corresponding pop names
head(ssr_dat_fmt_1820nocr)
length(colnames(ssr_dat_fmt_1820nocr))-2 #13 loci
length(unique(ssr_dat_fmt_1820nocr$id)) #2648 samples
unique(ssr_dat_fmt_1820nocr$pop) #50 populations
ssr_dat_fmt_1820nocr_sample_n <- ssr_dat_fmt_1820nocr %>% 
  group_by(pop) %>% tally() %>% ungroup() %>% 
  select(n,pop)
write_csv(ssr_dat_fmt_1820nocr_sample_n,"output/ssr_dat_fmt_1820nocr_sample_n.csv") #open in excel, copy, and transpose paste into top 2 rows 
# for samples nos per pop and corresponding pop names

# save with "_genalex" pasted at end of file name

