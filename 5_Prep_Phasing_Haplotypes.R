
## Add Allele Sizes
## SSR Het
## Sample Numbers for Ch.2 
## Prep to Phase Haplotypes


#Load packages
library(tidyverse)
library(readxl)
library(lubridate)


## 2015-2020 & Control dataset: cleaned, filtered, mat, clonality
## Read in or continue
dat_joined_clonal <- readRDS("r_data/filtered_data_mat_clonal.rds")
str(dat_joined_clonal)
length(unique(dat_joined_clonal$loci))

## Calculate SSR allele sizes and add to our dataset 
## LATER: Go through MSAs and fix allele sizes!
allseq <- read_delim("../../../../Desktop/202203_combined_ampseq_database/database/pipeline/output/output_h/HaplotypeAllele.fasta", delim = "\t", col_names = FALSE)
head(allseq)

allseqa <- allseq %>% 
  filter((str_detect(X1,">")==TRUE)) %>% 
  rowid_to_column() %>% 
  mutate(loc1 = word(X1,1,1,"#")) %>% 
  mutate(loc2 = str_sub(loc1,2,-1)) %>% 
  mutate(all = word(X1,2,-1,"#")) %>% 
  select(rowid,loci=loc2,all) %>% 
  distinct()
unique(allseqa$loci)
unique(allseqa$all)
head(allseqa)

allseqb <- allseq %>% 
  filter(str_detect(X1,">")==FALSE) %>% 
  rowid_to_column() %>% 
  select(rowid,seq=X1)
head(allseqb)

allseq2 <- full_join(allseqa,allseqb,by=c('rowid')) %>% 
  mutate(all_len = str_length(seq)) %>% 
  select(-rowid,-seq) %>% 
  distinct()
head(allseq2)
range(allseq2$all_len)
write.csv(allseq2,"output/allele_lengths.csv")
#nrow(allseq2)
#view(allseq2)
#allseq2 %>% select(-rowid) %>% distinct() %>% nrow() #526 rows
#allseq2 %>% select(-rowid,-seq) %>% distinct() %>% nrow() #526 rows

unique(allseq2$loci) #83 loci
unique(dat_joined_clonal$loci) #69 loci (some primers only used in collab samples?)




## Add allele length to our dataset and save
## Rename the dataset 
full_data_clonal <- dat_joined_clonal
str(full_data_clonal)

full_data_clonal1 <- full_data_clonal %>% 
  mutate(all = as.character(all)) %>% 
  left_join(allseq2, by = c('loci','all'))
nrow(full_data_clonal) == nrow(full_data_clonal1)
head(full_data_clonal1)

write_rds(full_data_clonal1,"r_data/filtered_data_mat_clonal_allsize.rds")






## Phasing for SSR Analyses: 
# phase out haplotype IDs using sequencing reads method
full_data_clonal1 <- readRDS("r_data/filtered_data_mat_clonal_allsize.rds")

unique(full_data_clonal1$loci) #69 loci

ssr_loci <- full_data_clonal1 %>% 
  filter(str_detect(loci,"EnMS")|str_detect(loci,"EnCDnew")) %>% 
  pull(loci) %>% 
  unique()
ssr_loci #23 loci

#Determine SSR loci efficiency and remove ones that don't work?
#Note: even samples with het will only be count once!!

# Remove samples that won't be used in Ch. 2
full_data_clonal1 <- full_data_clonal1 %>% 
  filter(!run %in% c("201512","201611-AV9")) #remove runs 201512 & 201611-AV9 (lon664 plate)

# How many samples total?
length(unique(full_data_clonal1$sample_id)) #3171 samples
unique(full_data_clonal1$run)
unique(full_data_clonal1$year)
unique(full_data_clonal1$site)
unique(full_data_clonal1$variety)

# Make table of sample numbers
sample_no <- full_data_clonal1 %>% 
  select(run,year,site,variety,trt,timepoint,sample,sample_id) %>% 
  distinct() %>% 
  group_by(year,run,site,variety,trt,timepoint) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(year,run,site,variety,trt,timepoint)
sample_no
sum(sample_no$n) #3171
write_csv(sample_no,"output/sample_no.csv")

sample_no2 <- full_data_clonal1 %>% 
  select(run,year,site,variety,trt,timepoint,sample,sample_id) %>% 
  distinct() %>% 
  group_by(year,run,site,variety) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(year,run,site,variety)
sample_no2
write_csv(sample_no2,"output/sample_no_simple.csv")

sample_no3 <- full_data_clonal1 %>% 
  select(run,year,site,variety,trt,timepoint,sample,sample_id) %>% 
  distinct()
head(sample_no3)
unique(sample_no3$run)
unique(sample_no3$year)


# make sample number table for table 1 ch. 2 diss
sample_no4 <- sample_no3 %>% 
  filter(year != "Control" & !is.na(variety)) %>% 
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
  group_by(year,site_new,variety_new) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(year,site_new,variety_new)
length(unique(sample_no3$sample_id)) #3170 samples total, including controls
sum(sample_no4$n) #3125 samples total (45 controls)
sample_no5 <- sample_no4 %>%
  pivot_wider(id_cols = c('site_new','variety_new'), names_from = "year", values_from = "n")
sample_no5
write_csv(sample_no5,"output/sample_no_diss_ch2.csv")


## same but detailed
sample_no6 <- sample_no3 %>% 
  filter(year != "Control" & !is.na(variety)) %>% 
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
  group_by(year,site_new,variety_new,trt,timepoint) %>% 
  tally() %>% 
  ungroup() %>% 
  arrange(year,site_new,variety_new,trt,timepoint)
head(sample_no6)
unique(sample_no6$variety_new)
length(unique(sample_no3$sample_id)) #3170 samples total, including controls
sum(sample_no6$n) #3125 samples total (45 controls)
sample_no7 <- sample_no6 %>%
  #mutate_at('timepoint', replace_na, "") %>%  # change timepoint column only: na -> blank
  pivot_wider(id_cols = c('year','site_new','variety_new','trt'), names_from = c("timepoint"), values_from = "n")
sample_no7


sample_no8 <- sample_no6 %>% mutate_at(c('trt','timepoint'),replace_na,"-") %>% arrange(year,site_new,variety_new,timepoint)
sample_no8
sample_no_2015 <- sample_no8 %>% filter(year=="2015") %>% select(-timepoint,-year)
sample_no_2015
sample_no_2016 <- sample_no8 %>% filter(year=="2016") %>% select(-year) %>% select(-timepoint)
sample_no_2016
sample_no_2017 <- sample_no8 %>% filter(year=="2017") %>% select(-year) %>% select(-timepoint)
sample_no_2017
sample_no_2018 <- sample_no8 %>% filter(year=="2018") %>% select(-year) %>% pivot_wider(names_from = "timepoint",values_from = "n")
sample_no_2018
sample_no_2019 <- sample_no8 %>% filter(year=="2019") %>% select(-year) %>% pivot_wider(names_from = "timepoint",values_from = "n")
sample_no_2019
sample_no_2020 <- sample_no8 %>% select(-trt) %>% filter(year=="2020") %>% pivot_wider(names_from = "timepoint",values_from = "n")
sample_no_2020

write_csv(sample_no_2015,"output/sample_no_2015.csv")
write_csv(sample_no_2016,"output/sample_no_2016.csv")
write_csv(sample_no_2017,"output/sample_no_2017.csv")
write_csv(sample_no_2018,"output/sample_no_2018.csv")
write_csv(sample_no_2019,"output/sample_no_2019.csv")
write_csv(sample_no_2020,"output/sample_no_2020.csv")






#create table showing how many samples each SSR amplified in
ssr_effic <- full_data_clonal1 %>% 
  filter(loci %in% ssr_loci) %>% 
  select(loci,sample_id) %>% 
  distinct() %>% 
  group_by(loci) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(perc = (n/length(unique(full_data_clonal1$sample_id)))*100) %>% 
  arrange(perc)
#need to group by loci and everything except sample to count samples its working in
head(ssr_effic)
range(ssr_effic$n)
range(ssr_effic$perc)

write_csv(ssr_effic,"output/ssr_loci_effic.csv")






# only use ssr loci that amplify in at least 25% of all samples 
# Note: (not all loci were used every year)
ssr_loci_keep <- ssr_effic %>% 
  filter(perc >= 25) %>% 
  pull(loci) %>% 
  unique()
ssr_loci #23
ssr_loci_keep #16

#now look at the heterozygosity for these loci across all samples
ssr_het <- full_data_clonal1 %>% 
  filter(loci %in% ssr_loci) %>% 
  select(loci,sample_id,het) %>% 
  distinct() %>% #this removes a TON of rows
  group_by(loci,het) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(het = case_when(
    het == TRUE ~ "het",
    het == FALSE ~ "hom"
  )) %>% 
  pivot_wider(names_from = "het", values_from = "n") %>% 
  mutate(het = as.numeric(het)) %>% 
  mutate(het = case_when(
    is.na(het) ~ 0,
    TRUE  ~ het
  )) %>% 
  mutate(tot = het+hom) %>% 
  mutate(perc_het = (het/tot)*100) %>% 
  arrange(perc_het)
ssr_het
range(ssr_het$perc_het)
view(ssr_het)
write_csv(ssr_het,"output/ssr_het_all_samples.csv")

## Now look at the heterozygosity for these loci again
## but only across high confidence clonal samples
ssr_het2 <- full_data_clonal1 %>% 
  filter(loci %in% ssr_loci & high_conf_clonal == "1") %>% 
  select(loci,sample_id,het) %>% 
  distinct() %>% #this removes a TON of rows
  group_by(loci,het) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(het = case_when(
    het == TRUE ~ "het",
    het == FALSE ~ "hom"
  )) %>% 
  pivot_wider(names_from = "het", values_from = "n") %>% 
  mutate(het = as.numeric(het)) %>% 
  mutate(het = case_when(
    is.na(het) ~ 0,
    TRUE  ~ het
  )) %>% 
  mutate(tot = het+hom) %>% 
  mutate(perc_het = (het/tot)*100) %>% 
  arrange(perc_het)
ssr_het2
range(ssr_het2$perc_het)
view(ssr_het2)
write_csv(ssr_het2,"output/ssr_het_clonal_samples.csv")

#Which SSR loci are het in >50% of all samples?
ssr_loci_het <- ssr_het %>% filter(perc_het>50) %>% pull(loci) %>% unique()
ssr_loci_het #4 loci: EnMS3, EnMS7, EnMS10, EnCDnew_31

#Which SSR loci are het in >50% of conf clonal samples?
ssr_loci_het2 <- ssr_het2 %>% filter(perc_het>50) %>% pull(loci) %>% unique()
ssr_loci_het2 #2 loci: EnMS3, EnMS7


