## Clean data and Create Database
## December 4, 2022


#Load packages
library(tidyverse)
library(readxl)
library(lubridate)

setwd("..")
#setwd("Tutoring/")

# Read in Full AmpSeq data
dat <- read_xlsx("../../../Desktop/202203_combined_ampseq_database/database/pipeline/output/output_h/hap_genotype_final.xlsx")
head(dat)

# Step 1: Clean Data

# Remove rows where haplotype column is NA
dat <- dat %>% filter(!is.na(Haplotypes))
head(dat)

# Get rid of Haplotype column
dat <- dat[,-2]
head(dat)

# Gather column 2 to the end, not the first column "Locus"
# The key is going to be the name of the sample, and the value is going to be the number of reads.
ncol_dat <- ncol(dat)
new_dat <- gather(dat[,2:(ncol(dat))], key = "sample_id", value = "reads")
# Check that dimensions are correct; should equal to dat's (obvervations * (variables - 1)):
dim(new_dat)
nrow(dat) * (ncol(dat) - 1)

# Now we want to bind back the column with the loci. 
# Notice that all the samples with the same name appear consecutively in the gathered table:
head(new_dat, 20)
# This means that we don't need to change the order of the Locus column. 
# We can just repeat it for as many times as the number of samples.
sample_number <- ncol_dat - 1
loci <- rep(dat$Locus, (sample_number))

# Now we can bind this column back to the data:
new_dat <- cbind.data.frame(loci, new_dat)
head(new_dat)
nrow(new_dat)


# Add columns containing the R1 and R2 ID's from the sample file
r1r2_dat <- read_delim("../../../Desktop/202203_combined_ampseq_database/database/pipeline/input/2022_combined_sample_file_final.txt", delim = "\t", col_names = FALSE)
head(r1r2_dat)
r1r2_dat <- r1r2_dat %>%
  select(sample_id_raw = X1, plate_well_raw = X2, r1_raw = X3, r2_raw = X4) %>%
  mutate(sample_id = paste0(sample_id_raw,"__",plate_well_raw))
head(r1r2_dat)
nrow(r1r2_dat) == length(unique(r1r2_dat$sample_id_raw)) #TRUE means no duplicate names in the sample file
nrow(r1r2_dat) == length(unique(new_dat$sample_id)) #TRUE means same number of sample ids in both files

new_dat <- new_dat %>%
  left_join(r1r2_dat)
head(new_dat)

str(new_dat)

test <- new_dat %>% mutate(sample_id=word(sample_id,1,1,"__")) %>% select(sample_id)
test2 <- new_dat %>% select(sample_id_raw) 
all.equal(test,test2)

test <- new_dat %>%
  mutate(sample_id = paste0(word(sample_id,1,1,"__")))
length(unique(new_dat$sample_id)) #5082
length(unique(test$sample_id)) #5082


# Split up the sample ID into run, year, site, variety, trt, timepoint, and sample number
head(new_dat$sample_id)

new_dat2 <- new_dat %>%
  mutate(run = str_extract(sample_id, "20[0-9]{4,}"))
head(new_dat2)

sort(unique(new_dat2$run))
#202203 (newest to oldest)
dat202203 <- new_dat2 %>% 
  filter(run == "202203") %>% 
  mutate(year = word(sample_id, 2, 2, "_")) %>%
  mutate(site = word(sample_id, 3, 3, "_")) %>%
  mutate(variety = word(sample_id, 4, 4, "_")) %>%
  mutate(trt = word(sample_id, 5, 5, "_")) %>%
  mutate(timepoint = word(sample_id, 6, 6, "_")) %>% 
  mutate(sample = word(sample_id, 7, 7, "_")) %>%
  mutate(plate = word(sample_id, 8, 8, "_")) %>%
  mutate(well = word(sample_id, 9, 9, "_")) %>% 
  mutate(year = case_when(
    site == "Empty" ~ "Empty",
    TRUE ~ year
  )) 
head(dat202203)
head(unique(dat202203$sample_id))
unique(dat202203$run)
unique(dat202203$year)
unique(dat202203$site)
unique(dat202203$variety)
dat202203 %>% filter(variety=="NA") %>% pull(sample_id) %>% unique() #HC, Empty, NR, W, R1
dat202203 %>% filter(variety=="NA") %>% pull(site) %>% unique() #later: need to input W sample varieties
unique(dat202203$trt)
unique(dat202203$timepoint)
unique(dat202203$sample)
unique(dat202203$plate)
unique(dat202203$well)
dat202203 %>% filter(sample=="NA") %>% pull(sample_id) %>% unique()
dat202203 %>% filter(sample=="A") %>% pull(sample_id) %>% unique()
dat202203 %>% filter(sample=="B") %>% pull(sample_id) %>% unique()
dat202203 %>% filter(sample=="C") %>% pull(sample_id) %>% unique()


sort(unique(new_dat2$run))
#202103 (newest to oldest)
dat202103 <- new_dat2 %>%
  filter(run == "202103") %>%
  mutate(year = word(sample_id, 2, 2, "_")) %>%
  mutate(site = word(sample_id, 3, 3, "_")) %>%
  mutate(variety = word(sample_id, 4, 4, "_")) %>%
  mutate(trt = word(sample_id, 5, 5, "_")) %>%
  mutate(timepoint = word(sample_id, 6, 6, "_")) %>% 
  mutate(sample = word(sample_id, 7, 7, "_")) %>%
  mutate(plate = word(sample_id, 8, 8, "_")) %>%
  mutate(well = word(sample_id, 9, 9, "_")) %>% 
  mutate(variety = case_when(
    variety == "CN5" ~ "CN",
    TRUE ~ variety
  )) %>% 
  mutate(plate = case_when(
    plate == "12" ~ "Plate12",
    TRUE ~ plate
  ))
str(dat202103)
head(dat202103)
unique(dat202103$run)
unique(dat202103$year)
unique(dat202103$site)
unique(dat202103$variety) #fix BH CN5 to BH CN
dat202103 %>% filter(variety=="NA") %>% pull(site) %>% unique() # W, Empty, Control
#dat202103 %>% filter(variety=="NA" & site == "W") %>% pull(sample_id) %>% unique() #  1 wild vine (T1_78)
unique(dat202103$trt)
unique(dat202103$timepoint)
unique(dat202103$sample)
unique(dat202103$plate)
unique(dat202103$well)
#dat202103 %>% filter(plate=="Plate12") %>% pull(sample_id) %>% unique() #96


sort(unique(new_dat2$run))
#201911
dat201911 <- new_dat2 %>%
  filter(run == "201911") %>%
  mutate(year = word(sample_id, 2, 2, "_")) %>%
  mutate(site = word(sample_id, 3, 3, "_")) %>%
  mutate(variety = word(sample_id, 4, 4, "_")) %>%
  mutate(trt = word(sample_id, 5, 5, "_")) %>%
  mutate(timepoint = word(sample_id, 6, 6, "_")) %>% 
  mutate(sample = word(sample_id, 7, 7, "_")) %>%
  mutate(plate = word(sample_id, 8, 8, "_")) %>%
  mutate(plate = case_when(
    plate == "4" ~ "Plate4",
    TRUE ~ plate
  )) %>% 
  mutate(well = word(sample_id, 9, 9, "_"))
  
str(dat201911)
head(dat201911)
unique(dat201911$run)
unique(dat201911$year)
unique(dat201911$site)
dat201911 %>% filter(site == "TBD") %>% pull(sample_id) %>% unique()
unique(dat201911$variety)
dat201911 %>% filter(variety=="NA") %>% pull(site) %>% unique()
unique(dat201911$trt)
#dat201911 %>% filter(!is.na(trt)) %>% pull(site) %>% unique() #CR only
unique(dat201911$timepoint)
unique(dat201911$sample)
unique(dat201911$plate)
dat201911 %>% filter(plate=="Plate4") %>% pull(sample_id) %>% unique() #96
unique(dat201911$well)



unique(new_dat2$run)
#201904
dat201904 <- new_dat2 %>%
  filter(run == "201904") %>%
  mutate(year = word(sample_id, 2, 2, "_")) %>%
  mutate(site = word(sample_id, 3, 3, "_")) %>%
  mutate(variety = word(sample_id, 4, 4, "_")) %>%
  mutate(trt = word(sample_id, 5, 5, "_")) %>%
  mutate(timepoint = word(sample_id, 6, 6, "_")) %>% 
  mutate(sample = word(sample_id, 7, 7, "_")) %>%
  mutate(plate = word(sample_id, 8, 8, "_")) %>%
  mutate(well = word(sample_id, 9, 9, "_"))

str(dat201904)
unique(dat201904$run)
unique(dat201904$year)
unique(dat201904$site)
unique(dat201904$variety) #Change CN5 to CN
dat201904 %>% filter(variety=="NA") %>% pull(sample_id) %>% unique() # LL, LLO
unique(dat201904$trt)
unique(dat201904$timepoint)
unique(dat201904$sample)
unique(dat201904$plate)
unique(dat201904$well)

unique(new_dat2$run)
#201803
dat201803 <- new_dat2 %>%
  filter(run == "201803") %>%
  mutate(year = word(sample_id, 2, 2, "_")) %>%
  mutate(site = word(sample_id, 3, 3, "_")) %>%
  mutate(variety = word(sample_id, 4, 4, "_")) %>%
  mutate(trt = word(sample_id, 5, 5, "_")) %>%
  mutate(timepoint = word(sample_id, 6, 6, "_")) %>% 
  mutate(sample = word(sample_id, 7, 7, "_")) %>%
  mutate(plate = word(sample_id, 8, 8, "_")) %>%
  mutate(well = word(sample_id, 9, 9, "_")) %>%
  mutate(plate = case_when(
    plate == "1" ~ "Plate1",
    plate == "2" ~ "Plate2",
    plate == "3" ~ "Plate3",
    plate == "4" ~ "Plate4",
    TRUE ~ plate
  )) %>%
  mutate(well = case_when(
    site == "Omer" ~ plate,
    TRUE ~ well
  )) %>%
  mutate(plate = case_when(
    site == "Omer" ~ "Omer",
    TRUE ~ plate
  ))
str(dat201803)
unique(dat201803$run)
unique(dat201803$year)
unique(dat201803$site) #P is??
unique(dat201803$variety) #Change CNB and CNBH to CN
dat201803 %>% filter(variety == "NA") %>% pull(site) %>% unique() #P, Omer, HPM
unique(dat201803$trt)
unique(dat201803$timepoint)
unique(dat201803$sample) #Change BLANK to Empty
unique(dat201803$plate)
dat201803 %>% filter(plate=="Plate1") %>% pull(sample_id) %>% unique() %>% length() #Omer
dat201803 %>% filter(plate=="Plate2") %>% pull(sample_id) %>% unique() %>% length() #Omer, 97
dat201803 %>% filter(plate=="Plate3") %>% pull(sample_id) %>% unique() %>% length() #CR
dat201803 %>% filter(plate=="Plate4") %>% pull(sample_id) %>% unique() %>% length()#CR, 95
dat201803 %>% filter(plate=="Omer") %>% pull(sample_id) %>% unique() %>% length() #42
dat201803 %>% filter(plate=="HPM") %>% pull(sample_id) %>% unique() %>% length() #96
unique(dat201803$well)

sort(unique(new_dat2$run))
#20170221
dat20170221 <- new_dat2 %>%
  filter(run == "20170221") %>%
  mutate(year = word(sample_id, 2, 2, "_")) %>%
  mutate(site = word(sample_id, 3, 3, "_")) %>%
  mutate(variety = word(sample_id, 4, 4, "_")) %>%
  mutate(trt = word(sample_id, 5, 5, "_")) %>%
  mutate(timepoint = word(sample_id, 6, 6, "_")) %>% 
  mutate(sample = word(sample_id, 7, 7, "_")) %>%
  mutate(plate = word(sample_id, 8, 8, "_")) %>%
  mutate(well = word(sample_id, 9, 9, "_")) %>% 
  mutate(plate = case_when(
    plate == "1" ~ "Plate1",
    plate == "4" ~ "Plate4",
    TRUE ~ plate
  ))
str(dat20170221)
unique(dat20170221$run)
unique(dat20170221$year)
unique(dat20170221$site)
unique(dat20170221$variety)
dat20170221 %>% filter(variety == "NA") %>% pull(site) %>% unique()
unique(dat20170221$trt)
unique(dat20170221$timepoint)
unique(dat20170221$sample)
dat20170221 %>% filter(year == "2016") %>% pull(sample) %>% unique()
unique(dat20170221$plate)
sort(unique(dat20170221$well)) #A08 is missing


unique(new_dat2$run)
#201611
dat201611 <- new_dat2 %>%
  filter(run == "201611") %>%
  mutate(year = word(sample_id, 2, 2, "_")) %>%
  mutate(site = word(sample_id, 3, 3, "_")) %>%
  mutate(variety = word(sample_id, 4, 4, "_")) %>%
  mutate(trt = word(sample_id, 5, 5, "_")) %>%
  mutate(timepoint = word(sample_id, 6, 6, "_")) %>% 
  mutate(sample = word(sample_id, 7, 7, "_")) %>%
  mutate(plate = word(sample_id, 8, 8, "_")) %>%
  mutate(well = word(sample_id, 9, 9, "_")) %>% 
  mutate(plate = case_when(
    plate == "1" ~ "Plate1",
    plate %in% c('lon664','664') ~ "PlateX",
    plate == "666" ~ "PlateXX",
    TRUE ~ plate
  ))
str(dat201611)
unique(dat201611$run)
unique(dat201611$year)
unique(dat201611$site)
dat201611 %>% filter(variety=="NA") %>% pull(year) %>% unique()
unique(dat201611$variety)
dat201611 %>% filter(variety=="NA") %>% pull(year) %>% unique()
unique(dat201611$trt)
unique(dat201611$timepoint)
unique(dat201611$sample)
unique(dat201611$plate)
sort(unique(dat201611$well))


unique(new_dat2$run)
#201512
dat201512 <- new_dat2 %>%
  filter(run == "201512") %>%
  mutate(year = word(sample_id, 2, 2, "_")) %>%
  mutate(site = word(sample_id, 3, 3, "_")) %>%
  mutate(variety = word(sample_id, 4, 4, "_")) %>%
  mutate(trt = word(sample_id, 5, 5, "_")) %>%
  mutate(timepoint = word(sample_id, 6, 6, "_")) %>% 
  mutate(sample = word(sample_id, 7, 7, "_")) %>%
  mutate(plate = word(sample_id, 8, 8, "_")) %>%
  mutate(well = word(sample_id, 9, 9, "_")) %>% 
  mutate(plate = case_when(
    plate == "3" ~ "Plate3",
    plate == "4" ~ "Plate4",
    TRUE ~ plate
  ))
str(dat201512)
unique(dat201512$run)
unique(dat201512$year)
unique(dat201512$site)
unique(dat201512$variety)
dat201512 %>% filter(variety == "NA") %>% pull(site) %>% unique()
unique(dat201512$trt)
unique(dat201512$timepoint)
unique(dat201512$sample)
unique(dat201512$plate)
unique(dat201512$well)


# Re-combine data
unique(new_dat2$run)
new_dat3 <- rbind.data.frame(dat202203,dat202103,dat201911,dat201904,dat201803,dat20170221,dat201611,dat201512)
nrow(new_dat3) == nrow(dat202203)+nrow(dat202103)+nrow(dat201911)+nrow(dat201904)+nrow(dat201803)+nrow(dat20170221)+nrow(dat201611)+nrow(dat201512)
nrow(new_dat3) == nrow(new_dat2)

## Check variety listed for all wild samples
new_dat3 %>% 
  filter(site == "W") %>% 
  select(year,site,timepoint,variety,sample) %>% 
  unique() %>% 
  arrange(year,site,timepoint) %>% 
  view()



### Go to 1B_Check_Clean_Data to compare to Ch.1 and look for issues

## Additional Issues: 
# Add varieties for wild 2020

new_dat4 <- new_dat3 %>%
  mutate(run = case_when(
    str_detect(r1_raw,"7776_1570_42889_AVD") ~ paste0(run,"-AVD"),
    str_detect(r1_raw,"7776_1570_42889_AV9") ~ paste0(run,"-AV9"),
    TRUE ~ run
  )) %>% 
  mutate(variety = case_when(
    year == "2020" & site == "W" & sample %in% c(1:72) ~ "VR",
    year == "2020" & site == "W" & sample %in% c(73:117) ~ "VL",
    year == "2018" & site == "LL" & timepoint == "T1" ~ "CH",
    site == "LLO" & sample %in% c(7:10) ~ "AU",
    site == "LLO" & sample %in% c(1:6) ~ "CH",
    year == "2018" & site == "LL" & sample %in% c('1','2','3A','3B') & is.na(variety) ~ "CH",
    variety == "LB" ~ "LM",
    year == "2017" & site == "BH" & variety == "CNBH" & sample %in% c(1:14,'8-2') ~ "CN35",
    year == "2017" & site == "BH" & variety == "CNB" & sample %in% c(1:20) ~ "CN5",
    year == "2017" & site == "BH" & variety == "CNB" & sample %in% c(21:40,'37-2') ~ "CN19",
    year == "2019" & site == "BH" & variety == "CN" ~ "CN5",
    TRUE ~ variety
  )) %>% 
  mutate(plate = case_when(
    plate == "PlateX" ~ "lon664",
    plate == "PlateXX" ~ "lon666",
    run == "20170221" & site == "AR" & plate == "Plate4" ~ "Plate04_lon665",
    plate == "Plate1" ~ "Plate01",
    plate == "Plate2" ~ "Plate02",
    plate == "Plate3" ~ "Plate03",
    plate == "Plate4" ~ "Plate04",
    plate == "lon664" & str_detect(r1_raw,"7776_1570_42889_AV9") ~ "lon664-AV9dups",
    TRUE ~ plate
  )) %>% 
  mutate(plate = case_when(
    r1_raw == "8133_1570_47125_B3DRK_Plate_Mahfee_necator_mahafee_AB1-2016_A08_R1.fastq.gz" ~ "Mahafee",
    r1_raw == "9222_1570_68633_HM5FWBGX5_10392056_LCD_201712_GPM_Plate3_H07_1S-3_R1.fastq.gz" ~ "Plate03",
    r1_raw == "9222_1570_68633_HM5FWBGX5_10392056_LCD_201712_GPM_Plate4_H04_1S-4_R1.fastq.gz" ~ "Plate04",
    TRUE ~ plate
  )) %>% 
  mutate(year = case_when(
    run == "201911" & site == "TBD" & plate == "Plate04" ~ "Omer",
    TRUE ~ year
  )) %>% 
  mutate(well = case_when(
    r1_raw == "8133_1570_47125_B3DRK_Plate_Mahfee_necator_mahafee_AB1-2016_A08_R1.fastq.gz" ~ "A08",
    TRUE ~ well
  )) %>% 
  mutate(site = case_when(
    run == "201911" & site == "TBD" ~ "Omer",
    TRUE ~ site
  )) %>% 
  distinct() %>% 
  replace(.=="NA", NA_character_) %>% 
  mutate(sample_id=paste0(run,"_",year,"_",site,"_",variety,"_",trt,"_",timepoint,"_",sample,"_",plate,"_",well))


## All done?: Check that there's <96 samples per run plate & wells match r1_raw & 5082 sample id's and r1_raws
nrow(new_dat4)==nrow(new_dat3)
length(unique(new_dat4$r1_raw))
length(unique(new_dat4$sample_id))
unique(new_dat4$site)
unique(new_dat4$variety)
test <- new_dat4 %>% select(run,plate,sample_id) %>% distinct() %>% group_by(run,plate) %>% tally() %>% ungroup()
range(test$n) #8-96, all good!

# Save cleaned data
write_rds(new_dat4,"Dissertation_Final_Draft/r_data/cleaned_data.rds")


##Check that it all looks good! :)
new_dat4 %>% filter(year %in% c(2015:2020)) %>% select(site,variety) %>% distinct() %>% arrange(site,variety) 
new_dat4 %>% filter(site=="LLO") %>% select(sample,variety) %>% distinct()
new_dat4 %>% filter(site=="W"&is.na(variety)) %>% pull(sample_id) %>% unique()


## LATER:
# Check details of 2015-2018 collections and varieties (LL CH winery v. house, LLO AU/CH, etc.)
# What is P?? Use 2017 Erysiphe necator sampling (has dates!!) & 201712_AmpSeq_Enecator_Samples
# I thought there was 1 W sample that was NA variety?
# Find and add dates to sample disease summary sheet and/or timepoint?

## LATER LATER:
# Find rest of plate lon665?
# The other plates and samples (Order# file shows Plate3 HPM, and rest of Plate4 has HPM,CPM,SPM, and Plate5 SPM Belachew)
# Find out where 2017 Plates1-3 Bruna's SPM went
# Go through 2020 P3 
# Rename Blanks and Empties to be the same? 

