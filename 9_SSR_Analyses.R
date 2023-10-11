### Script 9: SSR Analyses
## 202309

## Find where we plotted alleles per locus and things like that
## Get nearest and farthest distance algorithms to actually run 

#Load packages
library(tidyverse)
library(readxl)
library(lubridate)

## Load in full data
ssr_dat <- read_csv("output/poppr/data/ssr_dat.csv", col_names = TRUE)
head(ssr_dat)

ssr_dat_long <- ssr_dat %>%
  pivot_longer(!c(sample_id,pop), names_to = "loci", values_to = "all")
head(ssr_dat_long)

write_rds(ssr_dat_long,"output/poppr/full/ssr_dat_long.rds")


## Loci Plots

## number of alleles per locus
allperloc <- ssr_dat_long %>% select(loci,all) %>% distinct() %>% 
  filter(!is.na(all)) %>% group_by(loci) %>% tally() %>% ungroup()
allperloc
range(allperloc$n)
mean(allperloc$n)
write_csv(allperloc,"output/poppr/full/ssr_dat_allperloc.csv")

ssr_dat_long %>% filter(loci == "EnCDnew_33") %>% select(loci,all) %>% distinct()

ggplot(allperloc, aes(x = loci, y = n)) + geom_col() +
  labs(x = "Loci", y = "Number of Alleles") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/loci_numofall.tiff", dpi = 300, compression = "lzw")



allplot <- ssr_dat_long %>% 
  filter(!is.na(all)) %>% 
  mutate(year=word(pop,1,1,"_")) %>% 
  group_by(loci,all,year) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(loci = as.factor(loci)) %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(all = as.factor(all))
head(allplot)

ggplot(allplot, aes(x = loci, y = n, fill = all))+
  facet_grid(~year, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "SSR Loci", fill = "Allele", y = "%") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/ssr_allele_freqs.tiff", dpi = 300, compression = "lzw")


## SINGLE YEAR PLOTS
## Change allele names!

all_len <- ssr_dat_long %>% filter(!is.na(all)) %>% group_by(loci) %>% summarise(min=min(all))
str(all_len)
ssr_dat_long %>% filter(loci == "EnCDnew_32") %>% pull(all)
all_len %>% filter(loci == "EnCDnew_32") 
205-180 #25

set.seed(1)
test1 <- sample(LETTERS)
test1 #Y, D, G

allplot2 <- ssr_dat_long %>% 
  left_join(all_len) %>% 
  filter(!is.na(all)) %>% 
  mutate(all = (all - min)) %>% 
  mutate(all = as.character(all)) %>% 
  mutate(all = case_when(
    all == "0" ~ test1[1],
    all == "1" ~ test1[2],
    all == "3" ~ test1[3],
    all == "4" ~ test1[4],
    all == "5" ~ test1[5],
    all == "6" ~ test1[6],
    all == "8" ~ test1[7],
    all == "9" ~ test1[8],
    all == "10" ~ test1[9],
    all == "11" ~ test1[10],
    all == "12" ~ test1[11],
    all == "14" ~ test1[12],
    all == "15" ~ test1[13],
    all == "16" ~ test1[14],
    all == "17" ~ test1[15],
    all == "19" ~ test1[16],
    all == "20" ~ test1[17],
    all == "21" ~ test1[18],
    all == "24" ~ test1[19],
    all == "25" ~ test1[20],
    all == "28" ~ test1[21],
    all == "29" ~ test1[22],
    all == "30" ~ test1[23],
    all == "33" ~ test1[24],
    all == "48" ~ test1[25],
    all == "7" ~ test1[26],
    TRUE ~ all
  )) %>% 
  mutate(year=word(pop,1,1,"_")) %>% 
  mutate(site=word(pop,3,3,"_")) %>% 
  group_by(year,site,all,loci) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(site = as.factor(site)) %>% 
  mutate(loci = as.factor(loci)) %>% 
  mutate(all = as.factor(all))
head(allplot2)
#range(allplot2$all)
unique(allplot2$all)

# single year - 2015
allplot3 <- allplot2 %>%
  filter(year == "2015") 
head(allplot3)
# sep by year and site 
ggplot(allplot3, aes(x = loci, y = n, fill = all))+
  #facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Loci", fill = "Allele", y = "Alelle Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/2015_allele_freq.tiff", dpi = 300, compression = "lzw")


# single year - 2016
allplot3 <- allplot2 %>%
  filter(year == "2016") 
head(allplot3)
# sep by year and site 
ggplot(allplot3, aes(x = loci, y = n, fill = all))+
  #facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Loci", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/2016_allele_freq.tiff", dpi = 300, compression = "lzw")


# single year - 2017
allplot3 <- allplot2 %>%
  filter(year == "2017") 
head(allplot3)
# sep by year and site 
ggplot(allplot3, aes(x = loci, y = n, fill = all))+
  #facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Loci", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/2017_allele_freq.tiff", dpi = 300, compression = "lzw")


# single year - 2018
allplot3 <- allplot2 %>%
  filter(year == "2018") 
head(allplot3)
# sep by year and site 
ggplot(allplot3, aes(x = loci, y = n, fill = all))+
  #facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Loci", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/2018_allele_freq.tiff", dpi = 300, compression = "lzw")

# single year - 2019
allplot3 <- allplot2 %>%
  filter(year == "2019") 
head(allplot3)
# sep by year and site 
ggplot(allplot3, aes(x = loci, y = n, fill = all))+
  #facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Loci", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/2019_allele_freq.tiff", dpi = 300, compression = "lzw")


# single year - 2020
allplot3 <- allplot2 %>%
  filter(year == "2020") 
head(allplot3)
# sep by year and site 
ggplot(allplot3, aes(x = loci, y = n, fill = all))+
  #facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Loci", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/2020_allele_freq.tiff", dpi = 300, compression = "lzw")





## SINGLE LOCI PLOTS

allplot2 <- ssr_dat_long %>% 
  filter(loci == "EnCDnew_33"&
           !is.na(all)) %>% 
  mutate(year=word(pop,1,1,"_")) %>% 
  mutate(site=word(pop,3,3,"_")) %>% 
  group_by(site,all,year) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(site = as.factor(site)) %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(all = as.factor(all))
head(allplot2)

# single locus
# sep by year and site 
ggplot(allplot2, aes(x = year, y = n, fill = all))+
  facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Site", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/encdnew33_allele_freq_by_year_site.tiff", dpi = 300, compression = "lzw")



allplot2a <- ssr_dat_long %>% 
  filter(loci == "EnCDnew_30"&
           !is.na(all)) %>% 
  mutate(year=word(pop,1,1,"_")) %>% 
  mutate(site=word(pop,3,3,"_")) %>% 
  group_by(site,all,year) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(site = as.factor(site)) %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(all = as.factor(all))
head(allplot2a)

# single locus
# sep by year and site 
ggplot(allplot2a, aes(x = year, y = n, fill = all))+
  facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Site", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/encdnew30_allele_freq_by_year_site.tiff", dpi = 300, compression = "lzw")


allplot2a <- ssr_dat_long %>% 
  filter(loci == "EnCDnew_33"&
           !is.na(all)) %>% 
  mutate(year=word(pop,1,1,"_")) %>% 
  mutate(site=word(pop,3,3,"_")) %>% 
  group_by(site,all,year) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(site = as.factor(site)) %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(all = as.factor(all))
head(allplot2a)

# single locus
# sep by year and site 
ggplot(allplot2a, aes(x = year, y = n, fill = all))+
  facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Site", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/encdnew33_allele_freq_by_year_site.tiff", dpi = 300, compression = "lzw")


allplot2b <- ssr_dat_long %>% 
  filter(loci == "EnCDnew_41"&
           !is.na(all)) %>% 
  mutate(year=word(pop,1,1,"_")) %>% 
  mutate(site=word(pop,3,3,"_")) %>% 
  group_by(site,all,year) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(site = as.factor(site)) %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(all = as.factor(all))
head(allplot2b)

# single locus
# sep by year and site 
ggplot(allplot2b, aes(x = year, y = n, fill = all))+
  facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Site", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/encdnew41_allele_freq_by_year_site.tiff", dpi = 300, compression = "lzw")



allplot2b <- ssr_dat_long %>% 
  filter(loci == "EnCDnew_35"&
           !is.na(all)) %>% 
  mutate(year=word(pop,1,1,"_")) %>% 
  mutate(site=word(pop,3,3,"_")) %>% 
  group_by(site,all,year) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(site = as.factor(site)) %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(all = as.factor(all))
head(allplot2b)

# single locus
# sep by year and site 
ggplot(allplot2b, aes(x = year, y = n, fill = all))+
  facet_grid(~site, scales = "free", space = "free") + 
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Site", fill = "Allele", y = "Allele Frequency") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/encdnew35_allele_freq_by_year_site.tiff", dpi = 300, compression = "lzw")





## make succces/fail rate plot
ssr_dat_long %>% filter(str_detect(sample_id,"-dup")) %>% pull(sample_id) %>% unique() %>% length() #1558

# as is -> treat het as separate samples
ssr_loc <- ssr_dat_long %>% 
  mutate(all_dat = case_when(
    !is.na(all) ~ "Y",
    is.na(all) ~ "N"
  )) %>% 
  select(loci,sample_id,all_dat) %>% 
  distinct() %>% #this line removes samples that are het for that allele
  group_by(loci,all_dat) %>% 
  tally() %>% 
  pivot_wider(names_from = all_dat, values_from = n) %>% 
  mutate(sum=Y+N,
         perc_pass = Y/sum,
         perc_fail = N/sum)
head(ssr_loc)

range(ssr_loc$perc_pass) #0.39 - 0.94
ssr_loc %>% filter(perc_pass>0.9) %>% pull(loci) #EnCDnew_30 + EnCDnew_33
write_csv(ssr_loc, "output/poppr/full/ssr_dat_loci_perc_pass_fail.csv")

# success rate

ggplot(ssr_loc, aes(x = loci, y = perc_pass)) + geom_col() +
  labs(x = "Loci", y = "Success Rate (%)") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/loci_success_rate.tiff", dpi = 300, compression = "lzw")

# failure rate

ggplot(ssr_loc, aes(x = loci, y = perc_fail)) + geom_col() +
  labs(x = "Loci", y = "Failure Rate (%)") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/loci_failure_rate.tiff", dpi = 300, compression = "lzw")




ssr_dat2 <- ssr_dat_long %>% 
  mutate(sample_id = case_when(
    str_detect(sample_id,"-dup") ~ word(sample_id,1,-2,"-"),
    TRUE ~ sample_id
  ))
ssr_dat2 %>% filter(str_detect(sample_id,"-dup"))
ssr_dat %>% filter(str_detect(sample_id,"-dup")&!is.na(EnCDnew_33)) %>% select(sample_id,EnCDnew_33)
#201611-AVD_2016_AR_CH_NA_NA_1-19_lon666_A08-dup
ssr_dat %>% filter(str_detect(sample_id,"201611-AVD_2016_AR_CH_NA_NA_1-19_lon666_A08")) %>% select(sample_id,EnCDnew_33)
ssr_dat2 %>% filter(str_detect(sample_id,"201611-AVD_2016_AR_CH_NA_NA_1-19_lon666_A08")&loci=="EnCDnew_33")

head(ssr_dat2)
ssr_loc2 <- ssr_dat2 %>% 
  mutate(all_dat = case_when(
    !is.na(all) ~ "Y",
    is.na(all) ~ "N"
  )) %>% 
  select(loci,sample_id,all_dat) %>% 
  distinct() %>% #this line removes samples that are het for that allele
  group_by(loci,all_dat) %>% 
  tally() %>% 
  pivot_wider(names_from = all_dat, values_from = n) %>% 
  mutate(sum=Y+N,
         perc_pass = Y/sum,
         perc_fail = N/sum)
head(ssr_loc2)

range(ssr_loc2$perc_pass) #0.37 - 0.92
ssr_loc2 %>% filter(perc_pass>0.9) %>% pull(loci) #EnCDnew_30 + EnCDnew_33
write_csv(ssr_loc2, "output/poppr/full/ssr_dat_loci_perc_pass_fail_remdupAKAhet.csv")


# success rate

ggplot(ssr_loc2, aes(x = loci, y = perc_pass)) + geom_col() +
  labs(x = "Loci", y = "Success Rate (%)") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/loci_success_rate_remdupAKAhet.tiff", dpi = 300, compression = "lzw")

# failure rate

ggplot(ssr_loc2, aes(x = loci, y = perc_fail)) + geom_col() +
  labs(x = "Loci", y = "Failure Rate (%)") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/loci_failure_rate_remdupAKAhet.tiff", dpi = 300, compression = "lzw")





## look at ssr loci het (treated as separate samples in SSR dataset for PoppR)
ssr_dat3 <- ssr_dat2 %>% 
  distinct() %>% 
  group_by(sample_id,loci) %>% 
  tally() %>% 
  ungroup() %>% #unique(ssr_dat3$n)
  mutate(het=case_when(
    n==1 ~ FALSE,
    n==2 ~ TRUE
  )) %>% 
  select(-n)
head(ssr_dat3)

ssr_het1 <- ssr_dat3 %>% filter(het==TRUE) %>% group_by(loci) %>% tally() %>% ungroup() %>% select(loci,het=n)
ssr_het2 <- ssr_dat3 %>% filter(het==FALSE) %>% group_by(loci) %>% tally() %>% ungroup() %>% select(loci,hom=n)
ssr_het <- full_join(ssr_het1,ssr_het2, by = c('loci')) %>% 
  mutate(sum=het+hom,
         perc_het=het/sum)
head(ssr_het)
range(ssr_het$het) #range of 143 to 1118 samples het for a given locus
unique(ssr_het$sum) #3058 samples
signif(range(ssr_het$perc_het), digits = 2) # 4.7 - 37 % het

write_csv(ssr_het,"output/poppr/full/ssr_het.csv")

ggplot(ssr_het, aes(x = loci, y = perc_het)) + geom_col() +
  labs(x = "Loci", y = "Heterozygosity Rate (%)") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/loci_het_rate.tiff", dpi = 300, compression = "lzw")


# SKIP: clonality: probably not important, right?










## Missingness Heatmap (keep dups as is as sep samples)

# Read in data or skip
#ssr_dat_long <- readRDS("output/poppr/full/ssr_dat_long.rds")

#1: group samples by year only
ssr_missing <- ssr_dat_long %>%
  mutate(all_dat = case_when(
    !is.na(all) ~ "Y",
    is.na(all) ~ "N"
  )) %>% 
  mutate(grp=word(pop,1,1,"_")) %>% 
  select(grp,loci,sample_id,all_dat) %>% 
  distinct() %>% #length(unique(ssr_missing$sample_id)) 4616 before and after
  group_by(loci,all_dat,grp) %>% 
  tally() %>% 
  pivot_wider(names_from = all_dat, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(sum=Y+N,
         perc_pass = round(Y/sum,digits=2),
         perc_fail = round(N/sum,digits=2))
head(ssr_missing)

range(ssr_missing$perc_fail)

ggplot(ssr_missing, aes(grp, loci, fill = perc_fail)) + geom_tile(color = "black") + 
  labs(x = "Year", y = "Locus", fill = "Failure Rate") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/heatmap_perc_fail_by_yr.tiff", dpi = 300, compression = "lzw")

ggplot(ssr_missing, aes(grp, loci, fill = perc_pass)) + geom_tile(color = "black") + 
  labs(x = "Year", y = "Locus", fill = "Success Rate") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/heatmap_perc_success_by_yr.tiff", dpi = 300, compression = "lzw")


# Heatmap 2: group samples by year and CTW
ssr_missing2 <- ssr_dat_long %>% 
  mutate(all_dat = case_when(
    !is.na(all) ~ "Y",
    is.na(all) ~ "N"
  )) %>% 
  mutate(grp=word(pop,1,2,"_")) %>% 
  select(grp,loci,sample_id,all_dat) %>% 
  distinct() %>%  #6747 before, 4888 after (worked correctly)
  group_by(loci,all_dat,grp) %>% 
  tally() %>% 
  pivot_wider(names_from = all_dat, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(sum=Y+N,
         perc_pass = Y/sum,
         perc_fail = N/sum)
head(ssr_missing2)

ggplot(ssr_missing2, aes(grp, loci, fill = perc_fail)) + geom_tile(color = "black") + 
  labs(x = "Population", y = "Locus", fill = "Failure Rate") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/heatmap_perc_fail_by_yr_ctw.tiff", dpi = 300, compression = "lzw")

ggplot(ssr_missing2, aes(grp, loci, fill = perc_pass)) + geom_tile(color = "black") + 
  labs(x = "Population", y = "Locus", fill = "Success Rate") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/heatmap_perc_success_by_yr_ctw.tiff", dpi = 300, compression = "lzw")



# Heatmap 3: group samples by year and CTW and site
ssr_missing3 <- ssr_dat_long %>% 
  mutate(all_dat = case_when(
    !is.na(all) ~ "Y",
    is.na(all) ~ "N"
  )) %>% 
  mutate(grp=word(pop,1,3,"_")) %>% 
  select(grp,loci,sample_id,all_dat) %>% 
  distinct() %>%  #6747 before, 4888 after (worked correctly)
  group_by(loci,all_dat,grp) %>% 
  tally() %>% 
  pivot_wider(names_from = all_dat, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(sum=Y+N,
         perc_pass = Y/sum,
         perc_fail = N/sum)
head(ssr_missing3)

ggplot(ssr_missing3, aes(grp, loci, fill = perc_fail)) + geom_tile(color = "black") + 
  labs(x = "Population", y = "Locus", fill = "Failure Rate") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/heatmap_perc_fail_by_yr_ctw_site.tiff", dpi = 300, compression = "lzw")

ggplot(ssr_missing3, aes(grp, loci, fill = perc_pass)) + geom_tile(color = "black") + 
  labs(x = "Population", y = "Locus", fill = "Success Rate") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/heatmap_perc_success_by_yr_ctw_site.tiff", dpi = 300, compression = "lzw")





