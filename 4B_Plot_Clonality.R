
## Plot Clonality vs. Incidence, Severity, Date

## 2015-2020 & Control dataset: cleaned, filtered, mat, clonality
## Read in or continue
dat_joined_clonal <- readRDS("r_data/filtered_data_mat_clonal.rds")


## Plot % clonal samples by run, year & run, year, site

## All samples (excpet 201512 run)
str(dat_joined_clonal)
plot_clonal_all <- dat_joined_clonal %>% 
  mutate(perc_clonal_round=case_when(
    perc_clonal <= 0.49999999999 ~ 0,
    perc_clonal >= 0.5 ~ 1
  )) %>% 
  select(run,year,site,variety,timepoint,sample_id,perc_clonal_round) %>% 
  distinct() %>% 
  filter(!is.na(perc_clonal_round) & run!="201512")
head(plot_clonal_all)
unique(plot_clonal_all$perc_clonal_round)

## By Year
plot_clonal_all2 <- plot_clonal_all %>% 
  select(year,sample_id) %>% 
  distinct() %>% 
  group_by(year) %>% 
  tally() %>% 
  ungroup() %>% 
  select(year,tot=n)
plot_clonal_all2  

plot_clonal_all3 <- plot_clonal_all %>% 
  filter(perc_clonal_round=="1") %>% 
  select(year,sample_id) %>% 
  distinct() %>% 
  group_by(year) %>% 
  tally() %>% 
  ungroup() %>% 
  select(year,clonal=n)
plot_clonal_all3  

plot_clonal_all4 <- plot_clonal_all2 %>% 
  left_join(plot_clonal_all3) %>% 
  mutate(perc_clonal=(clonal/tot)*100)
plot_clonal_all4

## Plot barplot of (ALL conf) percent clonality by year
ggplot(plot_clonal_all4, aes(x = year, y = perc_clonal)) +
  geom_bar(stat = 'identity') + ylim(0,100) +
  labs(x="Year",y="Clonal Samples (%)") + theme_bw()
ggsave("output/all_perc_clonal_by_year.tiff", dpi = 300, compression = "lzw")



## By Year and Site

plot_clonal_all2 <- plot_clonal_all %>% 
  select(year,site,sample_id) %>% 
  distinct() %>% 
  group_by(year,site) %>% 
  tally() %>% 
  ungroup() %>% 
  select(year,site,tot=n)
plot_clonal_all2  

plot_clonal_all3 <- plot_clonal_all %>% 
  filter(perc_clonal_round=="1") %>% 
  select(year,site,sample_id) %>% 
  distinct() %>% 
  group_by(year,site) %>% 
  tally() %>% 
  ungroup() %>% 
  select(year,site,clonal=n)
plot_clonal_all3  

plot_clonal_all4 <- plot_clonal_all2 %>% 
  left_join(plot_clonal_all3) %>% 
  mutate(perc_clonal=(clonal/tot)*100)
plot_clonal_all4

## Plot barplot of (ALL conf) percent clonality by year
ggplot(plot_clonal_all4, aes(x = site, y = perc_clonal)) +
  geom_bar(stat = 'identity') + facet_grid(~year, scales = "free", space = "free") + 
  ylim(0,100) + labs(x="Year",y="Clonal Samples (%)") + theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11,angle=45,hjust=1), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/all_perc_clonal_by_year_site.tiff", dpi = 300, compression = "lzw")

## Plot barplot of (ALL conf) percent clonality by year
ggplot(plot_clonal_all4, aes(x = year, y = perc_clonal)) +
  geom_bar(stat = 'identity') + facet_grid(~site, scales = "free", space = "free") + 
  ylim(0,100) + labs(x="Year",y="Clonal Samples (%)") + theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11,angle=90), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/all_perc_clonal_by_site_year.tiff", dpi = 300, compression = "lzw")





## HIGH CONF ONLY
str(dat_joined_clonal)
plot_clonal_conf <- dat_joined_clonal %>% 
  select(run,year,site,variety,timepoint,sample_id,high_conf_clonal) %>% 
  distinct() %>% 
  filter(!is.na(high_conf_clonal) & run!="201512")
head(plot_clonal_conf)
unique(plot_clonal_conf$high_conf_clonal)

## By Year and Site

plot_clonal_conf2 <- plot_clonal_conf %>% 
  select(year,site,sample_id) %>% 
  distinct() %>% 
  group_by(year,site) %>% 
  tally() %>% 
  ungroup() %>% 
  select(year,site,tot=n)
plot_clonal_conf2  

plot_clonal_conf3 <- plot_clonal_conf %>% 
  filter(high_conf_clonal=="1") %>% 
  select(year,site,sample_id) %>% 
  distinct() %>% 
  group_by(year,site) %>% 
  tally() %>% 
  ungroup() %>% 
  select(year,site,clonal=n)
plot_clonal_conf3  

plot_clonal_conf4 <- plot_clonal_conf2 %>% 
  left_join(plot_clonal_conf3) %>% 
  mutate(perc_clonal=(clonal/tot)*100)
plot_clonal_conf4

## Plot barplot of (ALL conf) percent clonality by year
ggplot(plot_clonal_conf4, aes(x = site, y = perc_clonal)) +
  geom_bar(stat = 'identity') + facet_grid(~year, scales = "free", space = "free") + 
  ylim(0,100) + labs(x="Year",y="Clonal Samples (%)") + theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11,angle=45,hjust=1), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/all_perc_clonal_by_year_site.tiff", dpi = 300, compression = "lzw")

## Plot barplot of (ALL conf) percent clonality by year
ggplot(plot_clonal_conf4, aes(x = year, y = perc_clonal)) +
  geom_bar(stat = 'identity') + facet_grid(~site, scales = "free", space = "free") + 
  ylim(0,100) + labs(x="Year",y="Clonal Samples (%)") + theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11,angle=90), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/all_perc_clonal_by_site_year.tiff", dpi = 300, compression = "lzw")






## (Rest is more useful for Ch. 1)

## Plot incidence, severity, day of year, against clonality
## Remove samples: controls or those missing inc, sev, date, or clonality
dat_joined_clonal_simple <- dat_joined_clonal %>% 
  filter(!is.na(avg_inc) & !is.na(avg_sev) & !is.na(date) &
           !is.na(perc_clonal) & year != "Control" & site != "Control") %>% 
  mutate(days=yday(date)) %>% 
  mutate_at(c('avg_inc','avg_sev'), ~ as.numeric(.)) %>% 
  select(sample_id,avg_inc,avg_sev,days,perc_clonal_simple,high_conf_clonal) %>% 
  distinct() %>% 
  mutate(perc_clonal_simple2 = as.factor(perc_clonal_simple)) %>% 
  mutate(high_conf_clonal2 = as.factor(high_conf_clonal))
str(dat_joined_clonal_simple)
nrow(dat_joined_clonal_simple) #1818, not 1138


## Save RDS file
# Removed any samples missing inc, sev, date, or clonality determination
# Removed controls
# (Going to be only 2018 + 2019 since that's when disease assessments were performed)
write_rds(dat_joined_clonal_simple,"r_data/clonality_data__simple.rds")


## Plot histogram of percent clonality by inc, sev, days
ggplot(dat_joined_clonal_simple, aes(x = avg_inc, fill = perc_clonal_simple2)) +
  geom_histogram() + theme_bw() +
  labs(x="Average Incidence (%)",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_inc_vs_clonal.tiff", dpi = 300, compression = "lzw")

ggplot(dat_joined_clonal_simple, aes(x = avg_sev, fill = perc_clonal_simple2)) +
  geom_histogram() + theme_bw() +
  labs(x="Average Severity (%)",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_sev_vs_clonal.tiff", dpi = 300, compression = "lzw")

ggplot(dat_joined_clonal_simple, aes(x = days, fill = perc_clonal_simple2)) +
  geom_histogram() + theme_bw() +
  labs(x="Day of the Year",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_days_vs_clonal.tiff", dpi = 300, compression = "lzw")





## Repeat plots but only for high confidence samples
dat_joined_clonal_simple2 <- dat_joined_clonal_simple %>% 
  filter(!is.na(high_conf_clonal2))

## Plot histogram of high conf percent clonality by inc, sev, days
ggplot(dat_joined_clonal_simple2, aes(x = avg_inc, fill = high_conf_clonal2)) +
  geom_histogram() + theme_bw() +
  labs(x="Average Incidence (%)",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_inc_vs_conf_clonal.tiff", dpi = 300, compression = "lzw")

ggplot(dat_joined_clonal_simple2, aes(x = avg_sev, fill = high_conf_clonal2)) +
  geom_histogram() + theme_bw() +
  labs(x="Average Severity (%)",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_sev_vs_conf_clonal.tiff", dpi = 300, compression = "lzw")

ggplot(dat_joined_clonal_simple2, aes(x = days, fill = high_conf_clonal2)) +
  geom_histogram() + theme_bw() +
  labs(x="Day of the Year",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_days_vs_conf_clonal.tiff", dpi = 300, compression = "lzw")



## Repeat high confidence plots but separate by year
dat_joined_clonal_simple2_yr <- dat_joined_clonal_simple2 %>% 
  mutate(year=word(sample_id,2,2,"_"))
unique(dat_joined_clonal_simple2_yr$year) # 2018 & 2019

## Plot histogram of high conf percent clonality by inc, sev, days
ggplot(dat_joined_clonal_simple2_yr, aes(x = avg_inc, fill = high_conf_clonal2)) +
  geom_histogram() + theme_bw() + facet_wrap(~year) +
  labs(x="Average Incidence (%)",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_inc_vs_conf_clonal_wrapyr.tiff", dpi = 300, compression = "lzw")

ggplot(dat_joined_clonal_simple2_yr, aes(x = avg_sev, fill = high_conf_clonal2)) +
  geom_histogram() + theme_bw() + facet_wrap(~year) +
  labs(x="Average Severity (%)",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_sev_vs_conf_clonal_wrapyr.tiff", dpi = 300, compression = "lzw")

ggplot(dat_joined_clonal_simple2_yr, aes(x = days, fill = high_conf_clonal2)) +
  geom_histogram() + theme_bw() + facet_wrap(~year) +
  labs(x="Day of the Year",y="Count",fill="Clonality (%)") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/obs_days_vs_conf_clonal_wrapyr.tiff", dpi = 300, compression = "lzw")




#####  SKIP  #####
## Which samples are said are aclonal when we expect all samples would be clonal?
str(dat_joined_clonal_simple)
test_clonality <- dat_joined_clonal_simple %>% 
  select(-perc_clonal_simple2) %>% 
  filter((avg_inc < 10 | avg_sev < 5 | days < 195) &
           perc_clonal_simple < .5) %>%
  arrange(avg_inc) %>% 
  view() #15 samples that are said to be aclonal very very very early in the epidemic...
unique(dat_joined_clonal_simple$perc_clonal_simple)
length(unique(dat_joined_clonal_simple$sample_id)) #1138 samples total
##### END SKIP #####








### SKIP: DO LATER
## LATER: Plot percent clonality instead of absolute values (number of samples) by inc, sev, day

## Original/All (Include all conf levels for clonality)
str(dat_joined_clonal_simple)
range(dat_joined_clonal_simple$avg_inc)
range(dat_joined_clonal_simple$avg_sev)
range(dat_joined_clonal_simple$days)
dat_joined_clonal_simple3 <- dat_joined_clonal_simple %>% 
  mutate(avg_inc2 = case_when(
    avg_inc < 10.1 ~ 5,
    avg_inc > 10.1 & avg_inc < 20.1 ~ 15,
    avg_inc > 20.1 & avg_inc < 30.1 ~ 25,
    avg_inc > 30.1 & avg_inc < 40.1 ~ 35,
    avg_inc > 40.1 & avg_inc < 50.1 ~ 45,
    avg_inc > 50.1 & avg_inc < 60.1 ~ 55,
    avg_inc > 60.1 & avg_inc < 70.1 ~ 65,
    avg_inc > 70.1 & avg_inc < 80.1 ~ 75,
    avg_inc > 80.1 & avg_inc < 90.1 ~ 85,
    avg_inc > 90.1 ~ 95
  )) %>% 
  mutate(avg_sev2 = case_when(
    avg_sev < 10.11 ~ 5,
    avg_sev > 10.11 & avg_sev < 20.1 ~ 15,
    avg_sev > 20.1 & avg_sev < 30.1 ~ 25,
    avg_sev > 30.1 & avg_sev < 40.1 ~ 35,
    avg_sev > 40.1 & avg_sev < 50.1 ~ 45,
    avg_sev > 50.1 & avg_sev < 60.1 ~ 55,
    avg_sev > 60.1 & avg_sev < 70.1 ~ 65,
    avg_sev > 70.1 & avg_sev < 80.1 ~ 75,
    avg_sev > 80.1 & avg_sev < 90.1 ~ 85,
    avg_sev > 90.1 ~ 95
  )) %>%  
  mutate(days2 = case_when(
    days < 200.1 ~ 195,
    days > 200.1 & days < 210.1 ~ 205,
    days > 210.1 & days < 220.1 ~ 215,
    days > 220.1 & days < 230.1 ~ 225,
    days > 230.1 & days < 240.1 ~ 235,
    days > 240.1 & days < 250.1 ~ 245,
    days > 250.1 & days < 260.1 ~ 255,
    days > 260.1 & days < 270.1 ~ 265,
    days > 270.1 & days < 280.1 ~ 275,
    days > 280.1 ~ 285
  ))

## Percent - All confidence levels

#inc
all_perc_inc <- dat_joined_clonal_simple3 %>% 
  filter(!is.na(perc_clonal_simple)) %>% 
  group_by(avg_inc2) %>% 
  summarise(mean(perc_clonal_simple)) %>% 
  ungroup()
all_perc_inc  
  
#sev
all_perc_sev <- dat_joined_clonal_simple3 %>% 
  filter(!is.na(perc_clonal_simple)) %>% 
  group_by(avg_sev2) %>% 
  summarise(mean(perc_clonal_simple)) %>% 
  ungroup()
all_perc_sev

#days
all_perc_days <- dat_joined_clonal_simple3 %>% 
  filter(!is.na(perc_clonal_simple)) %>% 
  group_by(days2) %>% 
  summarise(mean(perc_clonal_simple)) %>% 
  ungroup()
all_perc_days


## PLOTS
dat_joined_clonal_simple2 <- dat_joined_clonal_simple %>% 
  filter(!is.na(high_conf_clonal2))

## Plot barplot of (ALL conf) percent clonality by inc, sev, days
ggplot(all_perc_inc, aes(x = avg_inc2, y = `mean(perc_clonal_simple)`)) +
  geom_bar(stat = 'identity') + ylim(0,1) +
  labs(x="Average Incidence (%)",y="Clonal Samples (%)") +
  #theme(text = element_text(size=13), #x and y titles
  #      axis.text.x=element_text(size=11), #x axis markings
  #      axis.text.y=element_text(size=11)) + #y axis markings
  theme_bw()
ggsave("output/inc_vs_perc_clonal.tiff", dpi = 300, compression = "lzw")

ggplot(all_perc_sev, aes(x = avg_sev2, y = `mean(perc_clonal_simple)`)) +
  geom_bar(stat = 'identity') + ylim(0,1) +
  labs(x="Average Severity (%)",y="Clonal Samples (%)") +
  theme_bw()
ggsave("output/sev_vs_perc_clonal.tiff", dpi = 300, compression = "lzw")

ggplot(all_perc_days, aes(x = days2, y = `mean(perc_clonal_simple)`)) +
  geom_bar(stat = 'identity') + ylim(0,1) +
  labs(x="Day of Year",y="Clonal Samples (%)") +
  theme_bw()
ggsave("output/day_vs_perc_clonal.tiff", dpi = 300, compression = "lzw")



## High conf: Include only high conf levels for clonality!! (remove NA's)

#inc
conf_perc_inc <- dat_joined_clonal_simple3 %>% 
  filter(!is.na(high_conf_clonal)) %>% 
  mutate(high_conf_clonal=as.numeric(high_conf_clonal)) %>% 
  group_by(avg_inc2) %>% 
  summarise(mean(high_conf_clonal)) %>% 
  ungroup()
conf_perc_inc  

#sev
conf_perc_sev <- dat_joined_clonal_simple3 %>% 
  filter(!is.na(high_conf_clonal)) %>% 
  mutate(high_conf_clonal=as.numeric(high_conf_clonal)) %>% 
  group_by(avg_sev2) %>% 
  summarise(mean(high_conf_clonal)) %>% 
  ungroup()
conf_perc_sev

#days
conf_perc_days <- dat_joined_clonal_simple3 %>% 
  filter(!is.na(high_conf_clonal)) %>% 
  mutate(high_conf_clonal=as.numeric(high_conf_clonal)) %>% 
  group_by(days2) %>% 
  summarise(mean(high_conf_clonal)) %>% 
  ungroup()
conf_perc_days


## PLOTS
## Plot barplot of high conf percent clonality by inc, sev, days
ggplot(conf_perc_inc, aes(x = avg_inc2, y = `mean(high_conf_clonal)`)) +
  geom_bar(stat = 'identity') + ylim(0,1) +
  labs(x="Average Incidence (%)",y="Clonal Samples (%)") + 
  theme_bw()
ggsave("output/inc_vs_perc_conf_clonal.tiff", dpi = 300, compression = "lzw")

ggplot(conf_perc_sev, aes(x = avg_sev2, y = `mean(high_conf_clonal)`)) +
  geom_bar(stat = 'identity') + ylim(0,1) +
  labs(x="Average Severity (%)",y="Clonal Samples (%)") +
  theme_bw()
ggsave("output/sev_vs_perc_conf_clonal.tiff", dpi = 300, compression = "lzw")

ggplot(conf_perc_days, aes(x = days2, y = `mean(high_conf_clonal)`)) +
  geom_bar(stat = 'identity') + ylim(0,1) +
  labs(x="Day of Year",y="Clonal Samples (%)") +
  theme_bw()
ggsave("output/day_vs_perc_conf_clonal.tiff", dpi = 300, compression = "lzw")









