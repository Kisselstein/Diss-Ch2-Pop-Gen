
## PoppR Analyses

#Load packages
library(tidyverse)
#library(readxl)
#library(lubridate)
library(poppr)
library(vegan)
#library(treemap)
#library(vegan)
#library(ape)

#sort data by pop column?

# Load data
ssr1 <- read.genalex("output/poppr/data/ssr_dat_genalex.csv",ploidy=1)
splitStrata(ssr1) <- ~year/ctw/site/variety/trt/timepoint
setPop(ssr1) <- ~year/ctw/site/variety/trt/timepoint
#popNames(ssr1) #104

ssr1 #3706 mlg, 4616 ind 13 loci

#tab(ssr1)
nmll(ssr1) #3706
length(mll(ssr1)) #4616 MLLs
#mlg(ssr1) #3706 mlgs, 4616 ind
length(levels(ssr1$pop)) #104 pops
levels(ssr1$loc.fac) #13 loci (note order)
#[1] "EnMS2"      "EnCDnew_24" "EnCDnew_27" "EnCDnew_29" "EnCDnew_30" "EnCDnew_33"
#[7] "EnCDnew_35" "EnCDnew_39" "EnCDnew_41" "EnMS8"      "EnCDnew_25" "EnCDnew_32"
#[13] "EnMS6" 

locus_replengths <- c(3,4,4,4,4,4,4,4,4.5,3,4,4,3) #repeat motif lengths for all 13 loci, in order
ssr1reps <- fix_replen(ssr1, c(locus_replengths)) # #bp for each locus
#Warning: The repeat lengths for EnCDnew_29, EnCDnew_41 are not consistent.
#This might be due to inconsistent allele calls or repeat lengths that are too large.
#Check the alleles to make sure there are no duplicated or similar alleles that might end up being the same after division.
#Repeat lengths with some modification are being returned: EnCDnew_24


ssr1 #3706 mlg, 4616 ind, 13 loci
test <- missingno(ssr1,type="geno",cutoff=0.6) #removed a lot of 2018 samples
test #3149 mlg, 3641 ind
test2 <- missingno(test,type="loci",cutoff=0.5)
test2 #12 loci: removed EnCDnew_35
#ssr1_filtered <- filter_stats(test2, distance = bruvo.dist, replen = ssr1reps, 
#                              plot = TRUE, missing = "ignore")

# Use this function to find cutoff if you don't see a small initial peak!
#print(farthest_thresh <- cutoff_predictor(ssr1_filtered$farthest$THRESHOLDS)) #0.02083
#print(average_thresh  <- cutoff_predictor(ssr1_filtered$average$THRESHOLDS)) #0.0006076
average_thresh  <- 0.0006076
#print(nearest_thresh  <- cutoff_predictor(ssr1_filtered$nearest$THRESHOLDS)) #0

# Now define MLL for ssr1 with the following criteria:
#[t]	threshold	0.0006076 (chose average)
#[d]	distance	Bruvo’s Distance
#[a]	algorithm	- average (not Farthest neighbor)

#mlg.filter(ssr1, missing = "ignore", distance = bruvo.dist, replen = ssr1reps, algorithm = "f") <- farthest_thresh
value <- mlg.filter(test2, missing = "ignore", distance = bruvo.dist, replen = ssr1reps, algorithm = "f") <- average_thresh
value #0.0006076
#Evaluate this output very critically

#mlg.filter(test2, distance = bruvo.dist, replen = ssr1reps, algorithm = "f") <- farthest_thresh
#test2 #takes too long to run



## Genotype curve, pairwise IA, missingness, etc.

test2

## Genotype Accumulation Curve
#assess power to discriminate between ind. given a random sample of n loci
tiff("output/poppr/full/genotype_curve.tiff", width = 300, height = 200)
gac <- genotype_curve(ssr1, sample = 1000, quiet = TRUE)
p <- last_plot()
p + 
  labs(x = "Number of Loci", y = "MLGs") + 
  theme(text = element_text(size=14), #x and y titles
        axis.text.x=element_text(size=12), #x axis markings
        axis.text.y=element_text(size=12)) #y axis markings

ggsave("output/poppr/full/genotype_curve.tiff", dpi = 300, compression = "lzw")

dev.off()
#for each boxplot, n loci were randomly sampled 1000 times in order to create the distribution
#I think a plateau is beginning at 10 loci

## LATER: Make boxplot manually to increase axis title sizes and standardize plot sizes?


# Pairwise r¯d over all loci (Index of Association)
#resample.ia(ssr1,reps=100,hist=TRUE)
pair_data <- pair.ia(test2)
pair_data <- as.data.frame(pair_data)
head(pair_data)
write_csv(pair_data,"output/poppr/full/pairwise_rbard_ia.csv")

pair_plot <- pair_data %>% rownames_to_column() %>% 
  mutate(locus1=word(rowname,1,1,":"),
         locus2=word(rowname,2,2,":")) %>% 
  select(-rowname,-Ia) %>% 
  mutate(locus1 = factor(pair_plot$locus1,levels = unique(pair_plot$locus1),ordered = T)) %>% 
  mutate(locus2 = factor(pair_plot$locus2,levels = unique(pair_plot$locus2),ordered = T))
head(pair_plot)
#view(pair_plot)

ggplot(pair_plot, aes(locus2, locus1, fill = rbarD)) + geom_tile(color = "black") + 
  labs(x = "Locus", y = "Locus") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/pairwise_rbard_ia_heatmap.tiff", dpi = 300, compression = "lzw")





# Missing data in the population heat map and table
setPop(test2) <- ~year/ctw
test2 #14 pops

#plot missingness heatmap - not colorblind friendly at all!
#tiff("output/poppr/full/missing_heatmap_yc.tiff", width = 800, height = 800)
info_table(test2, type = "missing", plot = TRUE)
#dev.off() 

p <- last_plot()
p + 
  #labs(y = "value") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings

#ggsave("output/poppr/full/heatmap_perc_fail_by_yr_ctw_CB1.tiff", dpi = 300, compression = "lzw")

# Make and save missingness table
ssr_info_table <- as.data.frame(info_table(test2, type = "missing", plot = FALSE)) %>% 
  rownames_to_column() %>% arrange(rowname)
head(ssr_info_table)
write_csv(ssr_info_table,"output/poppr/full/missing_yctw.csv")

# Missingness Heatmap - manual
ssr_info_table_long <- ssr_info_table %>% mutate(E=Mean, Population=word(rowname,1,2,"_")) %>% select(-Mean,-rowname) %>% 
  pivot_longer(cols = starts_with("E"), names_to = "Locus", values_to = "value") %>% 
  mutate(Locus = case_when(
    Locus=="E" ~ "Mean",
    TRUE ~ Locus
  )) %>%
  mutate(Population = case_when(
    !is.na(Population) ~ Population,
    TRUE ~ "Total"
  ))
head(ssr_info_table_long)

ggplot(ssr_info_table_long, aes(Population, Locus, fill = value)) + geom_tile(color = "black") + 
  labs(fill = "Failure Rate") + 
  scale_fill_viridis_c() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/full/heatmap_perc_fail_by_yr_ctw.tiff", dpi = 300, compression = "lzw")






## DIVERSITY: CTW only

# ctw only - 3 pops
setPop(test2) <- ~ctw
popNames(test2)
test2tab <- mlg.table(test2, quiet = TRUE) 
test2tab[, 1:5] # Showing the first 5 columns and 3 rows of the table.
str(test2tab)
test2tab_df <- as.data.frame(test2tab) %>% rownames_to_column() %>% arrange(rowname) %>% column_to_rownames() %>% as.vector()
str(test2tab_df)
test2tab_df[, 1:5]

test2_table <- poppr(test2) 
test2_table2 <- test2_table %>% as.data.frame() %>% 
  arrange(Pop)
test2_table2
write_csv(test2_table2,"output/poppr/full/ctw_poppr.csv")
  

# CF and other diversity functions
# CF: clonal fraction statistic (# MLGs/N)
myCF <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){ # if it's a matrix
    res <- rowSums(x > 0)/rowSums(x)
  } else {                 # if it's a vector
    res <- sum(x > 0)/sum(x)
  }
  return(res)
}

test2stat2 <- diversity_stats(test2tab_df, CF = myCF)
test2stat2
test2stat2 <- as.data.frame(test2stat2)
write_csv(test2stat2,"output/poppr/full/ctw_div_stats_cf.csv")
# You can use filtered or custom MLGs to compare diversity


#rarefaction:eMLGs
rarecurve(test2tab_df, ylab="Number of expected MLGs", sample=min(rowSums(test2tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")
min(rowSums(test2tab_df)) #356
rowSums(test2tab_df)

# Div CI
test2_divc <- diversity_ci(test2tab_df, n = 100L, raw = FALSE)
test2_divc
write_csv(test2_divc,"output/poppr/full/ctw_div_ci.csv")
# boxplots biased; read doc at ?diversity_ci
# ci can exist outside of possible range, like BB_10 and BB_11.

## Jack-knife rarefaction:
# diversity_ci() and diversity_boot() have the option
# meaning your data will be randomly sub-sampled to either 
# the smallest population size, or your specified n.rare value, 
# whichever is bigger. Here’s an example with the previous data set:

#DIV CI RARE (minimum)
#table
test2rare <- diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)
test2rare # ctw, n=356
write_csv(test2rare,"output/poppr/full/ctw_div_ci_rare356.csv")

#plot
diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)
ggsave("output/poppr/full/ctw_div_ci_rare356.tiff", dpi = 300, compression = "lzw")








## DIVERSITY: Site 
# Q. Should I keep CTW to make it obvious for readers or put into figure description?

# site - 13 pops
setPop(test2) <- ~site
popNames(test2)
test2
test2tab <- mlg.table(test2, sublist = c(-12), quiet = TRUE) 
test2tab[, 1:5] # Showing the first 5 columns and all rows of the table.
#str(test2tab)
test2tab_df <- as.data.frame(test2tab) %>% rownames_to_column() %>% arrange(rowname) %>% column_to_rownames() %>% as.vector()
rownames(test2tab_df)

test2_table <- poppr(test2) 
test2_table2 <- test2_table %>% as.data.frame() %>% 
  arrange(Pop)
test2_table2
write_csv(test2_table2,"output/poppr/full/ctw_site_poppr.csv")

# Use CF diversity function
test2stat2 <- diversity_stats(test2tab_df, CF = myCF)
test2stat2
test2stat2 <- as.data.frame(test2stat2) %>% rownames_to_column()
write_csv(test2stat2,"output/poppr/full/ctw_site_div_stats_cf.csv")
# You can use filtered or custom MLGs to compare diversity

#rarefaction:eMLGs
rarecurve(test2tab_df, ylab="Number of expected MLGs", sample=min(rowSums(test2tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")
min(rowSums(test2tab_df)) #15
rowSums(test2tab_df)

# Div CI
test2_divc <- diversity_ci(test2tab_df, n = 100L, raw = FALSE)
test2_divc
write_csv(test2_divc,"output/poppr/full/ctw_site_div_ci.csv")
# boxplots biased; read doc at ?diversity_ci
# ci can exist outside of possible range, like BB_10 and BB_11.

## Jack-knife rarefaction:
# diversity_ci() and diversity_boot() have the option
# meaning your data will be randomly sub-sampled to either 
# the smallest population size, or your specified n.rare value, 
# whichever is bigger. Here’s an example with the previous data set:

#DIV CI RARE (minimum)
#table
test2rare <- diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)
test2rare # year, n=15
write_csv(test2rare,"output/poppr/full/ctw_site_div_ci_rare15.csv")
#plot
diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)
p <- last_plot()
p + 
  labs(y = "value") + 
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings

ggsave("output/poppr/full/ctw_site_div_ci_rare15.tiff", dpi = 300, compression = "lzw")







## DIVERSITY: Year & Trts 
# Note: Only want Trial (T) populations

# Yr & Trts - 29 pops
# Remove pops 1, 2, 13, 14, 21, 26 from list
setPop(test2) <- ~year/trt
popNames(test2)
test2tab <- mlg.table(test2, sublist = c(-1,-2,-10,-13,-14,-21,-26), quiet = TRUE) 
test2tab[, 1:5] # Showing the first 5 columns and all rows of the table.
#str(test2tab)

test2_table <- poppr(test2) 
test2_table2 <- test2_table %>% as.data.frame() %>% 
  mutate(year = word(Pop,1,1,"_")) %>% 
  mutate(trt = word(Pop,2,2,"_")) %>% #sort(unique(test2tab_df$trt))
  mutate(Pop_New=case_when(
    trt == "1" ~ paste0(year,"_3"),
    trt == "10" ~ paste0(year, "_3,7"),
    year == "2018" & trt == "11" ~ paste0(year, "_7,11"),
    year == "2019" & trt == "11" ~ paste0(year, "_3,7*"), # add * bc second one for 2019
    #trt == "11S" ~ paste0(year, "_3,7"),
    trt == "13" ~ paste0(year, "_7"),
    trt == "15" ~ paste0(year, "_3*"), # add * bc second one for 2018
    trt == "16Q" ~ paste0(year, "_11"),
    trt == "18" ~ paste0(year, "_UTC"),
    trt == "1D" ~ paste0(year, "_3"),
    trt == "2" ~ paste0(year, "_3"),
    year == "2016" & trt == "29" ~ paste0(year, "_UTC"),
    year == "2018" & trt == "29" ~ paste0(year, "_UTC"), #actually JMS stylet oil 
    trt == "2Q" ~ paste0(year, "_13"),
    trt == "3" ~ paste0(year, "_3"),
    trt == "5" ~ paste0(year, "_3,40"),
    trt == "Control" ~ paste0(year, "_UTC"),
    trt == "D" ~ paste0(year, "_3"),
    trt == "E" ~ paste0(year, "_3*"), # add * bc second one for 2019
    trt == "Q" ~ paste0(year, "_3,11"),
    trt == "S" ~ paste0(year, "_3,7"),
    TRUE ~ Pop
  )) %>% 
  arrange(Pop_New)
test2_table2
write_csv(test2_table2,"output/poppr/full/T_year_trt_poppr.csv")


test2tab_df <- as.data.frame(test2tab) %>% rownames_to_column() %>%
  arrange(rowname) %>% column_to_rownames() %>% as.vector()
rownames(test2tab_df)
test2tab_df[,1:5]

test2tab_rownames <- as.data.frame(test2tab_df) %>% rownames_to_column() %>% 
  mutate(year = word(rowname,1,1,"_")) %>% 
  mutate(trt = word(rowname,2,2,"_")) %>% #sort(unique(test2tab_df$trt))
  mutate(rowname=case_when(
    trt == "1" ~ paste0(year,"_3"),
    trt == "10" ~ paste0(year, "_3,7"),
    year == "2018" & trt == "11" ~ paste0(year, "_7,11"),
    year == "2019" & trt == "11" ~ paste0(year, "_3,7*"), # add * bc second one for 2019
    #trt == "11S" ~ paste0(year, "_3,7"),
    trt == "13" ~ paste0(year, "_7"),
    trt == "15" ~ paste0(year, "_3*"), # add * bc second one for 2018
    trt == "16Q" ~ paste0(year, "_11"),
    trt == "18" ~ paste0(year, "_UTC"),
    trt == "1D" ~ paste0(year, "_3"),
    trt == "2" ~ paste0(year, "_3"),
    year == "2016" & trt == "29" ~ paste0(year, "_UTC"),
    year == "2018" & trt == "29" ~ paste0(year, "_UTC"), #actually JMS stylet oil 
    trt == "2Q" ~ paste0(year, "_13"),
    trt == "3" ~ paste0(year, "_3"),
    trt == "5" ~ paste0(year, "_3,40"),
    trt == "Control" ~ paste0(year, "_UTC"),
    trt == "D" ~ paste0(year, "_3"),
    trt == "E" ~ paste0(year, "_3*"), # add * bc second one for 2019
    trt == "Q" ~ paste0(year, "_3,11"),
    trt == "S" ~ paste0(year, "_3,7"),
    TRUE ~ rowname
  )) %>% pull(rowname)
test2tab_rownames
length(test2tab_rownames)==length(unique(test2tab_rownames))
length(test2tab_rownames)==length(rownames(test2tab_df))
rownames(test2tab_df) <- test2tab_rownames
test2tab_df[,1:5]


#now finish the sort by rowname
#rownames(test2tab_df)
nrow(test2tab_df)
test2tab_df <- test2tab_df[c(2,1,3,4,5,7,6,#2015-2016
                             9,10,12,11,8,#2017
                             13,16,14,15,17,#2018
                             22,20,18,19,21),] #2019
nrow(test2tab_df)
test2tab_df[,1:5]


# Use CF diversity function
test2stat2 <- diversity_stats(test2tab_df, CF = myCF)
test2stat2 <- as.data.frame(test2stat2) %>% rownames_to_column()
test2stat2
write_csv(test2stat2,"output/poppr/full/T_year_trt_div_stats_cf.csv")
# You can use filtered or custom MLGs to compare diversity

#rarefaction:eMLGs
rarecurve(test2tab_df, ylab="Number of expected MLGs", sample=min(rowSums(test2tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")
min(rowSums(test2tab)) #19
rowSums(test2tab)

# Div CI
test2_divc <- diversity_ci(test2tab_df, n = 100L, raw = FALSE)
test2_divc
write_csv(test2_divc,"output/poppr/full/T_year_trt_div_ci.csv")
# boxplots biased; read doc at ?diversity_ci
# ci can exist outside of possible range, like BB_10 and BB_11.

## Jack-knife rarefaction:
# diversity_ci() and diversity_boot() have the option
# meaning your data will be randomly sub-sampled to either 
# the smallest population size, or your specified n.rare value, 
# whichever is bigger. Here’s an example with the previous data set:

#DIV CI RARE (minimum)
#table
test2rare <- diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)
test2rare # year, n=19
write_csv(test2rare,"output/poppr/full/T_year_trt_div_ci_rare19.csv")
#plot
diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)

p <- last_plot()
p + 
  labs(y = "value") + 
  theme(text = element_text(size=14), #x and y titles
        axis.text.x=element_text(size=12), #x axis markings
        axis.text.y=element_text(size=12)) #y axis markings

ggsave("output/poppr/full/T_year_trt_div_ci_rare19.tiff", dpi = 300, compression = "lzw")

















## Create file and simplified table containing CR trt information
trial_trts <- readRDS("r_data/database.rds") %>% 
  filter(!is.na(trt)&year!="2020"&site=="CR") %>% select(year,trt,trt_desc) %>% 
  distinct() %>% arrange(year) %>% 
  mutate(fung_name = case_when(
    trt == "1" ~ "Rhyme",
    trt == "10" ~ "Luna Experience", #same for 2018 and 2019
    year == "2018" & trt == "11" ~ "Luna Sensation",
    year == "2019" & trt == "11" ~ "Fervent",
    trt == "11S" ~ NA_character_, #2015
    trt == "13" ~ "Kenja",
    trt == "15" ~ "Mettle",
    trt == "16Q" ~ NA_character_, #2015
    trt == "18" ~ "UTC",
    trt == "1D" ~ NA_character_, #2015
    trt == "2" ~ "UBI-4319_4SC",
    trt == "29" ~ "UTC-JMS Stylet Oil",
    trt == "2Q" ~ "Quintec", #2015 - solid guess here
    trt == "3" ~ "Rhyme",
    trt == "5" ~ "Revus Top",
    trt == "Control" ~ "UTC",
    trt == "D" ~ "Mettle", #2017
    trt == "E" ~ "Rhyme", #2017
    trt == "Q" ~ "Topguard", #2017
    trt == "S" ~ "Luna Experience", #2017
    year == "2018" & trt == "UTC" ~ "UTC-JMS Stylet Oil",
    year != "2018" & trt == "UTC" ~ "UTC",
    TRUE ~ trt
  )) %>% 
  mutate(fung_chem = case_when(
    trt == "1" ~ "Flutriafol",
    trt == "10" ~ "Fluopyram and Tebuconazole", #same for 2018 and 2019
    year == "2018" & trt == "11" ~ "Fluopyram and Trifloxystrobin",
    year == "2019" & trt == "11" ~ "Isofetamid and Tebuconazole",
    trt == "11S" ~ "NA",
    trt == "13" ~ "Isofetamid",
    trt == "15" ~ "Tetraconazole",
    trt == "16Q" ~ NA_character_,
    trt == "18" ~ "UTC",
    trt == "1D" ~ NA_character_,
    trt == "2" ~ "UBI-4319_4SC",
    trt == "29" ~ "UTC-JMS Stylet Oil",
    trt == "2Q" ~ "Quinoxyfen",
    trt == "3" ~ "Flutriafol",
    trt == "5" ~ "Difenoconazole and Mandipropamid",
    trt == "Control" ~ "UTC",
    trt == "D" ~ "Tetraconazole",
    trt == "E" ~ "Flutriafol",
    trt == "Q" ~ "Azoxystrobin and Flutriafol",
    trt == "S" ~ "Fluopyram and Tebuconazole",
    year == "2018" & trt == "UTC" ~ "JMS Stylet Oil",
    year != "2018" & trt == "UTC" ~ NA_character_,
    TRUE ~ trt
  )) %>% 
  mutate(trt_simple=case_when(
    trt == "1" ~ paste0(year,"_3"),
    trt == "10" ~ paste0(year, "_3,7"),
    year == "2018" & trt == "11" ~ paste0(year, "_7,11"),
    year == "2019" & trt == "11" ~ paste0(year, "_3,7*"), # add * bc second one for 2019
    trt == "11S" ~ paste0(year, "_3,7"),
    trt == "13" ~ paste0(year, "_7"),
    trt == "15" ~ paste0(year, "_3*"), # add * bc second one for 2018
    trt == "16Q" ~ paste0(year, "_11"),
    trt == "18" ~ paste0(year, "_UTC"),
    trt == "1D" ~ paste0(year, "_3"),
    trt == "2" ~ paste0(year, "_3"),
    year == "2016" & trt == "29" ~ paste0(year, "_UTC"),
    year == "2018" & trt == "29" ~ paste0(year, "_UTC"), #actually JMS stylet oil 
    trt == "2Q" ~ paste0(year, "_13"),
    trt == "3" ~ paste0(year, "_3"),
    trt == "5" ~ paste0(year, "_3,40"),
    trt == "Control" ~ paste0(year, "_UTC"),
    trt == "D" ~ paste0(year, "_3"),
    trt == "E" ~ paste0(year, "_3*"), # add * bc second one for 2019
    trt == "Q" ~ paste0(year, "_3,11"),
    trt == "S" ~ paste0(year, "_3,7"),
    trt == "UTC" ~ paste0(year, "_UTC"),
    TRUE ~ trt
  )) %>% 
  mutate(wayne_trial_trt_number = case_when(
    trt == "11S" ~ NA_character_, #2015
    trt == "16Q" ~ NA_character_, #2015
    trt == "1D" ~ NA_character_, #2015
    trt == "2Q" ~ NA_character_, #2015
    trt == "UTC" ~ NA_character_, #2015
    trt == "Control" ~ "33", #2017
    trt == "D" ~ "25", #2017
    trt == "E" ~ "12", #2017
    trt == "Q" ~ "13", #2017
    trt == "S" ~ "10", #2017
    TRUE ~ trt
  ))
view(trial_trts)
unique(trial_trts$trt_simple)
#trial_trts %>% filter(year=="2015")

trial_trts_table <- trial_trts %>% 
  arrange(trt_simple) %>% 
  select(trt_simple,fung_name,fung_chem,wayne_trial_trt_number)
trial_trts_table

write_csv(trial_trts_table,"output/trial_trts_table.csv")


## Note: Did not perform clone corrections, AMOVA, etc.
