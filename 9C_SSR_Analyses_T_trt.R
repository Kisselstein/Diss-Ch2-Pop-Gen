
## break down as populations (remove asterisk):

# Option 1:
# 3
# 7
# 11
# 13
# 3, 7
# 3, 11
# 7, 11
# 3, 40
# UTC

# Option 2:
# 3
# 7
# 11
# 13 - remove?
# UTC
# 3, 7 - assign 1
# 3, 11 - assign 1
# 7, 11 - assign 1
# 3, 40 - assign 1

# Option 3:
# 3
# 7
# 11
# 13 - remove?
# UTC
# 3, 7 - remove
# 3, 11 - remove
# 7, 11 - remove
# 3, 40 - remove

# Option 4:
# 3
# 7
# 11
# 13 - remove?
# UTC
# 3, 7 - remove
# 3, 11 - remove
# 7, 11 - remove
# 3, 40 - remove



chemclass_names <- test2table2 %>% 
  
  test2_table2
unique(test2_table2$Pop_New)
test2_table2 %>% filter(year=="2015") %>% filter() %>% select(Pop,year,trt,Pop_New)
test2_table2 %>% filter(year=="2016") %>% filter() %>% select(Pop,year,trt,Pop_New)
test2_table2 %>% filter(year=="2017") %>% filter() %>% select(Pop,year,trt,Pop_New)
test2_table2 %>% filter(year=="2018") %>% filter() %>% select(Pop,year,trt,Pop_New)
test2_table2 %>% filter(year=="2019") %>% filter() %>% select(Pop,year,trt,Pop_New)

write_csv(test2_table2,"output/poppr/full/chem_class_poppr.csv")



## PoppR Analyses:Trial (T) populations

## 1: Chemical Classes
## 2: Fungicide Name

#Load packages
library(tidyverse)
#library(readxl)
#library(lubridate)
library(poppr)
library(vegan)
#library(treemap)
#library(vegan)
#library(ape)


## First need to create a new file to import into PoppR with new names
t_dat <- read_csv("output/poppr/data/ssr_dat.csv") %>% as.data.frame()
ncol(t_dat) #15 columns

t_dat <- t_dat %>% 
  mutate(year=word(pop,1,1,"_"),
         ctw=word(pop,2,2,"_"),
         site=word(pop,3,3,"_"),
         variety=word(pop,4,4,"_"),
         trt=word(pop,5,5,"_"),
         timepoint=word(pop,6,6,"_")) %>% 
  filter(ctw=="T")
head(t_dat)
unique(t_dat$trt)



# 1: Treatment numbers only (remove asterisk)
t_dat1 <- t_dat %>% 
  mutate(pop_new=case_when(
    trt == "1" ~ paste0(year,"_3"),
    trt == "10" ~ paste0(year, "_3,7"),
    year == "2018" & trt == "11" ~ paste0(year, "_7,11"),
    year == "2019" & trt == "11" ~ paste0(year, "_3,7"), # add * bc second one for 2019
    trt == "11S" ~ paste0(year, "_3,7"),
    trt == "13" ~ paste0(year, "_7"),
    trt == "15" ~ paste0(year, "_3"), # add * bc second one for 2018
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
    trt == "E" ~ paste0(year, "_3"), # add * bc second one for 2019
    trt == "Q" ~ paste0(year, "_3,11"),
    trt == "S" ~ paste0(year, "_3,7"),
    TRUE ~ paste0(year,"_",trt)
  )) %>% 
  select(sample_id, pop=pop_new, EnMS2, EnCDnew_24,
         EnCDnew_27, EnCDnew_29, EnCDnew_30, EnCDnew_33,
         EnCDnew_35, EnCDnew_39, EnCDnew_41, EnMS8,
         EnCDnew_25, EnCDnew_32, EnMS6)
#str(t_dat1)
sort(unique(t_dat1$pop))
#head(t_dat1)
ncol(t_dat1) #needs to be 15 columns
write_csv(t_dat1,"output/poppr/data/ssr_dat_T_trts_1.csv")

## create file with sample no. per pop in this dataset
t_dat1_sample_n <- t_dat1 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
write_csv(t_dat1_sample_n,"output/poppr/data/ssr_dat_T_trts_1_sample_n.csv")

#open both files, update and save with "_genalex" pasted at end of file name
#13 loci
length(unique(t_dat1$sample_id)) #1662 samples
length(unique(t_dat1$pop)) #20 populations

# Load data
test2 <- read.genalex("output/poppr/data/ssr_dat_T_trts_1_genalex.csv",ploidy=1)
splitStrata(test2) <- ~year/trt

# Trts only
setPop(test2) <- ~trt
test2 #1172 mlgs, 1662 ind, 13 loci, 9 pops
popNames(test2)

test2tab <- mlg.table(test2, quiet = TRUE) 
test2tab[, 1:5] # Showing the first 5 columns and all rows of the table.
#str(test2tab)

test2_table <- poppr(test2) 
test2_table2 <- test2_table %>% as.data.frame() %>% arrange(Pop)
#test2_table2
write_csv(test2_table2,"output/poppr/full/T_trts_1_poppr.csv")

test2tab_df <- as.data.frame(test2tab) %>% rownames_to_column() %>%
  arrange(rowname) %>% column_to_rownames() %>% as.vector()
rownames(test2tab_df)
test2tab_df[,1:5]

#now finish the sort by rowname
rownames(test2tab_df)
nrow(test2tab_df)
test2tab_df <- test2tab_df[c(3,7,1,2,6,4,5,8,9),]
rownames(test2tab_df)
nrow(test2tab_df)
test2tab_df[,1:5]


# Use CF diversity function
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
test2stat2 <- as.data.frame(test2stat2) %>% rownames_to_column()
test2stat2
write_csv(test2stat2,"output/poppr/full/T_trts_1_div_stats_cf.csv")
# You can use filtered or custom MLGs to compare diversity

#rarefaction:eMLGs
rarecurve(test2tab_df, ylab="Number of expected MLGs", sample=min(rowSums(test2tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")
min(rowSums(test2tab)) #22
rowSums(test2tab)

# Div CI
test2_divc <- diversity_ci(test2tab_df, n = 100L, raw = FALSE) %>% 
  as.data.frame() %>% rownames_to_column()
test2_divc
write_csv(test2_divc,"output/poppr/full/T_trts_1_div_ci.csv")
# boxplots biased; read doc at ?diversity_ci
# ci can exist outside of possible range, like BB_10 and BB_11.

## Jack-knife rarefaction:
# diversity_ci() and diversity_boot() have the option
# meaning your data will be randomly sub-sampled to either 
# the smallest population size, or your specified n.rare value, 
# whichever is bigger. Here’s an example with the previous data set:

#DIV CI RARE (minimum)
#table
test2rare <- diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE) %>% 
  as.data.frame() %>% rownames_to_column()
test2rare # year, n=22
write_csv(test2rare,"output/poppr/full/T_trts_1_div_ci_rare22.csv")
#plot
diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)

p <- last_plot()
p + 
  labs(x = "Chemical Class", y = "value") + 
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings

ggsave("output/poppr/full/T_trts_1_div_ci_rare22.tiff", dpi = 300, compression = "lzw")






## 2: Fungicide Name
t_dat1 <- t_dat %>% 
  mutate(pop_new=case_when(
    trt == "1" ~ paste0(year,"_Rhyme"),
    trt == "10" ~ paste0(year, "_LunaExp"),
    year == "2018" & trt == "11" ~ paste0(year, "_LunaSens"),
    year == "2019" & trt == "11" ~ paste0(year, "_Fervent"), # add * bc second one for 2019
    trt == "11S" ~ paste0(year, "_NA"),
    trt == "13" ~ paste0(year, "_Kenja"),
    trt == "15" ~ paste0(year, "_Mettle"), # add * bc second one for 2018
    trt == "16Q" ~ paste0(year, "_NA"),
    trt == "18" ~ paste0(year, "_UTC"),
    trt == "1D" ~ paste0(year, "_NA"),
    trt == "2" ~ paste0(year, "_UBI-4319"),
    year == "2016" & trt == "29" ~ paste0(year, "_UTC"),
    year == "2018" & trt == "29" ~ paste0(year, "_UTC"), #actually JMS stylet oil 
    trt == "2Q" ~ paste0(year, "_NA"), #Probably Quintec
    trt == "3" ~ paste0(year, "_Rhyme"),
    trt == "5" ~ paste0(year, "_RevusTop"),
    trt == "Control" ~ paste0(year, "_UTC"),
    trt == "D" ~ paste0(year, "_Mettle"),
    trt == "E" ~ paste0(year, "_Rhyme"), # add * bc second one for 2019
    trt == "Q" ~ paste0(year, "_TopGuard"),
    trt == "S" ~ paste0(year, "_LunaExp"),
    TRUE ~ paste0(year,"_",trt)
  )) %>% 
  select(sample_id, pop=pop_new, EnMS2, EnCDnew_24,
         EnCDnew_27, EnCDnew_29, EnCDnew_30, EnCDnew_33,
         EnCDnew_35, EnCDnew_39, EnCDnew_41, EnMS8,
         EnCDnew_25, EnCDnew_32, EnMS6) %>% 
  filter(pop!="2015_NA")
#str(t_dat1)
sort(unique(t_dat1$pop))
#head(t_dat1)
ncol(t_dat1) #needs to be 15 columns
write_csv(t_dat1,"output/poppr/data/ssr_dat_T_trts_2.csv")

## create file with sample no. per pop in this dataset
t_dat1_sample_n <- t_dat1 %>% 
  group_by(pop) %>% 
  tally() %>% 
  ungroup() %>% 
  select(n,pop)
write_csv(t_dat1_sample_n,"output/poppr/data/ssr_dat_T_trts_2_sample_n.csv")

#13 loci
length(unique(t_dat1$sample_id)) #1540 samples
length(unique(t_dat1$pop)) #19 populations

###open both files, update and save with "_genalex" pasted at end of file name




# Load this new data
test2 <- read.genalex("output/poppr/data/ssr_dat_T_trts_2_genalex.csv",ploidy=1)
splitStrata(test2) <- ~year/trt
test2

# Trts only
setPop(test2) <- ~trt
test2 #1057 mlgs, 1508 ind, 13 loci, 10 pops
popNames(test2)

test2tab <- mlg.table(test2, quiet = TRUE) 
test2tab[, 1:5] # Showing the first 5 columns and all rows of the table.
#str(test2tab)

test2_table <- poppr(test2) 
test2_table2 <- test2_table %>% as.data.frame() %>% arrange(Pop)
test2_table2
write_csv(test2_table2,"output/poppr/full/T_trts_2_poppr.csv")

test2tab_df <- as.data.frame(test2tab) %>% rownames_to_column() %>%
  arrange(rowname) %>% column_to_rownames() %>% as.vector()
rownames(test2tab_df)
test2tab_df[,1:5]


# Use CF diversity function
# CF: clonal fraction statistic (# MLGs/N)
test2stat2 <- diversity_stats(test2tab_df, CF = myCF)
test2stat2 <- as.data.frame(test2stat2) %>% rownames_to_column()
test2stat2
write_csv(test2stat2,"output/poppr/full/T_trts_2_div_stats_cf.csv")
# You can use filtered or custom MLGs to compare diversity

#rarefaction:eMLGs
rarecurve(test2tab_df, ylab="Number of expected MLGs", sample=min(rowSums(test2tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")
min(rowSums(test2tab)) #22
rowSums(test2tab)

# Div CI
test2_divc <- diversity_ci(test2tab_df, n = 100L, raw = FALSE) %>% 
  as.data.frame() %>% rownames_to_column()
test2_divc
write_csv(test2_divc,"output/poppr/full/T_trts_2_div_ci.csv")
# boxplots biased; read doc at ?diversity_ci
# ci can exist outside of possible range, like BB_10 and BB_11.

## Jack-knife rarefaction:
# diversity_ci() and diversity_boot() have the option
# meaning your data will be randomly sub-sampled to either 
# the smallest population size, or your specified n.rare value, 
# whichever is bigger. Here’s an example with the previous data set:

#DIV CI RARE (minimum)
#table
test2rare <- diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE) %>% 
  as.data.frame() %>% rownames_to_column()
test2rare # year, n=22
write_csv(test2rare,"output/poppr/full/T_trts_2_div_ci_rare22.csv")
#plot
diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)

p <- last_plot()
p + 
  labs(x = "Treatment", y = "value") + 
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings

ggsave("output/poppr/full/T_trts_2_div_ci_rare22.tiff", dpi = 300, compression = "lzw")






##3: remove fungicides with <50 individuals
popNames(test2)
rowSums(test2tab)
test2tab <- mlg.table(test2, sublist = c(-3,-5,-9,-10), quiet = TRUE) 
test2tab[, 1:5] # Showing the first 5 columns and all rows of the table.
#str(test2tab)

test2_table <- poppr(test2) 
test2_table2 <- test2_table %>% as.data.frame() %>% arrange(Pop)
test2_table2
write_csv(test2_table2,"output/poppr/full/T_trts_3_poppr.csv")

test2tab_df <- as.data.frame(test2tab) %>% rownames_to_column() %>%
  arrange(rowname) %>% column_to_rownames() %>% as.vector()
rownames(test2tab_df)
test2tab_df[,1:5]


# Use CF diversity function
# CF: clonal fraction statistic (# MLGs/N)
test2stat2 <- diversity_stats(test2tab_df, CF = myCF)
test2stat2 <- as.data.frame(test2stat2) %>% rownames_to_column()
test2stat2
write_csv(test2stat2,"output/poppr/full/T_trts_3_div_stats_cf.csv")
# You can use filtered or custom MLGs to compare diversity

#rarefaction:eMLGs
rarecurve(test2tab_df, ylab="Number of expected MLGs", sample=min(rowSums(test2tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")
min(rowSums(test2tab)) #54
rowSums(test2tab)

# Div CI
test2_divc <- diversity_ci(test2tab_df, n = 100L, raw = FALSE) %>% 
  as.data.frame() %>% rownames_to_column()
test2_divc
write_csv(test2_divc,"output/poppr/full/T_trts_3_div_ci.csv")
# boxplots biased; read doc at ?diversity_ci
# ci can exist outside of possible range, like BB_10 and BB_11.

## Jack-knife rarefaction:
# diversity_ci() and diversity_boot() have the option
# meaning your data will be randomly sub-sampled to either 
# the smallest population size, or your specified n.rare value, 
# whichever is bigger. Here’s an example with the previous data set:

#DIV CI RARE (minimum)
#table
test2rare <- diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE) %>% 
  as.data.frame() %>% rownames_to_column()
test2rare # year, n=54
write_csv(test2rare,"output/poppr/full/T_trts_3_div_ci_rare54.csv")
#plot
diversity_ci(test2tab_df, n = 100L, rarefy = TRUE, raw = FALSE)

p <- last_plot()
p + 
  labs(x = "Treatment", y = "value") + 
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings

ggsave("output/poppr/full/T_trts_3_div_ci_rare54.tiff", dpi = 300, compression = "lzw")
