### Script 8: Subset Test SSR Analyses (TESTING!)
## 202308

## Find where we plotted alleles per locus and things like that
## Get nearest and farthest distance algorithms to actually run 

#Load packages
library(tidyverse)
library(readxl)
library(lubridate)

## Load in subset data
ssr_dat <- read_csv("output/poppr/data/ssr_dat_subset.csv", col_names = FALSE)
head(ssr_dat)

names(ssr_dat) <- ssr_dat[3,]
ssr_dat <- ssr_dat[c(-1,-2,-3),]
ssr_dat <- ssr_dat %>%
  pivot_longer(!c(sample_id,pop), names_to = "loci", values_to = "all")
head(ssr_dat)



## Loci Plots

## number of alleles per locus
allperloc <- ssr_dat %>% select(loci,all) %>% distinct() %>% 
  group_by(loci) %>% tally() %>% ungroup()
allperloc
range(allperloc$n)
mean(allperloc$n)

ggplot(allperloc, aes(x = loci, y = n)) + geom_col() +
  labs(x = "Loci", y = "Number of Alleles") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/subset/loci_numofall.tiff", dpi = 300, compression = "lzw")



allplot <- ssr_dat %>% 
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
ggsave("output/poppr/subset/ssr_allele_freqs.tiff", dpi = 300, compression = "lzw")


allplot2 <- ssr_dat %>% 
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
  labs(x = "Site", fill = "Allele", y = "%") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/subset/encdnew33_allele_freq_by_year_site.tiff", dpi = 300, compression = "lzw")



allplot2a <- ssr_dat %>% 
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
  labs(x = "Site", fill = "Allele", y = "%") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/subset/encdnew30_allele_freq_by_year_site.tiff", dpi = 300, compression = "lzw")




allplot2b <- ssr_dat %>% 
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
  labs(x = "Site", fill = "Allele", y = "%") + 
  theme_bw() + scale_fill_viridis_d(option="turbo") +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/subset/encdnew41_allele_freq_by_year_site.tiff", dpi = 300, compression = "lzw")







# look at het and success rate for our ssr loci
#ssr_het <- write_csv(ssr_het,"output/1520_ssr_het_all_samples.csv")
dat <- readRDS("r_data/filtered_data_mat_clonal_allsize.rds")
head(dat)

ssr_dat2 <- as.data.frame(ssr_dat) %>% 
  filter(!is.na(all)) %>% mutate(sample_id_new = case_when(
    str_detect(sample_id,"-dup") ~ word(sample_id,1,-2,"-"),
    TRUE ~ sample_id
  )) %>% #ssr_dat3 %>% filter(str_detect(sample_id_new,"-dup"))
  select(sample_id=sample_id_new,loci,all) %>% 
  distinct() %>% 
  group_by(sample_id,loci) %>% 
  tally() %>% 
  ungroup() %>% #unique(ssr_dat3$n)
  mutate(het=case_when(
    n==1 ~ FALSE,
    n==2 ~ TRUE
  )) %>% 
  select(-n)
head(ssr_dat2)

ssr_het1 <- ssr_dat2 %>% filter(het==TRUE) %>% group_by(loci) %>% tally() %>% ungroup() %>% select(loci,het=n)
ssr_het2 <- ssr_dat2 %>% filter(het==FALSE) %>% group_by(loci) %>% tally() %>% ungroup() %>% select(loci,hom=n)
ssr_het <- full_join(ssr_het1,ssr_het2, by = c('loci')) %>% 
  mutate(sum=het+hom,
         perc_het=het/sum)
head(ssr_het)
range(ssr_het$sum)
signif(range(ssr_het$perc_het), digits = 2)


ggplot(ssr_het, aes(x = loci, y = perc_het)) + geom_col() +
  labs(x = "Loci", y = "Heterozygosity Rate (%)") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/poppr/subset/loci_het_rate.tiff", dpi = 300, compression = "lzw")


# SKIP: clonality: probably not important, right?


# locus: success and failure rate

ssr_loc <- as.data.frame(ssr_dat) %>% 
  mutate(sample_id_new = case_when(
    str_detect(sample_id,"-dup") ~ word(sample_id,1,-2,"-"),
    TRUE ~ sample_id
  )) %>% #ssr_loc_nodup %>% filter(str_detect(sample_id_new,"-dup"))
  mutate(all_dat = case_when(
    !is.na(all) ~ "Y",
    is.na(all) ~ "N"
  )) %>% 
  select(loci,sample_id_new,all_dat) %>% 
  distinct() %>%  
  group_by(loci,all_dat) %>% 
  tally() %>% 
  pivot_wider(names_from = all_dat, values_from = n) %>% 
  mutate(sum=Y+N,
         perc_pass = Y/sum,
         perc_fail = N/sum)
head(ssr_loc)


# success rate

ggplot(ssr_loc, aes(x = loci, y = perc_pass)) + geom_col() +
  labs(x = "Loci", y = "Success Rate (%)") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/loci_success_rate.tiff", dpi = 300, compression = "lzw")

# failure rate

ggplot(ssr_loc, aes(x = loci, y = perc_fail)) + geom_col() +
  labs(x = "Loci", y = "Failure Rate (%)") + 
  theme_bw() +
  theme(text = element_text(size=13), #x and y titles
        axis.text.x=element_text(angle=45,hjust=1,size=11), #x axis markings
        axis.text.y=element_text(size=11)) #y axis markings
ggsave("output/loci_failure_rate.tiff", dpi = 300, compression = "lzw")








## Missingness Heatmap

#1: group samples by year
ssr_missing <- as.data.frame(ssr_dat) %>% 
  mutate(sample_id_new = case_when(
    str_detect(sample_id,"-dup") ~ word(sample_id,1,-2,"-"),
    TRUE ~ sample_id
  )) %>% #ssr_loc_nodup %>% filter(str_detect(sample_id_new,"-dup"))
  mutate(all_dat = case_when(
    !is.na(all) ~ "Y",
    is.na(all) ~ "N"
  )) %>% 
  mutate(grp=word(pop,1,1,"_")) %>% 
  select(grp,loci,sample_id_new,all_dat) %>% 
  distinct() %>%  #6747 before, 4888 after (worked correctly)
  group_by(loci,all_dat,grp) %>% 
  tally() %>% 
  pivot_wider(names_from = all_dat, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(sum=Y+N,
         perc_pass = Y/sum,
         perc_fail = N/sum)
head(ssr_missing)

ggplot(ssr_missing, aes(grp, loci)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = perc_fail)) + 
  scale_fill_viridis_c()
ggsave("output/poppr/subset/heatmap_perc_fail_by_yr.tiff", dpi = 300, compression = "lzw")

ggplot(ssr_missing, aes(grp, loci)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = perc_pass)) + 
  scale_fill_viridis_c()
ggsave("output/poppr/subset/heatmap_perc_fail_by_yr.tiff", dpi = 300, compression = "lzw")


# Heatmap 2: group samples by year and CTW
ssr_missing <- as.data.frame(ssr_dat) %>% 
  mutate(sample_id_new = case_when(
    str_detect(sample_id,"-dup") ~ word(sample_id,1,-2,"-"),
    TRUE ~ sample_id
  )) %>% #ssr_loc_nodup %>% filter(str_detect(sample_id_new,"-dup"))
  mutate(all_dat = case_when(
    !is.na(all) ~ "Y",
    is.na(all) ~ "N"
  )) %>% 
  mutate(grp=word(pop,1,2,"_")) %>% 
  select(grp,loci,sample_id_new,all_dat) %>% 
  distinct() %>%  #6747 before, 4888 after (worked correctly)
  group_by(loci,all_dat,grp) %>% 
  tally() %>% 
  pivot_wider(names_from = all_dat, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(sum=Y+N,
         perc_pass = Y/sum,
         perc_fail = N/sum)
head(ssr_missing)

ggplot(ssr_missing, aes(grp, loci)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = perc_fail)) + 
  scale_fill_viridis_c()
ggsave("output/poppr/subset/heatmap_perc_fail_by_yr_ctw.tiff", dpi = 300, compression = "lzw")

ggplot(ssr_missing, aes(grp, loci)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = perc_pass)) + 
  scale_fill_viridis_c()
ggsave("output/poppr/subset/heatmap_perc_fail_by_yr_ctw.tiff", dpi = 300, compression = "lzw")








## PoppR

#Load packages
library(poppr)
#library(treemap)
#library(vegan)
#library(ape)

#sort data by pop column 

# Load data
ssr1 <- read.genalex("output/poppr/data/ssr_dat_subset.csv",ploidy=1)
splitStrata(ssr1) <- ~year/ctw/site/variety/trt/timepoint
setPop(ssr1) <- ~year/ctw/site/variety/trt/timepoint
popNames(ssr1) #11

ssr1 #491 mlg, 519 ind 13 loci

#tab(ssr1)
#nmll(ssr1) #491
#mll(ssr1) #519 MLLs
mlg(ssr1) #491 mlgs, 519 ind
levels(ssr1$pop) #11 pops
levels(ssr1$loc.fac) #13 loci (note order)

locus_replengths <- c(4,3,3,4,4,4,4,4,4,4,4.5,3,4) #repeat motif lengths for all 13 loci, in order
ssr1reps <- fix_replen(ssr1, c(locus_replengths)) # #bp for each locus
#The repeat lengths for EnCDnew_41, EnCDnew_29 are not consistent.
#Repeat lengths with some modification are being returned: EnCDnew_24



ssr1 #491 mlg, 519 ind 13 loci
test <- missingno(ssr1,type="geno",cutoff=0.5)
test #359 mlg, 378 ind 13 loci
test2 <- missingno(test,type="loci",cutoff=0.5)
test2 #12 loci (removed EnCDnew_35)
ssr1_filtered <- filter_stats(test2, distance = bruvo.dist, replen = ssr1reps, 
                              plot = TRUE, missing = "ignore")

# Use this function to find cutoff if you don't see a small initial peak!
print(farthest_thresh <- cutoff_predictor(ssr1_filtered$farthest$THRESHOLDS))
print(average_thresh  <- cutoff_predictor(ssr1_filtered$average$THRESHOLDS))
print(nearest_thresh  <- cutoff_predictor(ssr1_filtered$nearest$THRESHOLDS))

# Now define MLL for ssr1 with the following criteria:
#[t]	threshold	0.0208 (farthest or nearest is 0.0208, avg is 0.0104)
#[d]	distance	Bruvo’s Distance
#[a]	algorithm	Farthest neighbor

#mlg.filter(ssr1, missing = "ignore", distance = bruvo.dist, replen = ssr1reps, algorithm = "f") <- farthest_thresh
value <- mlg.filter(test2, missing = "ignore", distance = bruvo.dist, replen = ssr1reps, algorithm = "f") <- farthest_thresh
value #w late comm 0.016
#wout late comm 0.016
#Evaluate this output very critically

mlg.filter(test2, distance = bruvo.dist, replen = ssr1reps, algorithm = "f") <- farthest_thresh
test2 #311 contracted mlgs




## DIVERSITY

# c only
setPop(test2) <- ~ctw
popNames(test2)
test2tab <- mlg.table(test2)
test2tab[1:3, 1:5] # Showing the first 5 columns and 3 rows of the table.

#ycs 
setPop(test2) <- ~year/ctw/site
popNames(test2)
test2tab <- mlg.table(test2)
# remove 2019 C B1
test2tab <- mlg.table(test2, sublist = c(-3))
test2tab[,1:5] # Showing the first 5 columns and all rows of the table.


#rarefaction
library(vegan)
rarecurve(test2tab, ylab="Number of expected MLGs", sample=min(rowSums(test2tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")
min(rowSums(test2tab)) 
rowSums(test2tab)
# remove 2018 samples for this analysis since low numbers?

mlg.table(test2, sublist = "2019_W_W") #barchart for 1 pop
mlg.table(test2, strata = ~ctw, sublist = "W", plot = TRUE) #barchart for 1 pop, diff strata

v.tab <- mlg.table(test2, plot = FALSE)
v.tab[1:5, 1:5] # Showing the first 5 columns and rows of the table.

test2_divc <- diversity_ci(test2tab, n = 100L, raw = FALSE)
test2_divc
write_csv(test2_divc,"output/poppr/subset/div_ci.csv")
# boxplots biased; read doc at ?diversity_ci
# ci can exist outside of possible range, like BB_10 and BB_11.


# lots of diversity_* functions to choose from
# here's how to build custom (MLG/n)
# write it for both a matrix and a vector of counts 
# if you want to be able to bootstrap it
myCF <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){ # if it's a matrix
    res <- rowSums(x > 0)/rowSums(x)
  } else {                 # if it's a vector
    res <- sum(x > 0)/sum(x)
  }
  return(res)
}

#test2stat2 <- diversity_stats(test2tab, CF = myCF)
#test2stat2
#test2stat2 <- as.data.frame(test2stat2)
#write_csv(test2stat2,"output/poppr/subset/div_stats_cf.csv")
# You can use filtered or custom MLGs to compare diversity




## Jack-knife rarefaction

# diversity_ci() and diversity_boot() have the option
# meaning your data will be randomly sub-sampled to either 
# the smallest population size, or your specified n.rare value, 
# whichever is bigger. Here’s an example with the previous data set:

#minimum
test2rare <- diversity_ci(test2tab, n = 100L, rarefy = TRUE, raw = FALSE)
test2rare # ycs, n=10
write_csv(test2rare,"output/poppr/subset/div_ci_rare10.csv")

test2rare <- diversity_ci(test2tab, n = 100L, rarefy = TRUE, n.rare = 5, raw = FALSE)
test2rare # ycs, n = 9 is the minimum
write_csv(test2rare,"output/poppr/subset/div_ci_rare9.csv")

test2rare <- diversity_ci(test2tab, n = 100L, rarefy = TRUE, n.rare = 20, raw = FALSE)
test2rare # ycs, n = 20
write_csv(test2rare,"output/poppr/subset/div_ci_rare20.csv")


test2_poppr <- poppr(test2, CF = myCF) #wait look back at sample sizes?
test2_poppr
write_csv(test2_poppr, "output/poppr/subset/poppr.csv")
#This can give you comparable estimates of diversity when not all samples are of equal size.




?resample.ia
resample.ia(test2,reps=99,quiet=TRUE,use_psex=FALSE,plot=TRUE) #no plot??
pair.ia(test2)
pair.ia(test2,index="Ia",quiet=TRUE,plot=TRUE)

#mlg.filter(test2, missing = "ignore", distance = bruvo.dist, replen = ssr1reps, algorithm = "f") <- farthest_thresh

test2amova <- poppr.amova(test2, ~ctw, filter = TRUE, missing = "ignore", 
                          dist = bruvo.dist, algorithm = "f") #rare10
test2amova







## Clone correction - addtl code

#order of the samples affect the sampling, 
#we will take the sum of all pairwise distances 
#between clone-corrected samples 
#(corrected without respect to populations):

ssr1 %>%
  clonecorrect(strata = ~ctw) %>%  # 1. clone correct whole data set
  dist() %>%                     # 2. calculate distance
  sum()  
set.seed(999)
ssr1[sample(nInd(ssr1))] %>% # 1. shuffle samples
  clonecorrect(strata = ~ctw) %>%  # 2. clone correct whole data set
  dist() %>%                     # 3. calculate distance
  sum()                          # 4. take the sum of the distance





## Permutations and bootstrap resampling
# A common Ho for pops with mixed rep. modes is panmixia (lots of sex)
# Poppr randomly shuffles data sets in order to calculate P-values 
# for the index of association (𝐼𝐴) (Agapow and Burt 2001) 
# using 4 different methods:
#   method	strategy	units sampled
#1	permutation	alleles
#2	simulation	alleles
#3	simulation	alleles
#4	permutation	genotypes

#Function: shufflepop
shufflepop(test2, method = 1)
#pop - a genind object.
#method - a number indicating the method of sampling you wish to use
#The following methods are available for use:
#1.Permute Alleles (default) 
#2.Parametric Bootstrap 
#4.Multilocus permutation 
#Named this bc same method as the permutation analysis in program (Agapow and Burt 2001). 
#This will shuffle the genotypes at each locus.
#Note you have the same genotypes after shuffling, so at each locus, 
#you will maintain the same allelic freqs and heterozygosity. 
#So, in this sample, you will only see a homozygote with allele 2. 
#This ensures that the P-values for 𝐼𝐴 and 𝑟¯𝑑 are exactly the same. 
#This method assumes that alleles are not independently assorting within individuals. 
#This strategy is useful if you suspect the pop is inbreeding (Jerome Goudet, personal communication).

#These shuffling schemes have been implemented for the index of association, 
#but there may be other summary statistics you can use shufflepop for. 
#All you have to do is use the function replicate. 
#ex. use average Bruvo's distance with the first pop of nancycats

test2_sub <- popsub(test2, 1) #pick certain pops to keep
test2_sub
test2
observed <- mean(bruvo.dist(test2_sub, replen = ssr1reps))
observed #0.33

set.seed(9999)
bd.test <- replicate(99, mean(bruvo.dist(shufflepop(test2_sub, method = 2), replen = ssr1reps)))

hist(bd.test, xlab = "Bruvo's Distance", main = "Average Bruvo's distance over 99 randomizations")
abline(v = observed, col = "red")
legend('topleft', legend="observed", col="red", lty = 1)


## Removing uninformative loci
#Phylogenetically uninformative loci have only one sample differentiating from the rest. 
#This can lead to biased results when using multilocus analyses such as the index of association (Brown, Feldman, and Nevo 1980; Smith et al. 1993). 
#These nuisance loci can be removed with the following function.

test2
informloci(test2, cutoff = 0.1, MAF = 0.01, quiet = FALSE) #all sites polymorphic
#cutoff 0.1 = 38 samples


## MLG CROSSPOP ANALYSIS: same MLG across pops? Track like this!!!
test2.dup <- mlg.crosspop(test2, quiet = FALSE)
head(test2.dup) #number of copies in each population

#count the num of pops each MLG crosses using sapply to loop over the data with length.
test2.num <- sapply(test2.dup, length) # count the number of populations each MLG crosses.
head(test2.num)


## Combining MLG functions
#Let's find which MLGs were dup across pops


#setPop(test2) <- ~year/ctw/site
UGNN.list <- c(popNames(test2))
UGNN.list
UGNN <- mlg.crosspop(test2, sublist=UGNN.list, indexreturn=TRUE)
UGNN #shows number MLGs crossing between these pops

#but how many are in each? subset our original table, v.tab.
# and whats the incidence of these MLGs throughout our data set?
#Let's use mlg.vector to find individuals corresponding to the MLGs. 
#First we”ll investigate what the output of this function looks like.

test2.vec <- mlg.vector(test2)
str(test2.vec) # Analyze the structure.

#integers are MLG assignment of each ind in the same order as dataset. 
#first two individuals have same alleles at each locus, so both are MLG: 605. 

length(unique(test2.vec)) # 1905 MLGs
test2 # 1905 MLGs

#We will take UGNN (MLGs crossing UK, Germany, Netherlands, and Norway) 
#and compare its elements to the MLG vector (v.vec) to see where else they occur.

UGNN # Show what we are looking for
UGNN_match <- test2.vec %in% UGNN
table(UGNN_match) # How many individuals matched to those crosspop MLGs?
## UGNN_match
## FALSE  TRUE 
##  1881    22
#22 inds matched to those three MLGs

view(indNames(test2)[UGNN_match])
length(indNames(test2)[UGNN_match]) #73

#alternative way to list individuals matching specific MLGs using the function mlg.id. 
#Each element in the list is named with the MLG, but the index does not necessarily match up, 
#so it is important to convert your query MLGs to strings:

test2.id <- mlg.id(test2)
test2.id[as.character(UGNN)]

#use the vector of MLGs to subset mlg.table() with the mlgsub flag.
mlg_cross <- mlg.table(test2, mlgsub = UGNN)
mlg_cross
write_csv(mlg_cross,"output/poppr/subset/crosspop.csv")
#That showed us exactly which pops these MLGs came from in our data set.



## Appendix
## Manipulating Graphics (ggplot)
## Exporting Graphics
## References

?poppr.amova
test2amova <- poppr.amova(test2, ~year/site, filter = TRUE, missing = "ignore", 
                          dist = bruvo.dist, algorithm = "f") #rare10
test2amova
#This can give you comparable estimates of diversity when not all samples are of equal size.

