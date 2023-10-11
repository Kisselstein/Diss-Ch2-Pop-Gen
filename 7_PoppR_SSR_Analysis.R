
#Load packages
library(tidyverse)
library(readxl)
library(lubridate)
library(poppr)
library(treemap)
library(vegan)
library(ape)


# Load data
ssr1 <- read.genalex("output/ssr_dat_final.csv",ploidy=1)
splitStrata(ssr1) <- ~year/ctw/site/variety/trt/timepoint
ssr1 
#ssr1 unf: 3706 mlg 4616 ind 13 loci 6 strata 104 pops


#resample.ia(ssr1,reps=100,hist=TRUE)
pair.ia(ssr1)


#Genotype Accumulation Curve
#assess power to discriminate between ind. given a random sample of n loci
tiff("output/poppr/sp_genotype_curve.tiff", width = 400, height = 400)
gac <- genotype_curve(ssr1, sample = 1000, quiet = TRUE)
dev.off()
#for each boxplot, n loci were randomly sampled 1000 times in order to create the distribution
#we see plateau at 12 loci?


# Missing data in the population heat map and table
setPop(ssr1) <- ~year/ctw
ssr1 #14 pops

tiff("output/poppr/missing_heatmap_yc.tiff", width = 800, height = 800)
info_table(ssr1, type = "missing", plot = TRUE)
dev.off() 
#table
ssr_info_table <- as.data.frame(info_table(ssr1, type = "missing", plot = FALSE)) %>% 
  rownames_to_column()
write_csv(ssr_info_table,"output/poppr/missing_yc.csv")


#Allele frequencies, missing data, and ploidy
#look for missing data, rare alleles and overall data quality
(ssr1lt <- locus_table(ssr1))
# have anywhere from 3-14 microsat alleles per locus
ssr1lt <- as.data.frame(ssr1lt) %>% rownames_to_column()
head(ssr1lt)
write_csv(ssr1lt,"output/poppr/sp_all_freq.csv")





## Fix replen (repeat motif length)

unique(ssr1$loc.fac) #loci and their order
#MS2 3, NEW24 4, NEW27 4, NEW29 4, NEW30 4, NEW33 4, NEW35 4, NEW39 4, NEW41 5, MS8 3, NEW25 4, NEW32 4, MS6 3

locus_replengths <- c(3,4,4,4,4,4,4,4,5,3,4,4,3) #repeat motif lengths for all 16 loci
ssr1reps <- fix_replen(ssr1, c(locus_replengths)) # #bp for each locus
#The repeat lengths for EnCDnew_29, EnCDnew_41 are not consistent.
#Repeat lengths with some modification are being returned: EnCDnew_24



ssr1 <- missingno(ssr1,type="geno",cutoff=0.5)
ssr1 
#ssr1 unf: 3706 mlg 4616 ind
#ssr1 fil: 2815 mlg 3166 ind

ssr1 <- missingno(ssr1,type="loci",cutoff=0.5)
ssr1 #no change




ssr1_filtered <- filter_stats(ssr1, distance = bruvo.dist, replen = ssr1reps, plot = FALSE, missing = "ignore")
#plot_filter_stats(x, fstats, distmat, cols = NULL, nclone = NULL, breaks = NULL)


####################
bres <- filter_stats(Pinf, distance = bruvo.dist, replen = pinfreps, plot = TRUE, threads = 1L)
print(bres) # shows all of the statistics

# Use these results with cutoff_filter()
print(thresh <- cutoff_predictor(bres$farthest$THRESHOLDS))
mlg.filter(Pinf, distance = bruvo.dist, replen = pinfreps) <- thresh
Pinf 

# Different distances will give different results -----------------------
nres <- filter_stats(Pinf, distance = nei.dist, plot = TRUE, threads = 1L, missing = "mean")
print(thresh <- cutoff_predictor(nres$farthest$THRESHOLDS))
mlg.filter(Pinf, distance = nei.dist, missing = "mean") <- thresh
Pinf 
###################



# Use this function to find cutoff if you don't see a small initial peak!
print(farthest_thresh <- cutoff_predictor(ssr1_filtered$farthest$THRESHOLDS))
print(average_thresh  <- cutoff_predictor(ssr1_filtered$average$THRESHOLDS))
print(nearest_thresh  <- cutoff_predictor(ssr1_filtered$nearest$THRESHOLDS))

# Now define MLL for ssr1 with the following criteria:
#[t]	threshold	0.017 doesn't change y / ctw / s
#[d]	distance	Bruvo’s Distance
#[a]	algorithm	Farthest neighbor

#doesn't work, I have to treat missing first for some reason!?
ssr1_filter_value <- mlg.filter(ssr1, missing = "ignore", distance = bruvo.dist, replen = ssr1reps, 
                                algorithm = "f") <- 0.017
ssr1_filter_value #0.017
is output very critically

# ok jumping around
# now create a minimum spanning network for MLLs
#mll(ssr1)
#LETTERS[mll(ssr1)]
#mll.custom(ssr1) <- LETTERS[mll(ssr1)]
#mlg.table(ssr1)
#ssr1pal <- colorRampPalette(c("blue", "gold"))
#set.seed(9001)
#ssr1msn <- bruvo.msn(ssr1, replen = c(1,1,1,1), palette = ssr1pal,
#                     vertex.label.color = "firebrick", vertex.label.font = 2,
#                     vertex.label.cex = 1.5)
#mll.levels(ssr1)[mll.levels(ssr1) == "C"] <- "B" #Manually make X3+X5 same too


## DIVERSITY

popNames(ssr1)

ssr1tab <- mlg.table(ssr1) #save as tiff

rarecurve(ssr1tab, ylab="Number of eMLGs", sample=min(rowSums(ssr1tab)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")


## Remove certain sites/datapoints
min(rowSums(ssr1tab)) 
rowSums(ssr1tab) 
#ssr1tab <- mlg.table(ssr1, sublist = c(-1,-2,-12,-13,-14,-15,-26,-27))
#ssr1tab <- mlg.table(ssr1, sublist = c(-1,-2,-3,-9,-10,-11,-(15:20),-(26:30),-(36:41),-(47:53))) #save as tiff

#min value from rowSums() of the table is the min common sample size of all pops 
#Setting the “sample” flag draws the horizontal and vertical lines you see on the graph. 
#The intersections of these lines = eMLG when running poppr function

# Bar chart for 1 pop
#mlg.table(ssr1, sublist = "2020_W", plot = TRUE) #barchart for 1 pop
#mlg.table(ssr1, strata = ~year/ctw/site, sublist = "2020_W_W", plot = TRUE) #barchart for 1 pop

#v.tab <- mlg.table(ssr1, plot = FALSE)
#v.tab[1:5, 1:5] # Showing the first 5 columns and rows of the table.
#mlg.table(ssr1, sublist = "2020_W",plot = TRUE)

ssr1stat <- as.data.frame(diversity_stats(ssr1tab))
head(ssr1stat)
ssr1_divci <- diversity_ci(ssr1tab, n = 100L, raw = FALSE)
head(ssr1_divci)
write_csv(ssr1_divci,"output/poppr/1520_ctw_s_div_ci.csv")

# boxplots biased; read doc at ?diversity_ci
# be careful when interpreting bc confidence the ci can 
# exist outside of possible range, like BB_10 and BB_11.


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

ssr1stat_CF <- as.data.frame(diversity_stats(ssr1tab, CF = myCF))
head(ssr1stat_CF)
write_csv(ssr1stat_CF,"output/poppr/div_s_cf.csv")
# You can use filtered or custom MLGs to compare diversity
# Here, I’m filtering genotypes in monpop that are different by only a single mutational step (Bruvo et al. 2004).




## Jack-knife rarefaction

# diversity_ci() and diversity_boot() have the option
# meaning your data will be randomly sub-sampled to either 
# the smallest population size, or your specified n.rare value, 
# whichever is bigger. Here’s an example with the previous data set:

#minimum
ssr1rare <- diversity_ci(ssr1tab, n = 100L, rarefy = TRUE, raw = FALSE)
#13
head(ssr1rare)
write_csv(ssr1rare,"output/poppr/div_ci_rare.csv")
#save tiff 1520_ctw_s_noA_div_ci_rare209

tiff("output/poppr/cr_yst_div_ci_rare.tiff", width = 1600, height = 1000)
ssr1rare <- diversity_ci(ssr1tab, n = 100L, rarefy = TRUE, raw = FALSE)
dev.off()







ssr1amova <- poppr.amova(ssr1, ~ctw/site, filter = TRUE, missing = "ignore", 
                            dist = bruvo.dist, algorithm = "f") #rare10
ssr1amova

