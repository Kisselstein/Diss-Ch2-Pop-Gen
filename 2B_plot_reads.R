

## Export and plot all reads from all loci and samples

#Load data
new_dat7b <- read_rds("output/database_clean.rds")
str(new_dat7b)

# Cleaned data 
reads1 <- new_dat7b %>%
  select(reads = reads1)
reads2 <- new_dat7b %>%
  select(reads = reads2)
reads_plot <- rbind.data.frame(reads1,reads2) %>% filter(!is.na(reads)) %>% arrange(reads)
#all reads
reads_plot2 <- reads_plot 
range(reads_plot2) #1 293507
ggplot(reads_plot2, aes(x = seq(1,nrow(reads_plot2)), y = reads)) +
  geom_point() +
  labs(title = "Reads: Cleaned Data", y = "Reads", x = "Samples") +
  theme_classic()
ggsave("output/reads_cleaned_all.tiff", dpi = 300, compression = "lzw")
#<50 reads
reads_plot2 <- reads_plot %>% filter(reads <= 50)
ggplot(reads_plot2, aes(x = seq(1,nrow(reads_plot2)), y = reads)) +
  geom_point() +
  labs(title = "Reads: Cleaned Data", y = "Reads (max 50)", x = "Samples") +
  theme_classic()
ggsave("output/reads_cleaned_max50.tiff", dpi = 300, compression = "lzw")

# Filtered data 
reads_plot <- dat_filtered2 %>% select(reads)  %>% filter(!is.na(reads)) %>% arrange(reads)
range(reads_plot) #10 293507
#all reads
ggplot(reads_plot, aes(x = seq(1,nrow(reads_plot)), y = reads)) +
  geom_point() +
  labs(title = "Reads: Filtered Data", y = "Reads", x = "Sample") +
  theme_classic()
ggsave("output/reads_filtered_all.tiff", dpi = 300, compression = "lzw")
#<50 reads
reads_plot2 <- reads_plot %>% filter(reads <= 50)
ggplot(reads_plot2, aes(x = seq(1,nrow(reads_plot2)), y = reads)) +
  geom_point() +
  labs(title = "Reads: Filtered Data", y = "Reads (max 50)", x = "Sample") +
  theme_classic()
ggsave("output/reads_filtered_max50.tiff", dpi = 300, compression = "lzw")



## Break down read stats by loci

#Cleaned data
new_dat7b <- read_rds("output/database_clean.rds")
str(new_dat7b)
cln_rds1 <- new_dat7b %>% 
  select(loci,reads=reads1)
head(cln_rds1)
cln_rds2 <- new_dat7b %>% 
  select(loci,reads=reads2)
head(cln_rds2)
cln_rds <- rbind.data.frame(cln_rds1,cln_rds2) %>%
  filter(!is.na(reads)) %>%
  group_by(loci) %>% 
  summarise(min(reads),
            max(reads),
            sum(reads)) %>% 
  arrange(loci)
head(cln_rds)
write_csv(cln_rds,"output/cln_reads.csv")


#Filtered data
# Filtered data 
dat_filtered2 <- readRDS("r_data/filtered_data.rds")
f_rds <- dat_filtered2 %>% 
  select(loci,reads) %>%
  filter(!is.na(reads)) %>%
  group_by(loci) %>% 
  summarise(min(reads),
            max(reads),
            sum(reads)) %>% 
  arrange(loci)
head(f_rds)
write_csv(f_rds,"output/fltr_reads.csv")


sum(cln_rds$`sum(reads)`) #271,660,739
sum(f_rds$`sum(reads)`)   #257,100,018


## LATER: Repeat but only for samples analyzed in Ch. 2

