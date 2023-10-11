
## Phase out Haplotypes using the remaining SSR loci from Prep
# Continue with data from 5_Prep_Phasing_Haplotypes.R

## Q. Wait, why can't we keep loci that are het in more than half of samples?
## Q. Isn't it still possible to perform phasing??

## Phase out haplotypes for aclonal samples, and
## Fix het loci for clonal samples

str(full_data_clonal1)
ssr_loci_het
ssr_dat <- full_data_clonal1 %>% 
  filter(loci %in% ssr_loci_keep & #keep only ssr loci
           !loci %in% ssr_loci_het) #remove 4 loci with >50% het across all samples
head(ssr_dat)
unique(ssr_dat$loci) #13 remaining 

#pt 1 - clonal samples, hom loci
ssr_dat_pt1 <- ssr_dat %>% 
  filter(perc_clonal > 0.75 & het == "FALSE") %>% 
  mutate(len=as.character(all_len)) %>% 
  select(loci,sample_id,all,len)
head(ssr_dat_pt1)  

#pt 2 - clonal samples, het loci
ssr_dat_pt2 <- ssr_dat %>% 
  filter(perc_clonal > 0.75 & het == "TRUE") %>% 
  select(loci,sample_id,all,reads,all_len) %>% 
  arrange(sample_id,loci)
head(ssr_dat_pt2)

ssr_dat_pt2a <- ssr_dat_pt2[seq(1, nrow(ssr_dat_pt2), 2), ]
ssr_dat_pt2b <- ssr_dat_pt2[- 1, ]
ssr_dat_pt2b <- ssr_dat_pt2b[seq(1, nrow(ssr_dat_pt2b), 2), ]
ssr_dat_pt2a <- ssr_dat_pt2a %>% 
  select(loci,sample_id,all1=all,reads1=reads,len1=all_len)
ssr_dat_pt2b <- ssr_dat_pt2b %>% 
  select(loci,sample_id,all2=all,reads2=reads,len2=all_len)
head(ssr_dat_pt2a)
head(ssr_dat_pt2b)

ssr_dat_pt2_final <- full_join(ssr_dat_pt2a,ssr_dat_pt2b)
ssr_dat_pt2_final <- ssr_dat_pt2_final %>% 
  mutate(all = case_when(
    reads1 > reads2 ~ all1,
    reads1 < reads2 ~ all2
  )) %>% 
  mutate(reads = case_when(
    reads1 > reads2 ~ reads1,
    reads1 < reads2 ~ reads2
  )) %>% 
  mutate(len = case_when(
    reads1 > reads2 ~ len1,
    reads1 < reads2 ~ len2
  )) %>% 
  mutate(len=as.character(len)) %>% 
  select(loci,sample_id,all,len)
head(ssr_dat_pt2_final)

# pt 3 - aclonal samples
# hom loci
ssr_dat_pt3a <- ssr_dat %>% 
  filter(perc_clonal <= 0.75  & het == "FALSE") %>% 
  select(loci,sample_id,all,len=all_len) %>% 
  arrange(sample_id,loci)
head(ssr_dat_pt3a)  
# hom loci again; duplicate with new sample names
ssr_dat_pt3b <- ssr_dat %>% 
  filter(perc_clonal <= 0.75 & het == "FALSE") %>% 
  mutate(id = paste0(sample_id,"-","dup")) %>% 
  select(loci,sample_id=id,all,len=all_len) %>% 
  arrange(sample_id,loci)
head(ssr_dat_pt3b)  
# het loci
ssr_dat_pt3c <- ssr_dat %>% 
  filter(perc_clonal <= 0.75 & het == "TRUE") %>% 
  select(loci,sample_id,all,reads,len=all_len) %>% 
  arrange(loci,sample_id)
head(ssr_dat_pt3c)  

ssr_dat_pt3d <- ssr_dat_pt3c[seq(1, nrow(ssr_dat_pt3c), 2), ]
ssr_dat_pt3e <- ssr_dat_pt3c[- 1, ]
ssr_dat_pt3e <- ssr_dat_pt3e[seq(1, nrow(ssr_dat_pt3e), 2), ]

ssr_dat_pt3d <- ssr_dat_pt3d %>% 
  select(loci,sample_id,all,len)
ssr_dat_pt3e <- ssr_dat_pt3e %>% 
  mutate(id=paste0(sample_id,"-dup")) %>% 
  select(loci,sample_id=id,all,len)
head(ssr_dat_pt3d)
head(ssr_dat_pt3e)

# Combine all data
head(ssr_dat_pt1,2) 
head(ssr_dat_pt2_final,2)
head(ssr_dat_pt3a,2) 
head(ssr_dat_pt3b,2) 
head(ssr_dat_pt3d,2) 
head(ssr_dat_pt3e,2)

ssr_dat_combined <- ssr_dat_pt1 %>% 
  rbind.data.frame(ssr_dat_pt2_final) %>% 
  rbind.data.frame(ssr_dat_pt3a) %>% 
  rbind.data.frame(ssr_dat_pt3b) %>% 
  rbind.data.frame(ssr_dat_pt3d) %>% 
  rbind.data.frame(ssr_dat_pt3e)
head(ssr_dat_combined)


#sample number should be the same!
ssr_dat %>% select(sample_id) %>% 
  filter(str_detect(sample_id,"-dup")==FALSE) %>% 
  distinct() %>% nrow() #3104
ssr_dat_combined %>% select(sample_id) %>% 
  filter(str_detect(sample_id,"-dup")==FALSE) %>%
  distinct() %>% nrow() #3104

write_rds(ssr_dat_combined,"r_data/ssr_dat.rds")




