######### branch-specific differential expression analysis

#need Monrun object and full data matrix (all mets) from metabolomics_monocle_sex_specific.R
#start with female trajectory
MonRun <- readRDS(file="data_objects/MonRunF_mets.RDS")
temp <- readRDS(file="data_objects/metsF_matrix_all.rds")
gene_short_name <- rownames(temp)
table(MonRun$State2)
#there are 8 states in the female tree
#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))
l2$p_7 <- rep(0,length(gene_short_name))
l2$p_8 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))
l2$d_7 <- rep(0,length(gene_short_name))
l2$d_8 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(MonRun$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$p_3[i] <- tk$s[2,4]
  l2$p_4[i] <- tk$s[3,4]
  l2$p_5[i] <- tk$s[4,4]
  l2$p_6[i] <- tk$s[5,4]
  l2$p_7[i] <- tk$s[6,4]
  l2$p_8[i] <- tk$s[7,4]
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
  l2$d_7[i] <- tk$s[6,1]
  l2$d_8[i] <- tk$s[7,1]
}

#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)
names(df3)[names(df3) == 'gene_names'] <- 'metabolite'

write.csv(df3,file='data_objects/FEMALE_mets_DEanova_stats.csv',quote=F,row.names=F)
file <- synapser::File(path='data_objects/FEMALE_mets_DEanova_stats.csv', parentId='syn45147359')
file <- synapser::synStore(file)


#Male trajectory
MonRun <- readRDS(file="data_objects/MonRunM_mets.RDS")
temp <- readRDS(file="data_objects/metsM_matrix_all.rds")
gene_short_name <- rownames(temp)
table(MonRun$State2)
#there are 6 states in the male tree
#pre-process data for ANOVA test
l2 <- list()
l2$gene_names <- gene_short_name
l2$p_2 <- rep(0,length(gene_short_name))
l2$p_3 <- rep(0,length(gene_short_name))
l2$p_4 <- rep(0,length(gene_short_name))
l2$p_5 <- rep(0,length(gene_short_name))
l2$p_6 <- rep(0,length(gene_short_name))

l2$d_2 <- rep(0,length(gene_short_name))
l2$d_3 <- rep(0,length(gene_short_name))
l2$d_4 <- rep(0,length(gene_short_name))
l2$d_5 <- rep(0,length(gene_short_name))
l2$d_6 <- rep(0,length(gene_short_name))

for (i in 1:length(gene_short_name)){
  l <- list()
  l$x <- as.vector(t(temp[i,]))
  l$s <- as.character(MonRun$State)
  
  df <- as.data.frame(l)
  res.aov <- aov(x ~ s, data = df)
  tk <- TukeyHSD(res.aov)
  
  l2$p_2[i] <- tk$s[1,4]
  l2$p_3[i] <- tk$s[2,4]
  l2$p_4[i] <- tk$s[3,4]
  l2$p_5[i] <- tk$s[4,4]
  l2$p_6[i] <- tk$s[5,4]
  
  l2$d_2[i] <- tk$s[1,1]
  l2$d_3[i] <- tk$s[2,1]
  l2$d_4[i] <- tk$s[3,1]
  l2$d_5[i] <- tk$s[4,1]
  l2$d_6[i] <- tk$s[5,1]
}

#save the data
df2 <- as.data.frame(l2)
dfa <- dplyr::select(df2,gene_names,dplyr::starts_with('p'))
dfb <- dplyr::select(df2,gene_names,dplyr::starts_with('d'))

dfa1 <- tidyr::gather(dfa,'state','pvalue',-gene_names)
dfb1 <- tidyr::gather(dfb,'state','effect',-gene_names)

dfa1$state <- sapply(dfa1$state,function(x) strsplit(x,'p_')[[1]][2])
dfb1$state <- sapply(dfb1$state,function(x) strsplit(x,'d_')[[1]][2])

df3 <- dplyr::left_join(dfa1,dfb1)
names(df3)[names(df3) == 'gene_names'] <- 'metabolite'

write.csv(df3,file='data_objects/MALE_mets_DEanova_stats.csv',quote=F,row.names=F)
file <- synapser::File(path='data_objects/MALE_mets_DEanova_stats.csv', parentId='syn45147359')
file <- synapser::synStore(file)


