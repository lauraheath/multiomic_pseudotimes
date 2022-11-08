install.packages('ggpubr')
install.packages('VennDiagram')
library(VennDiagram)
library(gridExtra)
synapser::synLogin()

#set working directory and create folders to stash data objects and figures
setwd("~/multiomic_pseudotimes/")
dir.create("data_objects")
dir.create("figures")


#upload proteomics pseudotimes & states by sex:
p1 <- synapser::synGet('syn35317137')
female_prots <- read.csv(p1$path)

p1 <- synapser::synGet('syn35317790')
male_prots <- read.csv(p1$path)


#change column names to differentiate from RNAseq
names(female_prots)[names(female_prots) == 'Pseudotime'] <- 'Pseudotime_prot'
names(female_prots)[names(female_prots) == 'State2'] <- 'State_prot'
names(female_prots)[names(female_prots) == 'SampleID'] <- 'SampleID_prot'
names(female_prots)[names(female_prots) == 'pseudotime_sc'] <- 'pseudotime_sc_prot'
female_prots$msex[female_prots$msex == 0] <- 'female'
female_prots$diagnosis[female_prots$diagnosis == 'control'] <- 'CT'
female_prots$diagnosis[female_prots$diagnosis == 'other'] <- 'OTHER'
female_prots$batch<-NULL


names(male_prots)[names(male_prots) == 'Pseudotime'] <- 'Pseudotime_prot'
names(male_prots)[names(male_prots) == 'State2'] <- 'State_prot'
names(male_prots)[names(male_prots) == 'SampleID'] <- 'SampleID_prot'
names(male_prots)[names(male_prots) == 'pseudotime_sc'] <- 'pseudotime_sc_prot'
male_prots$msex[male_prots$msex == 1] <- 'male'
male_prots$diagnosis[male_prots$diagnosis == 'control'] <- 'CT'
male_prots$diagnosis[male_prots$diagnosis == 'other'] <- 'OTHER'
male_prots$batch<-NULL

#combine 
allprots <- rbind(female_prots, male_prots)


#upload rnaseq pseudotimes & states
p2 <- synapser::synGet('syn38354143')
female_rnaseq <- read.csv(p2$path)


p2 <- synapser::synGet('syn38348029')
male_rnaseq <- read.csv(p2$path)


#female_rnaseq <- subset(female_rnaseq, select=c(individualID, SampleID, State, Pseudotime, pseudotime_sc))
names(female_rnaseq)[names(female_rnaseq) == 'Pseudotime'] <- 'Pseudotime_rnaseq'
names(female_rnaseq)[names(female_rnaseq) == 'State'] <- 'State_rnaseq'
names(female_rnaseq)[names(female_rnaseq) == 'SampleID'] <- 'SampleID_rnaseq'
names(female_rnaseq)[names(female_rnaseq) == 'pseudotime_sc'] <- 'pseudotime_sc_rnaseq'

names(male_rnaseq)[names(male_rnaseq) == 'Pseudotime'] <- 'Pseudotime_rnaseq'
names(male_rnaseq)[names(male_rnaseq) == 'State'] <- 'State_rnaseq'
names(male_rnaseq)[names(male_rnaseq) == 'SampleID'] <- 'SampleID_rnaseq'
names(male_rnaseq)[names(male_rnaseq) == 'pseudotime_sc'] <- 'pseudotime_sc_rnaseq'

#combine:
allRNAseq <- rbind(female_rnaseq, male_rnaseq)

#merge proteomics with rnaseq pseudotimes plus clinical metadata
pstimes <- dplyr::full_join(allprots, allRNAseq)

#create a flag for whether any given patient has a sample from each datatype 
pstimes$rna_samples = 0
pstimes$rna_samples[!is.na(pstimes$SampleID_rnaseq)] <- 1
pstimes$prot_samples = 0
pstimes$prot_samples[!is.na(pstimes$Pseudotime_prot)] <- 1

table(pstimes$prot_samples, pstimes$rna_samples)

#save to the workspace:
write.csv(pstimes, file="data_objects/Prot_and_RNA_pseudotimes.csv")
file <- synapser::File(path='data_objects/Prot_and_RNA_pseudotimes.csv', parentId='syn44292253')
file <- synapser::synStore(file)


#get gene lists that were used to create the monocle objects:
#differentially expressed proteins (not stratified by sex--sample sizes too small)
p3 <- synapser::synGet('syn35221005')
DEproteins <- read.csv(p3$path)
#keep proteins with unadjusted p value < 0.05
DEproteins <- subset(DEproteins, DEproteins$PVal<0.05)
dim(DEproteins)

#differentially expressed genes from AMP-AD 2.0 RNAseq harmonization project, stratified by sex and tissue (use DLPFC)
de_file <- synapser::synGet('syn26967458')
de1 <- read.delim(de_file$path)
#separate by sex and make short gene names unique
de_male <- dplyr::filter(de1,Comparison=='AD_male_DLPFC - CT_male_DLPFC')
de_male$hgnc_symbol<-make.unique(de_male$hgnc_symbol)
de_female <- dplyr::filter(de1,Comparison=='AD_female_DLPFC - CT_female_DLPFC')
de_female$hgnc_symbol<-make.unique(de_female$hgnc_symbol)

#pull all genes with FDR p-value<0.1
InF <- which(de_female$adj.P.Val<0.1)
DEFemaleGenes <- de_female[InF,]
InM <- which(de_male$adj.P.Val<0.1)
DEMaleGenes <- de_male[InM,]


#venn diagram of rnaseq genes and prot genes
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

#create venn diagram showing overlap of genes and proteins in female data sets
rna <- DEFemaleGenes$hgnc_symbol
prot <- DEproteins$GeneName


venn <- venn.diagram(
  x=list(rna, prot),
  category.names=c("RNAseq Genes","Proteins"),
  filename='figures/FEMALE_gene_prot_overlap.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer",
  cat.pos=c(-20,20),
  cat.dist=c(0.05,0.05)
)


#create venn diagram showing overlap of genes and proteins in male data sets
rna <- DEMaleGenes$hgnc_symbol
prot <- DEproteins$GeneName


venn <- venn.diagram(
  x=list(rna, prot),
  category.names=c("RNAseq Genes","Proteins"),
  filename='figures/MALE_gene_prot_overlap.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer",
  cat.pos=c(-20,20),
  cat.dist=c(0.05,0.05)
)






#venn diagram of which patients overlap between rnaseq and protein studies
rna_ps <- subset(pstimes, pstimes$rna_samples==1)
prot_ps <- subset(pstimes, pstimes$prot_samples==1)
rna_patients <- rna_ps$individualID
prot_patients <- prot_ps$individualID

venn.diagram(
  x=list(rna_patients, prot_patients),
  category.names=c("RNAseq IDs","Proteomics IDs"),
  filename='figures/Patient_overlap_venn.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer",
  cat.pos=c(-24,18),
  cat.dist=c(0.05,0.05)
)


#venn diagram of which patients overlap between rnaseq and protein studies BY SEX
rna_ps <- subset(pstimes, pstimes$rna_samples==1 & pstimes$msex=='female')
prot_ps <- subset(pstimes, pstimes$prot_samples==1 & pstimes$msex=='female')
rna_patients <- rna_ps$individualID
prot_patients <- prot_ps$individualID

venn.diagram(
  x=list(rna_patients, prot_patients),
  category.names=c("RNAseq IDs","Proteomics IDs"),
  filename='figures/FEMALE_Patient_overlap_venn.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer",
  cat.pos=c(-24,18),
  cat.dist=c(0.05,0.05)
)

#now the males
rna_ps <- subset(pstimes, pstimes$rna_samples==1 & pstimes$msex=='male')
prot_ps <- subset(pstimes, pstimes$prot_samples==1 & pstimes$msex=='male')
rna_patients <- rna_ps$individualID
prot_patients <- prot_ps$individualID

venn.diagram(
  x=list(rna_patients, prot_patients),
  category.names=c("RNAseq IDs","Proteomics IDs"),
  filename='figures/MALE_Patient_overlap_venn.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer",
  cat.pos=c(-24,18),
  cat.dist=c(0.05,0.05)
)

#run correlations between proteomics pseudotimes and rnaseq pseudotimes (use scaled pseudotimes)
cor(pstimes$pseudotime_sc_rnaseq, pstimes$pseudotime_sc_prot, method="pearson", use="pairwise.complete.obs")
cor(pstimes$pseudotime_sc_rnaseq, pstimes$pseudotime_sc_prot, method="spearman", use="pairwise.complete.obs")


plot(pstimes$pseudotime_sc_rnaseq, pstimes$pseudotime_sc_prot)
abline(lm(pstimes$pseudotime_sc_rnaseq~pstimes$pseudotime_sc_prot))



#The following monocle objects were calculated in the scripts UPDATE_rnaseq_lineage.r 
# (in AMP-AD-2.0-transcriptomics-lineage-update repo) and prot_lineage_monocle_rerun.R 
# (in prot-lineage repo), and stored in the shared working space

#female monocle objects (P for proteomics, R for RNAseq):
p4 <- synapser::synGet('syn44293198')
MonRun_FP <- readRDS(p4$path)

p5 <- synapser::synGet('syn44293242')
MonRun_FR <- readRDS(p5$path)

p6 <- synapser::synGet('syn44293269')
MonRun_MP <- readRDS(p6$path)

p7 <- synapser::synGet('syn44293357')
MonRun_MR <- readRDS(p7$path)



#add rna-seq pseudotimes to the monocle object
rnaseq_pstime <- subset(pstimes, select=c(SampleID_rnaseq, pseudotime_sc_rnaseq, pseudotime_sc_prot, rna_samples))
rnaseq_pstime <- subset(rnaseq_pstime, rnaseq_pstime$rna_samples==1)
rnaseq_pstime$rna_samples<-NULL
names(rnaseq_pstime)[names(rnaseq_pstime) == 'SampleID_rnaseq'] <- 'specimenID'



head(pData(F_MonRunT))
tail(pData(F_MonRunT))

#reorder rnaseq_pstime data frame to match the order of the specimenIDs in the monocle CDS:
sampleIDs <- as.data.frame(pData(F_MonRunT))
rnaseq_pstime <- rnaseq_pstime[order(match(rnaseq_pstime$specimenID, sampleIDs$specimenID)),]

F_MonRunT$rnaseq_pseudotime <- rnaseq_pstime$pseudotime_sc_rnaseq
F_MonRunT$prot_pseudotime <- rnaseq_pstime$pseudotime_sc_proteomics

#tiff(file='~/prot-lineage/figures/FEMALE_tree_rnaseq_pstime.tiff',height=150,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/figures/MALE_tree_rnaseq_pstime.tiff',height=150,width=100,units='mm',res=300)
g <- plot_cell_trajectory(F_MonRunT,color_by = "prot_pseudotime",show_branch_points=F,use_color_gradient = F,cell_size = 1) + scale_color_gradient (low="blue", high="red")
g
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(pstimes_combined, x="pseudotime_sc_proteomics", y="pseudotime_sc_rnaseq", 
                  add="reg.line",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_diagnosis.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_diagnosis.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(pstimes_combined, x="pseudotime_sc_proteomics", y="pseudotime_sc_rnaseq", color = "diagnosis",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

pstimes_combined$braaksc <- as.factor(pstimes_combined$braaksc)
#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_braak.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_braak.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(pstimes_combined, x="pseudotime_sc_proteomics", y="pseudotime_sc_rnaseq", color = "braaksc",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_cerad.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_cerad.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", color = "ceradsc",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr_cogdx.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr_cogdx.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="pseudotime_sc", y="rnaseq_pseudotime_sc", color = "cogdx",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()





#compare DE genes included as feature sets in bulk rnaseq and prot


inters <- intersect(rnaseq_genesF$gene_short_name,protgenes$GeneName)
diffs <- setdiff(rnaseq_genesF$gene_short_name,protgenes$GeneName)

#rerun lineage analysis on both rna-seq and proteomics data using only the shared genes as a feature set
inters <- as.data.frame(inters)
names(inters)[names(inters) == "inters"] <- "gene_short_name"
Log2_Normalized <- readRDS(file="~/prot-lineage/data/Log2_Normalized.rds")
Meta <- readRDS(file="~/prot-lineage/data/Meta.rds")
Log2_Normalized2 <- Log2_Normalized
Log2_Normalized2$proteins <- rownames(Log2_Normalized2)
Log2_Normalized2$proteins <- gsub("\\|.*", "", Log2_Normalized2$proteins)

#function to run monocle analysis
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}

#run get_proteins_script_JG.r to get metadata and log2 normalized protein matrix
synapser::synLogin()

Dat <- Log2_Normalized2
Dat[is.na(Dat)] <- 0

#select only the rows with gene short names that match the intersecting gene short names (some are repeated peptides, same gene)
genes2<-c()
for (gene in unique(c(as.vector(inters$gene_short_name)))){
  if (gene %in% Dat$proteins){
    genes2 <- c(genes2,which(Dat$proteins==gene))
  }
}
length(genes2)
Dat2 <- Dat[genes2,]
dim(Dat2)
Dat2$proteins<-NULL


#Keeping only female data (msex==0 is female, msex==1 is male; run separately for sex-specific analysis)
In_S <- which(Meta$msex == 0)
#In_S <- which(Meta$msex == 1)
Dat2 <- Dat2[,In_S]
Meta2 <- Meta[In_S,]

gene_short_name <- rownames(Dat2)
temp <- Dat2
temp2 <- Meta2


temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc <- factor(temp2$ceradsc,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx,levels = c(1:6))

rownames(temp)<-NULL
rownames(temp2)<-NULL

MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

x <- list()
x$SampleID <- MonRun$batchChannel
x$rnaseqID <- MonRun$rnaseq_id
x$State2 <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$diagnosis <- MonRun$diagnosis
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$apoe <- MonRun$APO
x$educ   <- MonRun$educ
x$pmi <- MonRun$pmi
x$batch <- MonRun$batch
x$mmse <- MonRun$cts_mmse30_lv
x$age_death <- MonRun$age_death
x$rna_seq_sample <- MonRun$rnaseq
x$SampleID <- as.character(x$SampleID)
F_inters_prots <- as.data.frame(x)
F_inters_prots$pseudotime_sc <- scale(F_inters_prots$Pseudotime, center=F)




#### now upload rna-seq data and rerun lineage analysis on intersecting proteins


#function to run monocle analysis
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}

dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- read.delim(dlpfcCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file (rosmap covariates): syn8466814
dlpfcCovObj <- synapser::synGet('syn11024258')
covars <- read.delim(dlpfcCovObj$path,stringsAsFactors = F)

#converting ENSG to gene symbols
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           #host='useast.ensembl.org')
                           host='uswest.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}
Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}
Dat$gene_short_name <- Make.Gene.Symb(Dat$ensembl_gene_id)

#select only the rows with gene short names that match the intersecting gene short names (some are repeated peptides, same gene)
genes2<-c()
for (gene in unique(c(as.vector(inters$gene_short_name)))){
  if (gene %in% Dat$gene_short_name){
    genes2 <- c(genes2,which(Dat$gene_short_name==gene))
  }
}
length(genes2)
Dat2 <- Dat[genes2,]
dim(Dat2)


Names <- colnames(Dat2)

for (i in 1:length(Names)){
  
  Names[i] <- substring(Names[i],2)
  
}


colnames(Dat2) <- Names
cNames <- covars$SampleID
l <- length(Names)

#deleting columns not in the covariate list
temp <- rep(T,l)
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    temp[i] <- F
  }
}

In <- which(temp)
#print(temp)
Dat2 <- Dat2[,In]

#deleting extra rows in covariate list
Names <- Names[In]
l <- length(cNames)
temp <- rep(T,l)
for (i in 1:l){
  if (!(cNames[i] %in% Names)){
    temp[i] <- F
  }
}
In <- which(temp)
covars <- covars[In,]

ColNorm <- function(Dat2){
  
  M = max(colSums(Dat2))
  l <- length(colnames(Dat2))
  
  for( i in 1:l){
    
    Dat2[,i] = Dat2[,i]*(M/sum(Dat2[,i]))
    
  }
  
  return(Dat2)
}

DatNorm <- ColNorm(Dat2)

#removing bad batches
DatNorm <- DatNorm[,covars$Batch<7]
covars <- covars[covars$Batch<7,] 


#Keeping only female data 
#Sex <- 'FEMALE'
In_S <- which(covars$msex == 0)
DatNorm2 <- DatNorm[,In_S]
covars2 <- covars[In_S,]

temp <- DatNorm2
temp2 <- covars2

rosmapObj <- synapser::synGet('syn3191087')
rosmap <- data.table::fread(rosmapObj$path,data.table=F)

#add in braak score & cerad score
rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)
rosmapRNAid<-dplyr::left_join(rosmapId,rosmap)
#remove duplicate rows
rosmapRNAid <- unique(rosmapRNAid)
rosmapRNAid2 <- subset(rosmapRNAid, select=c(rnaseq_id,braaksc,ceradsc))
names(rosmapRNAid2)[names(rosmapRNAid2) == "rnaseq_id"] <- "SampleID"

temp2<-dplyr::left_join(temp2,rosmapRNAid2, by="SampleID")

temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc <- factor(temp2$ceradsc,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx, levels = c(1:6))


gene_short_name <- inters$gene_short_name

rownames(temp)<-NULL
rownames(temp2)<-NULL

MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g

x <- list()
x$rnaseqID <- MonRun$SampleID
x$State_rnaseq <- MonRun$State
x$Pseudotime_rnaseq <- MonRun$Pseudotime
F_inters_rnaseq <- as.data.frame(x)
F_inters_rnaseq$pseudotime_scRNA <- scale(F_inters_rnaseq$Pseudotime, center=F)


F_intersecting_pseudotimes <- merge(F_inters_prots, F_inters_rnaseq, by="rnaseqID")

cor(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA, method="pearson")
cor(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA, method="spearman")
plot(F_intersecting_pseudotimes$pseudotime_sc, F_intersecting_pseudotimes$pseudotime_scRNA)
abline(lm(F_intersecting_pseudotimes$pseudotime_sc~F_intersecting_pseudotimes$pseudotime_scRNA))
lines(lowess(corrs$rnaseq_pseudotime_sc,corrs$pseudotime_sc))

#tiff(file='~/prot-lineage/figures/FEMALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_rnaseq_prot_corr.tiff',height=85,width=100,units='mm',res=300)
ggpubr::ggscatter(corrs, x="Pseudotime", y="rnaseq_Pseudotime", 
                  add="reg.line",
                  cor.coef=TRUE, cor.method="pearson",
                  xlab="Proteomics Pseudotime", ylab="RNA-seq Pseudotime")
dev.off()
