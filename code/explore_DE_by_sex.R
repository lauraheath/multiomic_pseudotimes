
setwd("~/multiomic_pseudotimes/")

library(pheatmap)

#load differential expression results (AD case vs control, by sex and tissue)
de_file <- synapser::synGet('syn26967458')
de1 <- read.delim(de_file$path)
de_male <- dplyr::filter(de1,Comparison=='AD_male_DLPFC - CT_male_DLPFC')
#make short gene names unique
de_male$hgnc_symbol<-make.unique(de_male$hgnc_symbol)
de_female <- dplyr::filter(de1,Comparison=='AD_female_DLPFC - CT_female_DLPFC')
de_female$hgnc_symbol<-make.unique(de_female$hgnc_symbol)


Fgenes <- subset(de_female, select=c(hgnc_symbol, logFC, adj.P.Val, Direction))
names(Fgenes)[names(Fgenes) == 'logFC'] <- 'logFC_female'
names(Fgenes)[names(Fgenes) == 'adj.P.Val'] <- 'adj.P.Val_female'
names(Fgenes)[names(Fgenes) == 'Direction'] <- 'Direction_female'
Mgenes <- subset(de_male, select=c(hgnc_symbol, logFC, adj.P.Val, Direction))
names(Mgenes)[names(Mgenes) == 'logFC'] <- 'logFC_male'
names(Mgenes)[names(Mgenes) == 'adj.P.Val'] <- 'adj.P.Val_male'
names(Mgenes)[names(Mgenes) == 'Direction'] <- 'Direction_male'
allgenes <- left_join(Fgenes, Mgenes)

plot(allgenes$logFC_female, allgenes$logFC_male)
abline(lm(allgenes$logFC_female~allgenes$logFC_male))

#select significant genes only
allgenes2 <- subset(allgenes, allgenes$adj.P.Val_female<0.1 | allgenes$adj.P.Val_male<0.1)
#allgenes2 <- subset(allgenes, allgenes$adj.P.Val_female<0.1)
plot(allgenes2$logFC_female, allgenes2$logFC_male)
abline(lm(allgenes2$logFC_female~allgenes2$logFC_male))
cor(allgenes2$logFC_female, allgenes2$logFC_male, method="pearson", use="pairwise.complete.obs")



#venn diagram of which patients overlap between direction of change (up, down, or none)
# female_up <- subset(allgenes2, allgenes2$Direction_female=='up')
# male_up <- subset(allgenes2, allgenes2$Direction_male=='up')
# female_down <- subset(allgenes2, allgenes2$Direction_female=='down')
# male_down <- subset(allgenes2, allgenes2$Direction_male=='down')
# female_none <- subset(allgenes2, allgenes2$Direction_female=='none')
# male_none <- subset(allgenes2, allgenes2$Direction_male=='none')
# upgenesF <- female_up$hgnc_symbol
# upgenesM <- male_up$hgnc_symbol
# downgenesF <- female_down$hgnc_symbol
# downgenesM <- male_down$hgnc_symbol
# nonegenesF <- female_none$hgnc_symbol
# nonegenesM <- male_none$hgnc_symbol

venn.diagram(
  x=list(upgenesF, upgenesM, downgenesF, downgenesM),
  category.names=c("Female Up","Male Up", "Female down", "Male down"),
  filename='figures/Gene direction overlap.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer"
  #cat.pos=c(-24,18,135),
  #cat.dist=c(0.05,0.05,0.05)
)


#overlap in significant genes
female_topgenes <- subset(allgenes2, allgenes2$adj.P.Val_female<0.1)
male_topgenes <- subset(allgenes2, allgenes2$adj.P.Val_male<0.1)
female_genes <- female_topgenes$hgnc_symbol
male_genes <- male_topgenes$hgnc_symbol
venn.diagram(
  x=list(female_genes, male_genes),
  category.names=c("Female Top Genes","Male Top Genes"),
  filename='figures/Gene by sex overlap.png',
  output=TRUE,
  imagetype="png",
  height=1500,
  width=2000,
  resolution=300,
  compression="lzw",
  lwd = 2,
  #col=c("#440154ff", '#21908dff', '#fde725ff'),
  #fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 2,
  fontface = "bold",
  cat.cex = 2,
  cat.default.post="outer"
  #cat.pos=c(-24,18,135),
  #cat.dist=c(0.05,0.05,0.05)
)
