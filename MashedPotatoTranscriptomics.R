#DE as FACTOR
# C elegans
setwd("~/Desktop/Honors/featureCounts/elegans/")
library(DESeq2)
library(tidyverse)

elegans <- read.table("counts_elegans", header=T)
elegans <- elegans[,-c(2:6)]
names(elegans) <- c("WBGENE","SRR11996348", "SRR11996349", "SRR11996350","SRR11996351","SRR11996352","SRR11996353","SRR11996354","SRR11996355","SRR11996356")
elegans$WBGENE <- str_remove(elegans$WBGENE, "^Gene:")
m <- as.matrix(elegans[,-1])
mode(m) <- "integer"
rownames(m) <- elegans$WBGENE
meta <- data.frame(colnames(m), Group = rep(c("Control","5ug_BaP","a20ug_BaP"), each=3)) #these are factors, numeric for an lm
dds <- DESeqDataSetFromMatrix(countData=m, 
                              colData=meta, 
                              design=~Group)
counts(dds)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
ddselegans1 <- dds[keep,]
elegansde1 <- DESeq(ddselegans1) #object used for results
res5 <- results(elegansde1, alpha = 0.05, lfcThreshold = 1, contrast = c("group","5ug_BaP","Control"))
res20 <- results(elegansde1, alpha = 0.05, lfcThreshold = 1, contrast = c("group","20ug_BaP","Control"))

summary(res5)
res5df <- as.data.frame(res5)
res5df1 <- subset(res5df, log2FoldChange > 1 & padj < 0.05) #just to remember where DEGs came from
rownames(res5df1)
summary(res20)
res20df <- as.data.frame(res20)
res20df1 <- subset(res20df, log2FoldChange > 1 & padj < 0.05)
res20df2 <- subset(res20df, log2FoldChange < -1 & padj < 0.05)
res20df3 <- rbind(res20df1, res20df2)
write.csv(res20df3, "Allsiggenesres20.csv")
view(res20df1)
view(res20df2)

sum(rownames(res20df) %in% rownames(res5df1)) # all in 5ug are in 20ug
table(rownames(res20df)[rownames(res20df) %in% rownames(res5df1)])

#worthless plots
#plotMA(res5)
#plotMA(res20)

vsd <- vst(ddselegans1)
?plotPCA
pcaData_elegans <- plotPCA(vsd, "Group", returnData = T)
percentVarEle <- round(100 * attr(pcaData_elegans, "percentVar"))
hull.data <- do.call(rbind, lapply(split(pcaData_elegans, list(pcaData_elegans$Group)), 
                                   function(x) x[chull(x[,c("PC1", "PC2")]),]))
ggplot(pcaData_elegans, aes(x = PC1, y= PC2, color=Group)) +
  geom_polygon(data = hull.data, alpha = 0) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarEle[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarEle[2],"% variance")) + 
  theme_classic() 


# plotPCA(vsd, "group", ntop = 50, pcsToUse = 1:2)
# plotPCA(vsd, "group", ntop = 50)

?plotPCA
# volcano plot
# library(EnhancedVolcano)
#EnhancedVolcano(res5,lab = "", x = 'log2FoldChange', y='padj', ylim = c(0,50))

#EnhancedVolcano(res20,lab = "", x = 'log2FoldChange', y='padj', ylim = c(0,50))

#DE as NUMERIC
meta_num <- data.frame(colnames(m), group = rep(c(0,5,20), each=3))
dds12 <- DESeqDataSetFromMatrix(countData=m,
                                colData=meta_num,
                                design=~group)
smallestGroupSizeelenum <- 3
keepelenum <- rowSums(counts(dds12) >= 10) >= smallestGroupSizeelenum
ddselegansnum <- dds12[keepelenum,]
elegansde2 <- DESeq(ddselegansnum)
resnum <- results(elegansde2)
resnumalpha <- subset(resnum, padj < 0.05)
resnumalphadf <- as.data.frame(resnumalpha)
resnumalphadfsub <- subset(resnumalphadf, log2FoldChange > 0)
resnumalphadfsubneg <- subset(resnumalphadf, log2FoldChange < 0)

write.csv(resnumalphadf, "resnumalphadf.csv")
view(resnumalphadf)
#rownames(resnumdf) %in% rownames(res20df) # all in 20ug are in num
#abdsig2 <- resnumdf %>% filter(rownames(resnumdf) %in% rownames(res20df))
#abdsig <- resnumdf %>% filter(rownames(resnumdf) %in% rownames(res5df1)) #THIS ONE
#Since all overlapped from 5 - 20 uM, drawing a conclusion is valid, but what is the conclusion?
#maybe just a table with FCs, or presenting a range of FC coefficients?

library(DEXSeq)
library(tidyverse)
setwd("~/Desktop/Honors/featureCounts/elegans/")
dexelegans <- read.table("dexseq_reads_elegans", header=T)
dexelegans <- dexelegans[,-c(2:6)]
names(dexelegans) <- c("WBGENE","SRR11996348", "SRR11996349", "SRR11996350","SRR11996351","SRR11996352","SRR11996353","SRR11996354","SRR11996355","SRR11996356")

dexm <- as.matrix(dexelegans[,-1])
mode(dexm) <- "integer"
rownames(dexm) <- dexelegans$WBGENE

sample.data.elegans <- read.table("transcript2gene.txt")
names(sample.data.elegans) <- c("WBGene","TranscriptID")
sample.data.elegans$TranscriptID <- str_remove(sample.data.elegans$TranscriptID, "^Name=")
sample.data.elegans$WBGene <- str_remove(sample.data.elegans$WBGene, "^Parent=Gene:")

dexm_16 <- dexm[,c(1:6)]
meta_16 <- meta[c(1:6),]
dxd05 <- DEXSeqDataSet(countData = dexm_16, #featureCounts output
                       sampleData = meta_16, #same as de
                       design=~sample + exon + group:exon, #but it's actually transcript
                       featureID = sample.data.elegans$TranscriptID, #the transcripts
                       groupID = sample.data.elegans$WBGene) #genes they belong to
system.time({
  dxd05 <- estimateSizeFactors(dxd05)
  dxd05 <- estimateDispersions(dxd05, quiet = T)
  dxd05 <- testForDEU(dxd05, reducedModel =~sample + exon) # also try coding as a number, getting a fc/unit 
}) # stratify model, only including 1 treatment and control

dx05red <- DEXSeqResults(dxd05, independentFiltering = F)
qval05 <- perGeneQValue(dx05red)
dxr.g05 <- data.frame(gene=names(qval05),qval05)
subset(dxr.g05, qval05 < 0.1)

dexm_020 <- dexm[,c(1:3, 7:9)]
meta_020 <- meta[c(1:3, 7:9),]
dxd020 <- DEXSeqDataSet(countData = dexm_020, #featureCounts output
                        sampleData = meta_020, #same as de
                        design=~sample + exon + group:exon, #but it's actually transcript
                        featureID = sample.data.elegans$TranscriptID, #the transcripts
                        groupID = sample.data.elegans$WBGene) #genes they belong to
system.time({
  dxd020 <- estimateSizeFactors(dxd020)
  dxd020 <- estimateDispersions(dxd020, quiet = T)
  dxd020 <- testForDEU(dxd020, reducedModel =~sample + exon) # also try coding as a number, getting a fc/unit 
}) # stratify model, only including 1 treatment and control

dx020red <- DEXSeqResults(dxd020, independentFiltering = F)
qval020 <- perGeneQValue(dx020red)
dxr.g020 <- data.frame(gene=names(qval020),qval020)
subset(dxr.g020, qval020 < 0.1)


# AS as num - nothing, same as above

meta_num <- data.frame(colnames(m), group = rep(c(0,5,20), each=3))
dxd_num <- DEXSeqDataSet(countData=dexm,
                         sampleData=meta_num,
                         design=~sample + exon + group:exon,
                         featureID = sample.data.elegans$TranscriptID,
                         groupID = sample.data.elegans$WBGene)

system.time({
  dxd_num <- estimateSizeFactors(dxd_num)
  dxd_num <- estimateDispersions(dxd_num, quiet = T)
  dxd_num <- testForDEU(dxd_num, reducedModel =~sample + exon) # also try coding as a number, getting a fc/unit 
}) # stratify model, only including 1 treatment and control

dx_numred <- DEXSeqResults(dxd_num, independentFiltering = F)
qvalnum <- perGeneQValue(dx_numred)
dxr.gnum <- data.frame(gene=names(qvalnum),qvalnum)
dxr.gnum_sig <- subset(dxr.gnum, qvalnum < 0.05)
write.csv(dxr.gnum_sig, "elegansDEXseqNUM.csv")

library(DEXSeq)
library(tidyverse)
setwd("~/Desktop/Honors/featureCounts/Elegans/")
exonelegans <- read.table("exoncountselegans2.txt", header=T)
exonelegans <- exonelegans[,-c(2:6)]
names(exonelegans) <- c("GeneID","SRR11996348", "SRR11996349", "SRR11996350","SRR11996351","SRR11996352","SRR11996353","SRR11996354","SRR11996355","SRR11996356")

exonm <- as.matrix(exonelegans[,-1])
mode(exonm) <- "integer"
rownames(exonm) <- exonelegans$GeneID
exons <- as.data.frame(exonm)
exons$WBGene <- gsub("\\..*$", "", rownames(exons))
exons$exonname <- rownames(exons)

exonm05 <- exonm[,c(1:6)]
# meta_16
dxdexon05 <- DEXSeqDataSet(countData = exonm05, #featureCounts output
                           sampleData = meta_16, #same as de
                           design=~sample + exon + Group:exon,
                           featureID = exons$exonname, #the exons
                           groupID = exons$WBGene) #genes they belong to

system.time({
  dxdexon05 <- estimateSizeFactors(dxdexon05)
  dxdexon05 <- estimateDispersions(dxdexon05, quiet = T)
  dxdexon05 <- testForDEU(dxdexon05, reducedModel =~sample + exon) 
}) 

dxexonred05 <- DEXSeqResults(dxdexon05, independentFiltering = F)

qvalex05 <- perGeneQValue(dxexonred05)
dxr.gex05 <- data.frame(gene=names(qvalex05),qvalex05)
dxrgex05sig <- subset(dxr.gex05, qvalex05 < 0.05) #19 sig genes
write.csv(dxrgex05sig, "elegansexon_05.csv")

exonm020 <- exonm[,c(1:3,7:9)]
# meta_020
dxdexon020 <- DEXSeqDataSet(countData = exonm020, #featureCounts output
                            sampleData = meta_020, #same as de
                            design=~sample + exon + group:exon,
                            featureID = exons$exonname, #the exons
                            groupID = exons$WBGene) #genes they belong to

system.time({
  dxdexon020 <- estimateSizeFactors(dxdexon020)
  dxdexon020 <- estimateDispersions(dxdexon020, quiet = T)
  dxdexon020 <- testForDEU(dxdexon020, reducedModel =~sample + exon) 
}) 

dxexonred020 <- DEXSeqResults(dxdexon020, independentFiltering = F)

qvalex020 <- perGeneQValue(dxexonred020)
dxr.gex020 <- data.frame(gene=names(qvalex020),qvalex020)
dxrgex020sig <- subset(dxr.gex020, qvalex020 < 0.05) #29 sig genes
write.csv(dxrgex020sig, "elegansexon_020.csv")
#DEU as NUMERIC
# meta_num is good for the sampleData
ddsexonnum <- DEXSeqDataSet(countData=exonm,
                            sampleData = meta_num,
                            design=~sample + exon + group:exon,
                            featureID = exons$exonname,
                            groupID = exons$WBGene)
system.time({
  ddsexonnum <- estimateSizeFactors(ddsexonnum)
  ddsexonnum <- estimateDispersions(ddsexonnum, quiet = T)
  ddsexonnum <- testForDEU(ddsexonnum, reducedModel =~sample + exon) 
}) 
ddsexonnumres <- DEXSeqResults(ddsexonnum, independentFiltering = F)

qvalexnum <- perGeneQValue(ddsexonnumres)
dxr.gexnum <- data.frame(gene=names(qvalexnum),qvalexnum)
dxrgexnumsig <- subset(dxr.gexnum, qvalexnum < 0.05)
write.csv(dxrgexnumsig, "elegansexon_NUM.csv")

#DE as FACTOR
# G morhua
setwd("~/Desktop/Honors/featureCounts/Elegans/")
library(DESeq2)
library(tidyverse)

morhua <- read.table("counts_cod_v1", header=T)
morhua <- morhua[,-c(2:6)]
names(morhua) <- c("Geneid","SRR6296791","SRR6296792","SRR6296793","SRR6296796","SRR6296797","SRR6296798","SRR6296802","SRR6296803","SRR6296804","SRR6296808","SRR6296809","SRR6296810", "SRR6296815","SRR6296816","SRR6296817","SRR6296821","SRR6296822","SRR6296823","SRR6296828","SRR6296829","SRR6296832","SRR6296833")
morhua$Geneid <- str_remove(morhua$Geneid, "^gene:")
mo <- as.matrix(morhua[,-1])
mode(mo) <- "integer"
rownames(mo) <- morhua$Geneid
meta2<- data.frame(colnames(mo),Group = c("Control","0.01ug_BaP","1ug_BaP","Control","0.01ug_BaP","1ug_BaP","Control","0.01ug_BaP","1ug_BaP","Control","0.01ug_BaP","1ug_BaP","Control","0.01ug_BaP","1ug_BaP","Control","0.01ug_BaP","1ug_BaP","Control","1ug_BaP","Control","0.01ug_BaP"), ind = c("a","a","a","b","b","b","c","c","c","d","d","d","e","e","e","f","f","f","g","g","h","h"))
dds2 <- DESeqDataSetFromMatrix(countData=mo,
                               colData=meta2,
                               design=~Group +ind)
smallestGroupSizeMo <- 3
keep <- rowSums(counts(dds2) >= 10) >= smallestGroupSizeMo
dds2 <- dds2[keep,]

pcaData <- plotPCA(vsdmo, "Group", returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar"))
hull.datamo <- do.call(rbind, lapply(split(pcaData, list(pcaData$Group)), 
                                     function(x) x[chull(x[,c("PC1", "PC2")]),]))
ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_classic() 
?vst


morhuade1 <- DESeq(dds2) #object stores model
res0.01 <- results(morhuade1, contrast = c("Group", "0.01ug_BaP","Control"), alpha = 0.05, lfcThreshold = 1)
res1 <- results(morhuade1, contrast = c("Group","1ug_BaP","Control"), alpha = 0.05, lfcThreshold = 1)
summary(res0.01)

summary(res1)
res0.01df <- as.data.frame(res0.01)
res0.01df <- subset(res0.01df, log2FoldChange > 1 & padj < 0.05)
write.csv(res0.01df1, "res0.01morhuaup.csv")
res1df <- as.data.frame(res1)
res1df1 <- subset(res1df, log2FoldChange > 1 & padj < 0.05)
write.csv(res1df1, "res1morhuaup.csv")

plotMA(res0.01) #no sig
plotMA(res1) #little to no sig
vsdmo <- vst(dds2, blind = T)
plotPCA(vsdmo, "group")
?plotPCA
dfcod1 <- results(morhuade1, contrast = c("group", "0.01ug_BaP", "Control"), alpha=0.05) #comparing 0.01uM to 0
dfcod1 <- as.data.frame(dfcod1)
dfcod2 <- results(morhuade1, contrast = c("group", "1ug_BaP", "Control"), alpha=1e-5) #comparing 1uM to 0 
dfcod2 <- as.data.frame(dfcod2)

# volcano plot
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(res0.01,lab = "", x = 'log2FoldChange', y='pvalue')
EnhancedVolcano(res1,lab = "", x = 'log2FoldChange', y='pvalue')

# Morhua DE as numeric
meta2_num<- data.frame(colnames(mo), group = c(0,0.01,1,0,0.01,1,0,0.01,1,0,0.01,1,0,0.01,1,0,0.01,1,0,1,0,0.01), ind = c("a","a","a","b","b","b","c","c","c","d","d","d","e","e","e","f","f","f","g","g","h","h"))
dds2_num <- DESeqDataSetFromMatrix(countData=mo,
                                   colData=meta2_num,
                                   design=~ind+group)
morhuade_num <- DESeq(dds2_num)
morhuade_numdf <- results(morhuade_num)
morhuade_numdf <- as.data.frame(morhuade_numdf)
morhuade_num_sig <- subset(morhuade_numdf, padj < 0.05)
# write.csv(morhuade_num_sig, "morhuade_num_sig.csv")

library(DEXSeq)
library(tidyverse)
setwd("~/Desktop/Honors/featureCounts/Elegans/")
dexmorhua <- read.table("dexseq_counts_cod_v1", header=T)
dexmorhua <- dexmorhua[,-c(2:6)]
names(dexmorhua) <- c("TranscriptID","SRR6296791","SRR6296792","SRR6296793","SRR6296796","SRR6296797","SRR6296798","SRR6296802","SRR6296803","SRR6296804","SRR6296808","SRR6296809","SRR6296810", "SRR6296815","SRR6296816","SRR6296817","SRR6296821","SRR6296822","SRR6296823","SRR6296828","SRR6296829","SRR6296832","SRR6296833")
dexmorhua$TranscriptID <- str_remove(dexmorhua$TranscriptID, "^transcript:") 
dexmorhua$Transcript <- dexmorhua$TranscriptID
sample.data.morhua <- read.table("morhua_transcript_metadata.tsv")
names(sample.data.morhua) <- c("TranscriptID","GeneID")
sample.data.morhua$TranscriptID <- str_remove(sample.data.morhua$TranscriptID, "^transcript:")
sample.data.morhua$GeneID <- str_remove(sample.data.morhua$GeneID, "^gene:")  
sample.data.morhua$TranscriptID <- str_remove(sample.data.morhua$TranscriptID, "^*;")
sample.data.morhua$GeneID <- str_remove(sample.data.morhua$GeneID, "^*;")  
dexmorhua <- merge(dexmorhua, sample.data.morhua, by.all = TranscriptID)
Transcript <- dexmorhua$TranscriptID
Gene <- dexmorhua$GeneID
dexmo <- as.matrix(dexmorhua[,-c(1,24,25)])
mode(dexmo) <- "integer"
rownames(dexmo) <- dexmorhua$TranscriptID

dxdmo0.01 <- DEXSeqDataSet(countData = dexmo[,c(1,2,4,5,7,8,10,11,13,14,16,17,19,21,22)], #featureCounts output
                           sampleData = meta2[c(1,2,4,5,7,8,10,11,13,14,16,17,19,21,22),-3], #same as de
                           design=~sample + exon + group:exon, #but it's actually transcript
                           featureID = Transcript, #the transcripts
                           groupID = Gene) #genes they belong to
system.time({
  dxdmo0.01 <- estimateSizeFactors(dxdmo0.01)
  dxdmo0.01 <- estimateDispersions(dxdmo0.01, quiet = T)
  dxdmo0.01 <- testForDEU(dxdmo0.01, reducedModel =~sample + exon)
})

dxredmo0.01 <- DEXSeqResults(dxdmo0.01, independentFiltering = F)

qval_mor_dex0.01 <- perGeneQValue(dxredmo0.01)
dxr.g_mor_dex0.01 <- data.frame(gene=names(qval_mor_dex0.01),qval_mor_dex0.01)
subset(dxr.g_mor_dex0.01, qval_mor_dex0.01 < 0.05)


#### now 0 and 1
dxdmo1 <- DEXSeqDataSet(countData = dexmo[,c(1,3,4,6,7,9,10,12,13,15,16,18,19,20,21)], #featureCounts output
                        sampleData = meta2[c(1,3,4,6,7,9,10,12,13,15,16,18,19,20,21),], #same as de
                        design=~sample + exon + group:exon, #but it's actually transcript
                        featureID = Transcript, #the transcripts
                        groupID = Gene) #genes they belong to

system.time({
  dxdmo1 <- estimateSizeFactors(dxdmo1)
  dxdmo1 <- estimateDispersions(dxdmo1, quiet = T)
  dxdmo1 <- testForDEU(dxdmo1, reducedModel =~sample + exon)
})

dxredmo1 <- DEXSeqResults(dxdmo1, independentFiltering = F)

qval_mor_dex1 <- perGeneQValue(dxredmo1)
dxr.g_mor_dex1 <- data.frame(gene=names(qval_mor_dex1),qval_mor_dex1)
subset(dxr.g_mor_dex1, qval_mor_dex1 < 0.05)
##### Morhua DEXSeq as Num

meta2_num
meta2_num_dex <- meta2_num[,c(1,2)]
dex_num_mo <- DEXSeqDataSet(countData=dexmo,
                            sampleData=meta2_num_dex,
                            design=~sample+exon+group:exon,
                            featureID = Transcript,
                            groupID = Gene)

system.time({
  dex_num_mo <- estimateSizeFactors(dex_num_mo)
  dex_num_mo <- estimateDispersions(dex_num_mo, quiet = T)
  dex_num_mo <- testForDEU(dex_num_mo, reducedModel =~sample + exon) 
})
dex_num_mo_res <- DEXSeqResults(dex_num_mo, independentFiltering = F)
qval_mo_dex <- perGeneQValue(dex_num_mo_res)
dxr.g_mo_dex <- data.frame(gene=names(qval_mo_dex),qval_mo_dex)
subset(dxr.g_mo_dex, qval_mo_dex < 0.05)

library(DEXSeq)
library(tidyverse)
setwd("~/Desktop/Honors/featureCounts/Elegans/")

exonmorhua <- read.table("exoncountsmorhua.txt", header = T)
exonmorhua <- exonmorhua[,-c(2:6)]
names(exonmorhua) <- c("Geneid","SRR6296791","SRR6296792","SRR6296793","SRR6296796","SRR6296797","SRR6296798","SRR6296802","SRR6296803","SRR6296804","SRR6296808","SRR6296809","SRR6296810", "SRR6296815","SRR6296816","SRR6296817","SRR6296821","SRR6296822","SRR6296823","SRR6296828","SRR6296829","SRR6296832","SRR6296833")
exonmo <- as.matrix(exonmorhua[,-1])
mode(exonmo) <- "integer"
rownames(exonmo) <- exonmorhua$Geneid

exonmo <- as.data.frame(exonmo)
exonmo$GeneID <- gsub("\\..*$", "", rownames(exonmo))
exonmo$exonname <- rownames(exonmo)
moexonname <- exonmo$exonname
moexongene <- exonmo$GeneID
exonmo <- exonmo[,-c(23,24)]

dxexonmo0.01 <- DEXSeqDataSet(countData = exonmo[,c(1,2,4,5,7,8,10,11,13,14,16,17,19,21,22)], #featureCounts output
                              sampleData = meta2[c(1,2,4,5,7,8,10,11,13,14,16,17,19,21,22),], #same as de
                              design=~sample + exon + group:exon,
                              featureID = moexonname, #the exons
                              groupID = moexongene) #genes they belong to
smallestGroupSizeMor <- 3
keep <- rowSums(counts(dxexonmo0.01) >= 10) >= smallestGroupSizeMor
dxexonmo0.01 <- dxexonmo0.01[keep,]

system.time({
  dxexonmo0.01 <- estimateSizeFactors(dxexonmo0.01)
  dxexonmo0.01 <- estimateDispersions(dxexonmo0.01, quiet = T)
  dxexonmo0.01 <- testForDEU(dxexonmo0.01, reducedModel =~sample + exon)
})

dxexonredmo0.01 <- DEXSeqResults(dxexonmo0.01, independentFiltering = F)
qval_mor_exon_dex0.01 <- perGeneQValue(dxexonredmo0.01)
dxr.g_mor_exon_dex0.01 <- data.frame(gene=names(qval_mor_exon_dex0.01),qval_mor_exon_dex0.01)
subset(dxr.g_mor_exon_dex0.01, qval_mor_exon_dex0.01 < 0.05)


###### 1 uM comparative #####
dxexonmo1 <- DEXSeqDataSet(countData = exonmo[,c(1,3,4,6,7,9,10,12,13,15,16,18,19,20,21)], #featureCounts output
                           sampleData = meta2[c(1,3,4,6,7,9,10,12,13,15,16,18,19,20,21),], #same as de
                           design=~sample + exon + group:exon,
                           featureID = moexonname, #the exons
                           groupID = moexongene) #genes they belong to
smallestGroupSizeMor <- 3
keep <- rowSums(counts(dxexonmo1) >= 10) >= smallestGroupSizeMor
dxexonmo1 <- dxexonmo1[keep,]

system.time({
  dxexonmo1 <- estimateSizeFactors(dxexonmo1)
  dxexonmo1 <- estimateDispersions(dxexonmo1, quiet = T)
  dxexonmo1 <- testForDEU(dxexonmo1, reducedModel =~sample + exon)
})

dxexonredmo1 <- DEXSeqResults(dxexonmo1, independentFiltering = F)
qval_mor_exon_dex1 <- perGeneQValue(dxexonredmo1)
dxr.g_mor_exon_dex1 <- data.frame(gene=names(qval_mor_exon_dex1),qval_mor_exon_dex1)
subset(dxr.g_mor_exon_dex1, qval_mor_exon_dex1 < 0.05)


exonmotest <- merge(exonmo, Transcript, by.x = rownames(exonmo), by.y = Transcript$x)
summary(Transcript)

###### exon as num #########
#meta2_num_dex
dexex_num_mo <- DEXSeqDataSet(countData=exonmo,
                              sampleData=meta2_num_dex,
                              design=~sample+exon+group:exon,
                              featureID = moexonname,
                              groupID = moexongene)

system.time({
  dexex_num_mo <- estimateSizeFactors(dexex_num_mo)
  dexex_num_mo <- estimateDispersions(dexex_num_mo, quiet = T)
  dexex_num_mo <- testForDEU(dexex_num_mo, reducedModel =~sample + exon) 
})
dexex_num_mo_res <- DEXSeqResults(dexex_num_mo, independentFiltering = F)
qvalex_mo_dex <- perGeneQValue(dexex_num_mo_res)
dxr.gex_mo_dex <- data.frame(gene=names(qvalex_mo_dex),qvalex_mo_dex)
abcd <- subset(dxr.gex_mo_dex, qvalex_mo_dex < 0.05)
write.csv(abcd, "EXMorhua_NUM.csv")

library(UpSetR)
library(tidyverse)



sample.data.elegans <- read.table("transcript2gene.txt")
names(sample.data.elegans) <- c("WBGene","TranscriptID")
sample.data.elegans$TranscriptID <- str_remove(sample.data.elegans$TranscriptID, "^Name=")
sample.data.elegans$WBGene <- str_remove(sample.data.elegans$WBGene, "^Parent=Gene:")
res5df1 <- read.csv("res5elegans.csv")
res20df3 <- read.csv("Allsiggenesres20.csv") 
resnumalphadf<- read.csv("resnumalphadf.csv")
dxr.gnum_sig <- read.csv("elegansDEXseqNUM.csv")
dxrgex05sig <- read.csv("elegansexon_05.csv")
dxrgex020sig <- read.csv("elegansexon_020.csv")
dxrgexnumsig <- read.csv("elegansexon_NUM.csv")


ups_elegans <- data.frame(gene=unique(sample.data.elegans$WBGene), 
                          `Gene_0vs5`=0, `Gene_0vs20`=0, `Gene_Linear`=0,`Transcript_Linear`=0,
                          `Exon_0vs5`=0,`Exon_0vs20`=0,`Exon_Linear`=0)
ups_elegans[ups_elegans$gene %in% res5df1$X,2] <- 1
ups_elegans[ups_elegans$gene %in% res20df3$X,3] <- 1
ups_elegans[ups_elegans$gene %in% resnumalphadf$X,4] <- 1
ups_elegans[ups_elegans$gene %in% dxr.gnum_sig$X,5] <- 1
ups_elegans[ups_elegans$gene %in% dxrgex05sig$X,6] <- 1
ups_elegans[ups_elegans$gene %in% dxrgex020sig$X,7] <- 1
ups_elegans[ups_elegans$gene %in% dxrgexnumsig$X,8] <- 1

upset(ups_elegans, nset=7, order.by = "freq", point.size = 3, text.scale = 1.3, 
      sets.bar.color = c("gold2","maroon1","gold2","maroon1",
      "gold2","maroon1","purple"),sets = c("Gene_0vs5","Gene_0vs20", 
      "Gene_Linear","Transcript_Linear","Exon_0vs5","Exon_0vs20","Exon_Linear"),
      queries = list(list(query = intersects, params = list("Exon_Linear",
      "Gene_Linear"),color = "skyblue", active = T), list(query = intersects, 
      params = list("Transcript_Linear","Exon_Linear"),color = "blue", active = T),
      list(query = intersects,params = list("Exon_0vs5","Gene_Linear"),color = "skyblue", active = T),
      list(query = intersects,params = list("Exon_0vs5","Gene_0vs5","Exon_0vs20","Gene_0vs20","Gene_Linear"),
      color ="skyblue", active = T), list(query = intersects, params = list("Gene_0vs20",
      "Gene_Linear"),color = "gold2", active = T), list(query = intersects, 
      params = list("Gene_0vs20","Gene_0vs5","Gene_Linear"),color = "gold2", active = T),
      list(query = intersects, params = list("Exon_0vs20","Exon_0vs5","Exon_Linear"),color = 
      "maroon1", active = T), list(query = intersects, params = list("Exon_0vs20",
      "Exon_Linear"),color = "maroon1", active = T), list(query = intersects, params = list("Exon_0vs5",
      "Exon_0vs20"), color = "maroon1", active = T), list(query = intersects, params = list("Gene_0vs5","Gene_0vs20"),
      color = "gold2", active = T)), set_size.show = T)


sample.data.morhua <- read.table("morhua_transcript_metadata.tsv")
names(sample.data.morhua) <- c("TranscriptID","GeneID")
sample.data.morhua$TranscriptID <- str_remove(sample.data.morhua$TranscriptID, "^transcript:")
sample.data.morhua$GeneID <- str_remove(sample.data.morhua$GeneID, "^gene:")  
sample.data.morhua$TranscriptID <- str_remove(sample.data.morhua$TranscriptID, "^*;")
sample.data.morhua$GeneID <- str_remove(sample.data.morhua$GeneID, "^*;")  
res1morhua <- read.csv("res1morhuaup.csv")
resdelin <- read.csv("morhuade_num_sig.csv")
exonlin <- read.csv("EXMorhua_NUM.csv")
view(resdelin)

ups_morhua <- data.frame(gene=unique(sample.data.morhua$GeneID), 
                         `Gene_Linear`=0,`Gene_0vs1`=0, `Exon_Linear`=0)
ups_morhua[ups_morhua$gene %in% resdelin$X,2] <- 1
ups_morhua[ups_morhua$gene %in% res1morhua$X,3] <- 1
ups_morhua[ups_morhua$gene %in% exonlin$X,4] <- 1

upset(ups_morhua, sets.bar.color=c("lightsalmon","firebrick1","lightsalmon"), order.by = "freq", 
      point.size = 3, text.scale = 1.3)

grid.arrange(upset_elegans, upset_morhua)

EnhancedVolcano(res5,lab = "", x = 'log2FoldChange', y='padj',
                gridlines.major = F, gridlines.minor = F, pCutoff = 0.05, labSize = 3)
EnhancedVolcano(res20,lab = "", x = 'log2FoldChange', y='padj', 
                gridlines.major = F, gridlines.minor = F, pCutoff = 0.05, ylim = c(0,100), xlim = c(-5,7),
                labSize = 5)
EnhancedVolcano(res1,lab = "", x = 'log2FoldChange', y='padj', ylim = c(0,25), xlim = c(-2.5,7.5),
                gridlines.major = F, gridlines.minor = F, pCutoff = 0.05, labSize = 3, legendPosition = "bottom")
