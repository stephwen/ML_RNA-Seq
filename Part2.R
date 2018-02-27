library(TCGAbiolinks)
library(dplyr)
library(DT)
library("DESeq2")
library(stringr)
library(ROCR)
library(caret)

baseDir<-"/Users/steph/github/ML_RNA-Seq"

setwd(baseDir)

# change accordingly if you're only interested in one of he datasets
datasets<-c("BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

for (dataset in datasets) {

	CancerProject <- paste0("TCGA-", dataset)
	DataDirectory <- paste0(baseDir, "/", gsub("-","_",CancerProject))
	load(paste0(DataDirectory, "/cond.robject"))
  	
  	sampleFiles <- grep("HTSeq_Counts",list.files(DataDirectory),value=TRUE)
	samples <- str_split_fixed(sampleFiles, "_", 3)[,2]

	# DESeq2 normalization
	df <- data.frame(samples, sampleFiles, cond)
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = df, directory = DataDirectory, design= ~ cond)
	dds2 = estimateSizeFactors(ddsHTSeq)
	normalizedCounts<-counts(dds2, normalized=TRUE)
	write.csv(normalizedCounts, file=paste0(DataDirectory, "/DESeq2-normalizedCounts.csv"))

	# split normalizedCounts into an LS and a VS
	# LS = 50% of all cancers + all healthy
	# VS = 50% of all cancers
	# gene selection will be performed on LS
	# survival analysis will be performed on VS

	LS<-normalizedCounts[,(sum(cond=="cancer")-sum(cond=="normal")):dim(normalizedCounts)[2]]	
	# start at end of all cancers minus nb of normal samples, up until the last (normal) sample

	VS<-normalizedCounts[,1:(sum(cond=="cancer")-sum(cond=="normal"))-1]
	# going from the 1st sample up until nb all cancer samples minus nb all healthy samples minus 1 (inverse of LS)

	# LS but in DESeq2 format
	dds <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, (sum(cond=="cancer")-sum(cond=="normal")):dim(normalizedCounts)[2]]

	# do the DESeq2 magic
	ddsHTSeq2 <- DESeq(dds)

	res <- results(ddsHTSeq2, alpha=0.05)
	resOrdered <- res[order(res$padj),]
	resSig <- subset(resOrdered, padj < 0.05)
	write.csv(as.data.frame(resOrdered), file=paste0(DataDirectory, "/DESeq2-results.csv"))

	# save the learning set and the validation set for the following steps
	save(LS, file = paste0(DataDirectory, "/LS.robject"))
	save(VS, file = paste0(DataDirectory, "/VS.robject"))

}