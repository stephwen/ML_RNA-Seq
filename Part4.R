library(TCGAbiolinks)
library(dplyr)
library(DT)
library("DESeq2")
library(stringr)
library(ROCR)
library(caret)
library(survival)

baseDir<-"/Users/steph/github/ML_RNA-Seq"

setwd(baseDir)

# change accordingly if you're only interested in one of he datasets
datasets<-c("BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

for (dataset in datasets) {

	CancerProject <- paste0("TCGA-", dataset)
	DataDirectory <- paste0(baseDir, "/", gsub("-","_",CancerProject))
	load(paste0(DataDirectory, "/cond.robject"))
	load(paste0(DataDirectory, "/RF-geneRanking.robject"))
	load(paste0(DataDirectory, "/DESeq2ResultsSig.robject"))
	load(paste0(DataDirectory, "/VS.robject"))

	topDESeq2Genes<-rownames(DESeq2ResultsSig)[1:20]
	topRFGenes<-as.character(RF_ResultsOrdered[1:20])



	# survival analysis on the VS

	clinical_patient_Cancer <- GDCquery_clinic(CancerProject,"clinical")
	survival<-clinical_patient_Cancer[c("submitter_id", "days_to_death", "days_to_last_follow_up")]

	surv2<-survival[survival[,1] %in% substr(colnames(VS), 0, 12),]
	colnames(VS)<-substr(colnames(VS), 0, 12)
	VS<-VS[,!duplicated(colnames(VS))]
	VS<-VS[,colnames(VS) %in% surv2[,1]]

	na.zero <- function (x) {
	    x[is.na(x)] <- 0
	    return(x)
	}

	survObj<-Surv(pmax(na.zero(surv2$days_to_death), na.zero(surv2$days_to_last_follow_up)), !(is.na(surv2$days_to_death)))

	VSdf<-as.data.frame(t(VS))

	# we test the association with survival between the 20 gene signatures
	# based on the top 20 genes from DESeq2
	DESeq2pvals<-NULL
	for (i in c(1:length(topDESeq2Genes))) {
		formula = as.formula(paste("survObj ~", paste(as.list(topDESeq2Genes[1:i]), collapse = " + ")))
		res.cox1 <- coxph(formula, data = VSdf)
		binaryVec <- res.cox1$linear.predictors < median(res.cox1$linear.predictors)
		if (sum(binaryVec)==0) { logRankPval<-"NA" } else {
			logrank<-survdiff(survObj ~ binaryVec)
			logRankPval<-1 - pchisq(logrank$chisq, 1)
			print(paste(i, round(logRankPval, 3)))
		}	
		DESeq2pvals<-rbind(DESeq2pvals, c(topDESeq2Genes[i], logRankPval))
	}
	# we test the association with survival between the 20 gene signatures
	# based on the top 20 genes ranked by the random forests permutation importance
	RFpvals<-NULL
	for (i in c(1:length(topRFGenes))) {
		formula = as.formula(paste("survObj ~", paste(as.list(topRFGenes[1:i]), collapse = " + ")))
		res.cox1 <- coxph(formula, data = VSdf)
		binaryVec <- res.cox1$linear.predictors < median(res.cox1$linear.predictors)
		if (sum(binaryVec)==0) { logRankPval<-"NA" } else {
			logrank<-survdiff(survObj ~ binaryVec)
			logRankPval<-1 - pchisq(logrank$chisq, 1)
			print(paste(i, round(logRankPval, 3)))
		}	
		RFpvals<-rbind(RFpvals, c(topRFGenes[i], logRankPval))
	}

	# output the results
	outputTable<-cbind(DESeq2pvals, RFpvals)
	colnames(outputTable)<-c("DESeq2 (gene)", "DESeq2 (multi-gene p-value)", "RF (gene)", "RF (multi-gene p-value)")
	write.table(outputTable, file=paste0(DataDirectory, "/outputTable.csv"), quote=F, sep=",", row.names=FALSE)
	save(outputTable, file = paste0(DataDirectory, "/outputTable.robject"))

}