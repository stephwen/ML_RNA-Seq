library(TCGAbiolinks)
library(dplyr)
library(DT)
library("DESeq2")
library(stringr)
library(randomForest)
library(ROCR)
library(caret)
library(survival)
library(ranger)

baseDir<-"/Users/steph/github/ML_RNA-Seq"

setwd(baseDir)

# change accordingly if you're only interested in one of he datasets
datasets<-c("BRCA", "COAD", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

for (dataset in datasets) {

	CancerProject <- paste0("TCGA-", dataset)
	DataDirectory <- paste0(baseDir, "/", gsub("-","_",CancerProject))
	load(paste0(DataDirectory, "/cond.robject"))
  	load(paste0(DataDirectory, "/LS.robject"))
  	load(paste0(DataDirectory, "/VS.robject"))

	RF_LS<-t(LS)	# transpose the LS matrix for use with ranger

	cond2<-cond[(sum(cond=="cancer")-sum(cond=="normal")):(dim(LS)[2]+dim(VS)[2])]		# conditional vector [ indexes of the LS ]

	# use the appropriate data format for ranger
	rangerCond<-as.factor(cond2)
	RF_LS<-cbind(RF_LS, rangerCond)

	# here we only extract the permutation importance, but it's very straightforward to extract the impurity importance too, see commented lines below
	rf.ranger <- ranger(dependent.variable.name="rangerCond", data = RF_LS, num.trees = 10001, mtry = 236, write.forest = FALSE, probability = FALSE, importance="permutation")
	permutationRanking<-rf.ranger$variable.importance[order(rf.ranger$variable.importance, decreasing=T)]

	# rf.ranger <- ranger(dependent.variable.name="rangerCond", data = RF_LS, num.trees = 10001, mtry = 236, write.forest = FALSE, probability = FALSE, importance="impurity")
	# impurityRanking<-rf.ranger$variable.importance[order(rf.ranger$variable.importance, decreasing=T)]

	RF_ResultsOrdered<-names(permutationRanking)
	write.table(RF_ResultsOrdered, file=paste0(DataDirectory, "/RF-geneRanking.csv"), row.names = FALSE, col.names = FALSE, quote=F)
	save(RF_ResultsOrdered, file = paste0(DataDirectory, "/RF-geneRanking.robject"))

}
  	