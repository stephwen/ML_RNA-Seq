library(TCGAbiolinks)
library(dplyr)
library(DT)
library("DESeq2")
library(stringr)
library(randomForest)
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
  FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")

  query <- GDCquery(project = CancerProject,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")

  samplesDown <- getResults(query,cols=c("cases"))

  # tumor samples
  dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "TP")

  # healthy samples
  dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown, typesample = "NT")

  queryDown <- GDCquery(project = CancerProject, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP, dataSmNT))

  # the download command may fail, for various reasons (connection timeout, etc)
  # so we-re-run it (it detects files already downloaded) until it outputs
  # a specific message on stderr (the function doesn't have output values)

  downloadResponse<-c("", "", "")
  while (downloadResponse[3] != "All samples have been already downloaded") {
    downloadResponse<-capture.output(GDCdownload(query = queryDown, directory = DataDirectory, files.per.chunk = 10), type = "message")
  }

  # we re-run it again, just in case
  GDCdownload(query = queryDown, directory = DataDirectory, files.per.chunk = 10)

  dataPrep <- GDCprepare(query = queryDown, 
                         save = TRUE, 
                         directory =  DataDirectory,
                         save.filename = FileNameData)

  dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep, datatype = "HTSeq - Counts")  

  samplesList<-list()
  cond<-NULL

  for (i in (1:dim(dataPrep)[2])) {
    sampleName<-substr(colnames(dataPrep)[i], 0, 15)
    if (is.null(samplesList[[sampleName]])) { # we do this because some samples might be present multiple times
      vector<-as.matrix(dataPrep[,i])
      write.table(vector, sep="\t", file=paste0(DataDirectory, "/",i,"_",sampleName,"_","HTSeq_Counts",".csv"), col.names = FALSE, quote=F)
      # we output the HTSeq-count data to individual files which will be easily loaded at the DESeq2 step
      samplesList[[sampleName]]=i
      if (colnames(dataPrep)[i] %in% dataSmNT) {
        cond<-c(cond, "normal")
      } else { cond<-c(cond, "cancer") }
    }
  }

  write.csv(cond, file=paste0(DataDirectory, "/cond.csv"))
  save(cond, file = paste0(DataDirectory, "/cond.robject"))

}
