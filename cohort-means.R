#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

segSites = read.table(args[1])
indStats = read.table(args[2],header=T)
cohortName = args[3]

cohortMeans = data.frame("Cohort"=cohortName)
indMeans = as.data.frame(t(colMeans(indStats[,2:ncol(indStats)])))
cohortMeans = cbind(cohortMeans,indMeans)
cohortMeans$SegregatingSites = as.numeric(segSites$V1)


if (file.exists("outputs/meanCohortStats.txt")) {
  write.table(cohortMeans,"outputs/meanCohortStats.txt",quote=F,row.names=F,col.names=F,append=T)
} else {
  write.table(cohortMeans,"outputs/meanCohortStats.txt",quote=F,row.names=F)
}
