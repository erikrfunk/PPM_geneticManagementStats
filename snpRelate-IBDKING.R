#!/usr/bin/env Rscript

library(SNPRelate)

args <- commandArgs(trailingOnly = TRUE)
if(length(args)<2){
  prefix = "output"
} else {
  prefix = args[2]
}


snpgdsVCF2GDS(args[1], "gdsin.gds", method="biallelic.only")
fin = snpgdsOpen("gdsin.gds")

ibd_robust = snpgdsIBDKING(fin, type="KING-robust", ,autosome.only=FALSE)
ibd_robust_mat <- snpgdsIBDSelection(ibd_robust)

ibd_homo = snpgdsIBDKING(fin, type="KING-homo", ,autosome.only=FALSE)
ibd_homo_mat = snpgdsIBDSelection(ibd_homo)

write.table(ibd_robust_mat,paste0(prefix,"_kingRobust.txt"),row.names=F,quote=F)
write.table(ibd_homo_mat,paste0(prefix,"_kingHomogenous.txt"),row.names=F,quote=F)

snpgdsClose(fin)
