#!/usr/bin/env Rscript

# Arguments are positional and should include:
# plink.eigenvec # The eigenvectors
# plink.imiss # The individual missingness
# optionally, an output prefix, to include path or file prefix
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)<3){
  prefix = "output"
} else {
  prefix = args[3]
}

vecs = read.table(args[1])
miss = read.table(args[2],header=T)

# Plink separates "family" from id, so paste these back together
vecs$ids = paste(vecs$V1,vecs$V2,sep="_")

# add missingness to the vecs data frame to color points by missingness
vecs$miss = miss$F_MISS

pdf(paste0(prefix,"_PCA.pdf"),width=5,height=4)
p=ggplot(vecs,aes(x=V3,y=V4,color=miss))+
  geom_point()+
  xlab("PC1")+
  ylab("PC2")+
  ggtitle(prefix)+
  theme_minimal()
print(p)
dev.off()

# A second plot that adds the individual label for assessment
pdf(paste0(prefix,"_PCA_labels.pdf"),width=5,height=4)
p=ggplot(vecs,aes(x=V3,y=V4,color=miss))+
  geom_point()+
  geom_label(aes(label=ids),nudge_y=0.1)+
  xlab("PC1")+
  ylab("PC2")+
  ggtitle(prefix)+
  theme_minimal()+
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_y_continuous(expand = expansion(mult = 0.25))
print(p)
dev.off()
