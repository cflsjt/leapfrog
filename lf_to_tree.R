#!/usr/bin/env Rscript
#library(ape)
#suppressPackageStartupMessages()
suppressMessages(library(ggplot2))
suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)

dd <- read.table(args[1])

nocol <- (length(colnames(dd)) - 4) / 2

bindata <- dd[colnames(dd)[5:(5+nocol-1)]]
scoredata <- dd[colnames(dd)[(5+nocol):(5+nocol+nocol-1)]]

data.bindist <- dist(t(bindata))
data.scoredist <- dist(t(scoredata))

dd$binsum <- rowSums(bindata)

chrnames <- c("chr1", "chr2", "chr3", "chr4",
              "chr5", "chr6", "chr7", "chr8",
              "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16",
              "chr17", "chr18", "chr19", "chr1_random",
              "chr3_random", "chr4_random",
              "chr5_random", "chr7_random",
              "chr9_random", "chr10_random",
              "chr11_random", "chr12_random",
              "chr13_random", "chr16_random",
              "chr17_random", "chr18_random",
              "chrUn")

#order chromosomes
## set the levels in order we want
dd <- within(dd, 
             chr <- factor(chr, levels=chrnames))


#sample stats
colSums(bindata)
png(paste(args[2], '.samples.png', sep=''), width=600, height=500)
qplot(y=colSums(bindata), x=colnames(bindata), geom="bar")
dummy<-dev.off()

png(paste(args[2], '.bin.tree.png', sep=''), width=600, height=500)
plot(nj(data.bindist), type='unro')
dummy<-dev.off()

png(paste(args[2], '.score.tree.png', sep=''), width=600, height=500)
plot(nj(data.scoredist), type='unro')
dummy<-dev.off()

#dive into families
png(paste(args[2], '.family.histo.png', sep=''), width=600, height=500)
qplot(dd$family) + opts(axis.text.x=theme_text(angle=-90, hjust=0))
dummy<-dev.off()

famtable <- table(factor(dd$family))
superfamilies <- rownames(famtable[famtable > 50])
superfamilies
ddd <- dd[dd$family %in% superfamilies,]

png(paste(args[2], '.samplefreq.png', sep=''), width=600, height=500)
ggplot(ddd, aes(binsum)) + geom_bar()  + scale_x_discrete(name="Present in X samples")
dummy<-dev.off()

png(paste(args[2], '.family.per.samplefreq.png', sep=''), width=600, height=500)
ggplot(ddd, aes(binsum, fill=family)) +
  geom_bar(position='fill') +
  scale_x_discrete(name="Present in X samples") +
  scale_y_continuous(name="Fraction")
dummy<-dev.off()

png(paste(args[2], '.family.per.chr.png', sep=''), width=600, height=500)
ggplot(ddd, aes(chr, fill=family)) +
  geom_bar(position='fill') +
  scale_x_discrete(name="Chromosome") +
  scale_y_continuous(name="Fraction") + opts(axis.text.x=theme_text(angle=-90, hjust=0))
dev.off()


#dd$binsum <- rowSums(bindata)


#no of REFS present per sample
#bin.colsums = colSums(bindata)
#bin.rowsums = colSums(bindata)
#png(paste(args[2], 'sums.per.sample.png', sep=''), width=600, height=500 )
#g <- qplot(data.colsums, geom="bar") #, stat="identity", geom='bar')
#dev.off()

## data.rowsums = rowSums(data)

