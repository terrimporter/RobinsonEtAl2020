# Teresita M. Porter, Jan. 2, 2020

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(dplyr) # group_by
library(ggrepel) #geom_text_repel

# Read in cat.csv
A <- read.csv(file="cat.csv", head=TRUE)

# fix sites numbers Pres2: S1, S2, S3; Pres3: S4, S5, S6 in the SampleName field
A$SampleName <- gsub("(Pres3_\\w+_)S1(_\\w*_\\w*_S\\d+)", "\\1S4\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S2(_\\w*_\\w*_S\\d+)", "\\1S5\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S3(_\\w*_\\w*_S\\d+)", "\\1S6\\2", A$SampleName)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")))
names(A.1)[32:37] <- c("Experiment","Treatment","Site", "Replicate","PrePost", "IlluminaSample")

# Focus on Arthropoda
A.2 <- A.1[A.1$Phylum=="Arthropoda",]

# Focus on Post-evaporation only
A.3 <- A.2[A.2$PrePost=="POST",]

# Only plot good genera gBP >= 0.30
A.4 <- A.3[A.3$gBP >= 0.30,]

# Create custom field for cast
A.4$OrderGenusGlobalESV <- paste(A.4$Order, A.4$Genus, A.4$GlobalESV, sep=";")

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
A.5.esv <- reshape2::dcast(A.4, SampleName ~ OrderGenusGlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.5.esv) <- A.5.esv$SampleName
A.5.esv$SampleName <- NULL

# remove columns with only zeros
esv.notnull<-A.5.esv[,colSums(A.5.esv) !=0]

# remove rows with only zeros & edit rownames
esv.notnull2 <- esv.notnull[rowSums(esv.notnull) !=0,]

# calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2), prob=0.15)
# 15% 
# 8932 

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# convert to presence absence
rare.mat[rare.mat>0] <-1

# convert to df
df <- as.data.frame(rare.mat)

# grab just ANTI and ETOH
anti <- df[grepl("ANTI", rownames(df)),]
etoh <- df[grepl("ETOH", rownames(df)),]

# Sum ESVs accross samples
anti_sums <- colSums(anti)
etoh_sums <- colSums(etoh)

# Convert to presence absence (i.e. pool accross samples above, then mark when detected with a 1)
anti_sums[anti_sums>0] <-1
etoh_sums[etoh_sums>0] <-1

# combine into df
anti_etoh <- data.frame(cbind(anti_sums, etoh_sums))

# move rownames to first column
setDT(anti_etoh, keep.rownames = TRUE)[]
names(anti_etoh)[1] <- "Taxon_GlobalESV"

# remove last to fields of Taxon
t <- data.frame(anti_etoh, do.call(rbind, str_split(anti_etoh$Taxon_GlobalESV,";")))
names(t)[4:6] <- c("Order","Genus","GlobalESV")

# remove Taxon_GlobalESV
t$Taxon_GlobalESV <- NULL

# Create Taxon
t$Taxon <- paste(t$Order, t$Genus, sep=";")

# group by taxon, then sum
t2 <- t %>%
  group_by(Taxon) %>%
  summarize(sumAnti=sum(anti_sums),
            sumEtoh=sum(etoh_sums))

# split order and genus into their own columns
t3 <- data.frame(t2, do.call(rbind, str_split(t2$Taxon,";")))
names(t3)[4:5] <- c("Order","Genus")

# Compare richness by site
p <- ggplot(t3, aes(x=sumAnti, y=sumEtoh)) +
  geom_point(aes(color=factor(Order)), position=position_jitter()) +
  labs(x="Antifreeze ESVs (log10)", y="Ethanol ESVs (log10)", colour = "Order") +
  geom_abline(intercept = 0, slope = 1, linetype = 3) +
  scale_x_continuous(trans='log10', expand=expand_scale(mult=c(0.05,0.05))) +
  scale_y_continuous(trans='log10', expand=expand_scale(mult=c(0.05,0.05))) +
  geom_text_repel(aes(x = sumAnti, 
                      y = sumEtoh, 
                      label = Genus),
                  data=t3[t3$sumAnti>=2 & t3$sumEtoh>=2,],
                  size=2) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=10),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.position="bottom")

ggsave("Order_xy.pdf", p, width = 8)
# based on normalized/rarefied data
# Arthropoda, Post evaporation only
# Confidently identified genera only
# controls excluded
 