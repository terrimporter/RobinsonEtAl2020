# Teresita M. Porter, Dec. 23, 2019

library(reshape2) # dcast
library(ggplot2) # ggplot

###################################################################
# Plot phyla vs ESVs/reads
# Calculate number of Arthropoda ESVs/reads for each amplicon
###################################################################

# Read infile
A <- read.csv(file="cat.csv", head=TRUE)

# fix sites numbers Pres2: S1, S2, S3; Pres3: S4, S5, S6 in the SampleName field
A$SampleName <- gsub("(Pres3_\\w+_)S1(_\\w*_\\w*_S\\d+)", "\\1S4\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S2(_\\w*_\\w*_S\\d+)", "\\1S5\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S3(_\\w*_\\w*_S\\d+)", "\\1S6\\2", A$SampleName)

# Summarize ESVs in all unique phyla (pooled across samples)
A.esv <- dcast(A, Phylum ~ . , value.var="GlobalESV", function (x) {length(unique(x))})
names(A.esv)<-c("Phyla","GlobalESV")

# Sort by descending ESVs
A.esv.desc<-A.esv[order(-A.esv$GlobalESV),]

# Summarize reads in all detected phyla
A.read <- dcast(A, Phylum ~ . , value.var="ESVsize", fun.aggregate=sum)
names(A.read)<-c("Phyla","ESVsize")

# Sort by descending reads
A.read.desc<-A.read[order(-A.read$ESVsize),]

# calc proportions
A.esv.phyla<-A.esv.desc[,1]
A.esvprop<-round(A.esv.desc[,2]/sum(A.esv.desc[,2])*100,digits=2)
A.read.phyla<-A.read.desc[,1]
A.readprop<-round(A.read.desc[,2]/sum(A.read.desc[,2])*100,digits=2)

# Create esv prop table
A.esv.table<-data.frame(phylum=A.esv.phyla, esv=A.esvprop)

# Create read prop table
A.read.table<-data.frame(phylum=A.read.phyla, read=A.readprop)

# Keep top 10
A.esv.top<-A.esv.table[1:10,]
A.read.top<-A.read.table[1:10,]

# Create 'other' df
A.esv.other<-A.esv.table[-(1:10),]
A.read.other<-A.read.table[-(1:10),]

# Sum 'other' line
A.esv.other.sum <- sum(A.esv.other$esv)
A.read.other.sum <- sum(A.read.other$read)

# Create df record for other
A.esv.other.sum.df <- data.frame("phylum"="Other", "esv"=A.esv.other.sum)
A.read.other.sum.df <- data.frame("phylum"="Other", "read"=A.read.other.sum)

# Add other to top df
A.esv.rbind <- rbind(A.esv.top, A.esv.other.sum.df)
A.read.rbind <- rbind(A.read.top, A.read.other.sum.df)

# outer join esvs and reads
A.table<-merge(A.esv.rbind, A.read.rbind, by="phylum", all=TRUE)

# Create long form form ggplot
A.long<-melt(A.table, id=c("phylum"))

# create factors
A.long$variable <- factor(A.long$variable,
                          levels=c("esv","read"),
                          labels=c("ESVs","Reads"))

# proprtions for top 10
p<-ggplot(data=A.long, aes(x=variable, y=value, fill=phylum, label=value)) +
  geom_bar(stat="identity") +
  geom_text(size = 2, position = position_stack(vjust = 0.5)) +
  labs(y="Proportions", x="Top 10 phyla") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=10),
        axis.title.x=element_blank(),
        plot.title=element_text(size=10),
        legend.position="bottom",
        legend.text=element_text(size=6),
        legend.title=element_blank(),
        legend.key.size = unit(0.45, "cm")) +
  guides(fill = guide_legend(nrow = 3))

ggsave("FigS3_phyla.pdf", p, width=8)

##############################################
# Supplementary Table Info
##############################################

# Get Arthropoda & BR5
BR5.Arth <- A[A$Phylum=="Arthropoda" &
            grepl("BR5", A$GlobalESV),]

# Sum Arthropoda & BR5 ESV counts
length(unique(BR5.Arth$GlobalESV))
# 1,711

# Sum Arthropoda & BR5 read counts
sum(BR5.Arth$ESVsize)
# 683,908

# Get Arthropoda & F230R
F230R.Arth<-A[A$Phylum=="Arthropoda" &
              grepl("F230R", A$GlobalESV),]

# Sum Arthropoda & F230R ESV counts
length(unique(F230R.Arth$GlobalESV))
# 1320

# Sum Arthropoda & F230R read counts
sum(F230R.Arth$ESVsize)
# 666,732

# Get Metazoa & ml-jg
mljg.Arth<-A[A$Phylum=="Arthropoda" &
                grepl("ml-jg", A$GlobalESV),]

# Sum Arthropoda & MiCOI ESV counts
length(unique(mljg.Arth$GlobalESV))
# 2,167

# Sum Arthropoda & MiCOI read counts
sum(mljg.Arth$ESVsize)
# 1,067,812
