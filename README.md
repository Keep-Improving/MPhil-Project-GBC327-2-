# MPhil-Project-GBC327-2-
MPhil project's R scripts
# Day 1 -------------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
library(tximport)
install.packages("locfit")
library(DESeq2)
library(tidyverse)
library(readr)
library(dplyr)

#2. Data import-----
#Reading in the sample metadata
sampleinfo <- read_tsv("D:/朱泓宇/MPhil in Obesity, Endocrinology and Metabolism/Project/Methods/Transcriptomics/Bulk RNA-Seq Cambridge Bioinformatic Training Program/data/samplesheet.tsv", col_types = c("cccc"))
arrange(sampleinfo, Status, TimePoint, Replicate)

#Reading in the count data
files <- file.path("salmon",sampleinfo$SampleName,"quant.sf") #have the path to all the samples' counts ?
files

library(rlang)
files <- set_names(files, sampleinfo$SampleName) # link the sample names to the transripts
files
#then need to link the transcripts to the gene names (gene IDs)

tx2gene <- read_tsv("D:/朱泓宇/MPhil in Obesity, Endocrinology and Metabolism/Project/Methods/Transcriptomics/Bulk RNA-Seq Cambridge Bioinformatic Training Program/references/tx2gene.tsv") # make a table of gene names

txi<- tximport(files, type = "salmon", tx2gene=tx2gene) #tximport: import the matrix data; tx2gene: the geneID
str(txi)

head(txi$counts)

#save the txi file by saveRDS()
dir.create("salmon_outputs") #first create the directory
saveRDS(txi, file = "salmon_outputs/txi.rds") 



#3. Prepare count matrix-----
#Create a raw counts matrix for data exploration
rawCounts <- round(txi$counts, 0) #round the count numbers
head(rawCounts)

#filter out low-count genes (which maybe noise, don't need them)
dim(rawCounts)

keep <- rowSums(rawCounts)>5 #here keep genes with counts >5
table(keep, useNA = "always")

filtCounts <- rawCounts[keep,]
dim(filtCounts)


#4. Count distribution and Data transformations-----
#log transformation (for later PCA and differential analysis)
summary(filtCounts)
boxplot(filtCounts, main='Raw counts',las = 2) #visualize the raw count distribution

plot(rowMeans(filtCounts),rowSds(filtCounts), # visualize the relationship between raw count mean and raw count SD
     main = 'Raw counts: sd vs mean',
     xlim= c(0,10000),
     ylim= c(0,5000)) 

#do the log transformation by easy log2 function; but this doesn't consider the variance-mean relationship
logCounts <- log2(filtCounts + 1)
boxplot(logCounts, main = 'Log2 counts', las = 2)

#therefore, use the DESeq2 library function: VST or rlog, to solve the relationship
rlogcounts <- rlog(filtCounts) #rlog is better for small-size data ; convert data to integer mode
boxplot(rlogcounts, main = 'rlog counts', las = 2)




#5 PCA -----
install.packages("purrr")
library(ggfortify)

rlogcounts <- rlog(filtCounts)

# run PCA
pcDat <- prcomp(t(rlogcounts)) #t: to separate samples, not genes
# plot PCA
autoplot(pcDat)

autoplot(pcDat,
         data = sampleinfo, 
         colour = "Status", 
         shape = "TimePoint",
         size = 5)


#Exercise: PCA showing the PC3 on x axis and PC2 on y axis
autoplot(pcDat,
         data = sampleinfo, 
         colour = "Status", 
         shape = "TimePoint",
         x=3,
         y=2,
         size = 5)


#found a wired sample in the opposite cluster --> its Status was found wrongly labelled

library(ggrepel) # add labels

autoplot(pcDat,
         data = sampleinfo,  
         colour = "Status", 
         shape = "TimePoint", # setting shape to FALSE causes the plot to default to using the labels instead of points
         size = 5) +
  geom_text_repel(aes(x = PC1, y = PC2, label = SampleName), box.padding = 0.8)

#change the status label
sampleinfo <- mutate(sampleinfo, Status = case_when(
  SampleName=="SRR7657882" ~ "Uninfected",
  SampleName=="SRR7657873" ~ "Infected", 
  TRUE ~ Status))

#output
dir.create("results")
write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")

#check PCA again
autoplot(pcDat,
         data = sampleinfo, 
         colour = "Status", 
         shape = "TimePoint",
         size = 5)



#Clustering
library(ggdendro)
hclDat <-  t(rlogcounts) %>%
  dist(method = "euclidean") %>%
  hclust()
ggdendrogram(hclDat, rotate = TRUE)

hclDat2 <- hclDat
hclDat2$labels <- str_c(sampleinfo$Status, ":", sampleinfo$TimePoint) #show the factors
ggdendrogram(hclDat2, rotate = TRUE)
