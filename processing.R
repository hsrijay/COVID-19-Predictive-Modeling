# MESSI RNA Seq Analysis #

# Both batches

batch_type = 'all'

# Path for output #
```{r}
path <- "D:/Research/Summer-2020-+DS/data/"
save.image(file = "processing.RData")
```



raw_data_path <- paste0(path, "data/raw/")
readable_data_path <- paste0(path, "data/processed/", batch_type, "/readable/")
binary_data_path <- paste0(path, "data/processed/", batch_type, "/binary/")
figure_path <- paste0(path, "figures/", batch_type, "/")


# Libraries #
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v79))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(NMF))

# Number of NMF factors #
```{r}
num_factors = 10
```



# Import Key #
```{r}
path <- "D:/Research/Summer-2020-+DS/data/"
key <- read.csv(paste0(path,"raw_data/qry_covid_rnaseq_key_demog_analysis(cn 20200512).csv"), check.names=FALSE, stringsAsFactors = FALSE) #qry_20200424_covid_rnaseq_key_demog_sxsum_wFastq(cn 20200507).csv
key <- key[, c('Subject_ID', 'Timepoint', 'age_at_enroll', 'gender', 'race', 'RNA_ID', 'cohort', 'Pathogen', 'studyname')]
colnames(key) <- c('subject_id', 'timepoint', 'age', 'gender', 'race', 'rna_id', 'cohort', 'pathogen', 'studyname')
```



# Import batch 1 data #
```{r}
batch1 <- read.csv(paste0(path, "raw_data/covid_starprocessed_count_matrix_batch1.csv"), check.names=FALSE, stringsAsFactors = FALSE)
rownames(batch1) <- batch1[, 1]
batch1 <- batch1[, -1]
```



# Import batch 2 data #
```{r}
batch2 <- read.csv(paste0(path, "raw_data/covid_starprocessed_count_matrix_batch2.csv"), check.names=FALSE, stringsAsFactors = FALSE)
rownames(batch2) <- batch2[, 1]
batch2 <- batch2[, -1]
```



# samples in common across batches
```{r}
colnames(batch2)[which(colnames(batch2) %in% colnames(batch1))]
```


# key for all repeat samples
```{r}
key_all_repeats = key[which(key$rna_id %in% colnames(batch2)[which(colnames(batch2) %in% colnames(batch1))]), ]
```


# stratify repeat samples
```{r}
healthy_repeats = key_all_repeats$rna_id[which(key_all_repeats$pathogen == 'healthy')]
other_repeats = key_all_repeats$rna_id[which(key_all_repeats$pathogen != 'healthy')]
```



# Remove non-healthy repeats from batch 1
```{r}
dim(batch1)
batch1 <- batch1[, -which(colnames(batch1) %in% other_repeats)]
dim(batch1)
```


# Put healthy repeats at end of batch 1
```{r}
batch1_h <- batch1[, which(colnames(batch1) %in% healthy_repeats)]
batch1 <- batch1[, -which(colnames(batch1) %in% healthy_repeats)]
batch1 <- cbind(batch1, batch1_h)
```


# Put healthy repeats at end of batch 2
```{r}
batch2_h <- batch2[, which(colnames(batch2) %in% healthy_repeats)]
colnames(batch2_h) <- paste0(colnames(batch2_h), "_2")
batch2 <- batch2[, -which(colnames(batch2) %in% healthy_repeats)]
```



# check gene names are identical
```{r}
identical(rownames(batch1), rownames(batch2)) # should be TRUE
```


# counts matrix
```{r}
counts <- cbind(batch1, batch2)
```



# re-order key by subject and time point
```{r}
key <- key[order(as.factor(key$subject_id), as.factor(key$timepoint)), ]
```


# remove samples in data that are not in key
```{r}
match_order <- match(key$rna_id, colnames(counts))
```



# re-organize remaining subjects in data matrix to match those in the key
```{r}
counts <- counts[, match_order]
identical(key$rna_id, colnames(counts))
```



# Add back in healthy repeat samples in batch2 to data and key
```{r}
counts_sub <- cbind(counts, batch2_h)
key_h <- key[which(key$rna_id %in% healthy_repeats), ]
key_h$rna_id <- paste0(key_h$rna_id, "_2")
key_h <- key_h[match(colnames(batch2_h), paste0(healthy_repeats, "_2")), ]
key <- rbind(key, key_h)
rownames(key) <- 1:nrow(key)
key$batch <- rep(2, nrow(key))
key$batch[1:ncol(batch1)] <- 1
identical(key$rna_id, colnames(counts_sub))
```



# ## Formatting data ##
```{r}
library(edgeR)
#xpr <- DGEList(counts = counts_sub, samples = keyz)
xpr <- DGEList(counts = f, samples = gr_truth)
```

# 
# 
# # Transcript names #
```{r}
ref_seq <- rownames(xpr)
```

# 
# 
# # Annotation to map from transcripts to gene symbols #
```{r}
library(EnsDb.Hsapiens.v79)
gene_id <- AnnotationDbi::select(EnsDb.Hsapiens.v79,
                                  key=ref_seq,
                                  columns=c("SYMBOL"),
                                  keytype="GENEID")
```

# 
# 
# # Not all transcripts have a gene symbol: different lengths #

# length(ref_seq)
# nrow(gene_id)
# 
# # Transcripts that do not map to a gene symbol #
```{r}
genes_nosymb = unique(ref_seq[which(!(ref_seq %in% gene_id$GENEID))])
gene_nosymb_mat = data.frame(GENEID = genes_nosymb, SYMBOL = genes_nosymb)
```

# 
# 
# # Add back in genes for which a symbol could not be found #
```{r}
nrow(gene_nosymb_mat) == length(ref_seq) - length(gene_id$GENEID) # should be TRUE
gene_id = rbind(gene_id, gene_nosymb_mat)
nrow(gene_id) == length(ref_seq) # should be TRUE
```

# 
# 
# # Re-order gene annotation to match order of transcripts in data #
```{r}
idx_match = match(rownames(xpr), gene_id$GENEID)
gene_id = gene_id[idx_match, ]
identical(rownames(xpr), gene_id$GENEID) # should be true
```

# 
# # Re-name row names of data from transcripts to gene symbols #
```{r}
rownames(xpr) = gene_id$SYMBOL
identical(rownames(xpr), gene_id$SYMBOL)
```

# 
# 
# # Genes symbols that come up multiple times #
```{r}
dup_uniq_genes = unique(rownames(xpr)[which(duplicated(rownames(xpr)) == T)])
# # Keep track of gene symbols for row names of data #
gene_symbs = rownames(xpr)
# # Sum up reads for gene symbols that are the same #
for(i in 1:length(dup_uniq_genes)){
#   # which rows correspond to gene symbol
   dup_symb = which(gene_symbs == dup_uniq_genes[i])
#   # sum up reads for rows that correspond to same gene symbol
   new_row = data.frame(t(colSums(xpr$counts[dup_symb, ])))
   rownames(new_row) = dup_uniq_genes[i]
   colnames(new_row) = colnames(xpr)
#   # delete rows with the same gene symbol
   xpr$counts = xpr$counts[-dup_symb, ]
#   # create new row that contains the summed reads for that gene symbol
   xpr$counts = rbind(xpr$counts, new_row)
#   # update row names for data
   gene_symbs = c(gene_symbs[-dup_symb], dup_uniq_genes[i])
 }
```

# 
# # Assign row names for data with appropriate gene symbols #
```{r}
rownames(xpr) <- gene_symbs
length(unique(rownames(xpr))) == length(rownames(xpr)) # should be TRUE
```

# 
# 
# # Save / Read data #
```{r}
saveRDS(xpr, file=paste0(path, "xpr.rds"))
xpr <- readRDS(paste0(path, "xpr.rds"))
dim(xpr)
```



# Exclude non-influenza samples
```{r}
non_flu_samps = key$rna_id[which(key$cohort == 'Viral' & !(key$pathogen %in% c("Influenza  A 2009 H1N1", "Influenza A")))]
length(non_flu_samps)
xpr <- xpr[, -which(colnames(xpr) %in% non_flu_samps)]
key <- key[-which(key$rna_id %in% non_flu_samps), ]
```


# Replace 'Viral' category with 'Influenza'
```{r}
key$cohort[which(key$cohort == 'Viral')] <- 'Influenza'
dim(xpr)
dim(key)
```



### Quality Control ###

# Remove genes with no expression values #
```{r}
dim(xpr$counts)
xpr$counts <- xpr$counts[which(rowSums(xpr$counts) != 0), ]
dim(xpr$counts)
```



# Keep genes if they have higher than 1 count per million in at least half the samples #
```{r}
cpm <- cpm(xpr)
keep_xpr <- rowSums(cpm > 1) >= (ncol(xpr$counts) * (1/2)) # 1/2 the number of samples
xpr <- xpr[keep_xpr, , keep.lib.sizes=FALSE]
dim(xpr)
```



### Quality Control ###

# Number of observations for each subject #
```{r}
subj_table <- table(key$subject_id)
```



# Summary of counts #
```{r}
summary(colSums(xpr$counts)/1e6)
```



# Visualize number of reads per sample #
```{r}
library(tidyverse)
boxplot(colSums(xpr$counts)/1e6,  xlab = "Samples", ylab = "Million reads per sample"); grid()
plotData <- data.frame(counts = (colSums(xpr$counts)/1e6))
p <- ggplot(plotData, aes(x="", y=counts))
p <- p + geom_boxplot()
p <- p + theme_bw(18)
p <- p + labs(title="Million Reads per Sample", x="samples", y="million reads per sample")
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(path, "processing/exploratory/boxplot_counts.png"), width=10, height=6, dpi=300, units="in", device="png")
```



# Visualize number of reads per sample: sorted by number of reads #
```{r}
plot(sort(colSums(xpr$counts)/1e6), main = "Samples sorted by Number of Reads", xlab = "Samples", ylab = "Million reads per sample")
plotData <- data.frame(samps = seq(1, length(sort(colSums(xpr$counts)/1e6)), 1), counts = sort(colSums(xpr$counts)/1e6))
p <- ggplot(plotData, aes(x=samps, y=counts))
p <- p + geom_point()
p <- p + theme_bw(18)
p <- p + labs(title="Samples sorted by Number of Reads", x="samples", y="million reads per sample")
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(figure_path, "processing/exploratory/boxplot_sort_counts.png"), width=10, height=6, dpi=300, units="in", device="png")
```



# Correlation between samples: boxplot #
# Computes log counts per million #
```{r}
xpr_lcpm <- cpm(xpr, log=TRUE)
```


# Find pairwise correlation of samples #
```{r}
cor_mat = cor(xpr_lcpm, method="spearman")
```


```{r}
library(tidyverse)
boxplot(apply(cor_mat, 1, mean), main = "Sample Correlation", ylab = "Correlation")

plotData <- data.frame(corr = apply(cor_mat, 1, mean))
p <- ggplot(plotData, aes(x="", y=corr))
p <- p + geom_boxplot()
p <- p + theme_bw(18)
p <- p + labs(title="Sample Correlation", x="samples", y="correlation")
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
```


ggsave(paste0(figure_path, "processing/exploratory/boxplot_samp_cor.png"), width=10, height=6, dpi=300, units="in", device="png")
# Samples with bad quality: spearman correlation less than 0.9 #
```{r}
cor_rm <- colnames(cor_mat)[which(apply(cor_mat, 1, mean) < 0.75)]

```

```{r}
keys <- key
xprs <- xpr
key <- keys
xpr <- xprs
xpr$samples <- xpr$samples[!rownames(xpr$samples) %in% key_all_repeats, ] 
xpr$counts <- xpr$counts[ ,!colnames(xpr$counts) %in% key_all_repeats]
keyz <- keyz[which(!(keyz$rna_id %in% key_all_repeats)), ] 
x <- keyz %>% anti_join(key_all_repeats)

dim(keyz)
dim(xpr)
```


# Correlation between samples: heatmap #
# Computes log counts per million #
```{r}
xpr_lcpm <- cpm(xpr, log=TRUE)
# Find pairwise correlation of samples #
corr <- cor(xpr_lcpm, method="spearman")
```


# Plot correlation #
```{r}
colfunc <- colorRampPalette(c("lightblue", "darkblue"))

heatmap.2(corr, trace="none", col=colfunc(10))
title(main = "Correlation", cex.main=1.5)
dev.off()
```



# Density plot: before sample removal #
```{r}
xpr_lcpm <- cpm(xpr, log=TRUE)
nsamples <- ncol(xpr_lcpm)
#png(paste0(figure_path, "processing/exploratory/density_before.png"), width=9, height=8, units="in", res=300) #width=10, height=6, dpi=300, units="in"
plot(density(xpr_lcpm[, 1]), lwd=2, las=2, main="", xlab="", col=1)#, xlim=c(-10, 15))
title(xlab="Log CPM")
rm_samp = c()
for (i in 2:nsamples)
{
  den <- density(xpr_lcpm[,i])
  lines(den$x, den$y,lwd=2,col=i)
  if(max(den$y[den$x< (0)]) > 0.07) # determine threshold after looking at plot without if statement
   {
     print(colnames(xpr_lcpm)[i])
     rm_samp = c(rm_samp, colnames(xpr_lcpm)[i])
   }
}
dev.off()
```



# # Delete significantly bimodal samples from analysis #
```{r}
dim(xpr)
dim(key)
xpr$counts <- xpr$counts[ ,!colnames(xpr$counts) %in% rm_samp]
xpr$samples <- xpr$samples[!rownames(xpr$samples) %in% rm_samp, ] 
key <- key[which(!(key$rna_id %in% rm_samp)), ] 
dim(xpr)
dim(key)
```


# 
# # Density plot: after sample removal #
# xpr_lcpm <- cpm(xpr, log=TRUE)
# nsamples <- ncol(xpr_lcpm)
# png(paste(figure_path, "processing/exploratory/density_after.png", sep=""), width=9, height=8, units="in", res=300) #width=10, height=6, dpi=300, units="in"
# plot(density(xpr_lcpm[,1]), lwd=2, las=2, main="", xlab="",col=1)
# title(xlab="Log CPM")
# for (i in 2:nsamples)
# {
#   den <- density(xpr_lcpm[,i])
#   lines(den$x, den$y,lwd=2,col=i)
# }
# dev.off()

# Apply trimmed mean normalization to gene expression distributions #
```{r}
xpr_norm <- calcNormFactors(xpr, method = "TMM")
xpr_norm2 <- xpr_norm
xpr_norm2$samples$norm.factors <- 1
```



# Unnormalized data #
```{r}
lcpm2 <- cpm(xpr_norm2, log=TRUE)
```


# samples to remove based on normalization plot
```{r}
norm_rm = colnames(lcpm2)[order(apply(lcpm2, 2, mean), decreasing=T)][(length(colnames(lcpm2))-2):length(colnames(lcpm2))]

# plot
#png(paste0(figure_path, "processing/exploratory/unnormalized_data.png"), width=9, height=8, units="in", res=300) #width=10, height=6, dpi=300, units="in"
par(mar = c(5,4,4,1))
boxplot(x, las=2, main="", xaxt = 'n')
title(main = "A. Unnormalised data", ylab = "Log CPM", xlab = NULL)
dev.off()
```

```{r}
xpr$samples <- xpr$samples[!rownames(xpr$samples) %in% norm_rm, ] 
xpr$counts <- xpr$counts[ ,!colnames(xpr$counts) %in% norm_rm]
keyz <- keyz[which(!(keyz$rna_id %in% norm_rm)), ] 
dim(keyz)
dim(xpr)
```

# Normalized data #
```{r}
lcpm <- cpm(xpr_norm, log=TRUE)
png(paste0(figure_path, "processing/exploratory/normalized_data.png"), width=9, height=8, units="in", res=300) #width=10, height=6, dpi=300, units="in"
boxplot(lcpm, las=2, main="")
title(main = "B. Normalised data", ylab = "Log CPM")
dev.off()
```



## PCA on RNA Data ##
# Perform PCA on log cpm of normalized counts data #
```{r}
xpr_norm <- calcNormFactors(xpr, method = "TMM")
xpr_nlcpm <- cpm(xpr_norm, log=TRUE)
PCA_data <- prcomp(t(xpr_nlcpm), center=TRUE, scale=TRUE)
```



# PCA colored by batch #
```{r}
groupLabel = as.factor(xpr$samples$batch)
plotData <- data.frame(pc1=PCA_data$x[, 1], pc2=PCA_data$x[, 2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of RNA Data colored by Batch\nwith outliers")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
#ggsave(paste0(figure_path, "processing/exploratory/pca_batch_with_outliers.png"), width=10, height=6, dpi=300, units="in", device="png")
```



# Outliers
```{r}
PCA_data$x[which(PCA_data$x[, 1] < -100), 1]
# from batch 1, batch 1, and batch 2 respectively (by pca_batch_with_outliers.png)
pca_rm = names(PCA_data$x[which(PCA_data$x[, 1] < -100), 1])
```


# Remove samples with bad quality from analysis
```{r}
dim(xpr)
dim(key)
xpr <- xpr[ ,!colnames(xpr) %in% pca_rm]
key <- key[which(!(key$rna_id %in% pca_rm)), ]
dim(xpr)
dim(key)
```



# Perform PCA on log cpm of normalized counts data #
```{r}
xpr_norm <- calcNormFactors(xpr, method = "TMM")
xpr_nlcpm <- cpm(xpr_norm, log=TRUE)
PCA_data <- prcomp(t(xpr_nlcpm), center=TRUE, scale=TRUE)
```



# PCA colored by batch without outliers #
```{r}
groupLabel = as.factor(xpr$samples$batch)
plotData <- data.frame(pc1=PCA_data$x[, 1], pc2=PCA_data$x[, 2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of RNA Data \n colored by Batch")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
#ggsave(paste0(figure_path, "processing/exploratory/pca_batch_without_outliers.png"), width=10, height=6, dpi=300, units="in", device="png")
```



# PCA colored by batch and healthy re-reruns #
```{r}
groupLabel = xpr$samples$batch
groupLabel[which(xpr$samples$rna_id %in% healthy_repeats)] = -1
groupLabel[which(xpr$samples$rna_id %in% paste0(healthy_repeats, "_2"))] = -2
groupLabel = as.factor(groupLabel)
plotData <- data.frame(pc1=PCA_data$x[, 1], pc2=PCA_data$x[, 2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + scale_color_manual(values=c("red", "darkblue", "lightblue", "pink"))
p <- p + labs(title="PCA of RNA Data colored by Batch\nshowing repeated healthy samples")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
#ggsave(paste0(figure_path, "processing/exploratory/pca_batch_with_healthy.png"), width=10, height=6, dpi=300, units="in", device="png")
```



# Remove repeat healthy samples from first batch
```{r}
dim(xpr)
dim(key)
xpr <- xpr[ ,!colnames(xpr) %in% healthy_repeats]
key <- key[which(!(key$rna_id %in% healthy_repeats)), ]
dim(xpr)
dim(key)

rm <- c("DU18-02S0011611")
#remove covid sample
key <- key[which(!(key$rna_id %in% rm)),]
xpr <- xpr[,!colnames(xpr) %in% rm]
dim(xpr)
dim(key)

#remove na genders
na <- which(!(key$gender != 'MALE' & key$gender != 'FEMALE'))
b <- key[na,]

key <- key[which(!(key$gender != 'MALE' & key$gender != 'FEMALE')),]
xpr <- xpr[,!colnames(xpr) %in% b$rna_id]
dim(key)
dim(xpr)

```



# Perform PCA on log cpm of normalized counts data #
xpr_norm <- calcNormFactors(xpr, method = "TMM")
xpr_nlcpm <- cpm(xpr_norm, log=TRUE)
PCA_data <- prcomp(t(xpr_nlcpm), center=TRUE, scale=TRUE)

# PCA colored by time point #
xpr$samples$timepoint[which(xpr$samples$timepoint == 'DAY 7')] = 'Day 7'
xpr$samples$timepoint[which(xpr$samples$timepoint == 'DAY 14')] = 'Day 14'
groupLabel <- factor(xpr$samples$timepoint, levels = c('T=0', 'Day 1', 'Day 2', 'Day 3', 'Day 7', 'Day 14', 'Day 21'))
plotData <- data.frame(pc1=PCA_data$x[which(xpr$samples$cohort=='COVID-19'),1], pc2=PCA_data$x[which(xpr$samples$cohort=='COVID-19'),2], group=groupLabel[which(xpr$samples$cohort=='COVID-19')])
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of COVID-19 Data \n colored by Time Point")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(figure_path, "processing/exploratory/pca_covid19_timepoints.png"), width=10, height=6, dpi=300, units="in", device="png")

# PCA colored by age #
groupLabel <- cut(xpr$samples$age, breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
plotData <- data.frame(pc1=PCA_data$x[,1], pc2=PCA_data$x[,2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of RNA Data \n colored by Age")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(figure_path, "processing/exploratory/pca_age.png"), width=10, height=6, dpi=300, units="in", device="png")

# PCA colored by gender #
xpr$samples$gender[which(xpr$samples$gender == 'FEMALE')] = 'Female'
xpr$samples$gender[which(xpr$samples$gender == 'MALE')] = 'Male'
xpr$samples$gender[which(xpr$samples$gender == "")] = NA
groupLabel <- as.factor(xpr$samples$gender)
plotData <- data.frame(pc1=PCA_data$x[,1], pc2=PCA_data$x[,2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of RNA Data \n colored by Gender")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(figure_path, "processing/exploratory/pca_gender.png"), width=10, height=6, dpi=300, units="in", device="png")

# PCA colored by race #
xpr$samples$race[which(xpr$samples$race == '\tAsian')] = 'Asian'
xpr$samples$race[which(xpr$samples$race == 'WHITE')] = 'White'
xpr$samples$race[which(xpr$samples$race == 'BLACK/AFRICAN AMERICAN')] = 'Black/African American'
xpr$samples$race[which(xpr$samples$race == 'UNKNOWN')] = 'Unknown/Not reported'
xpr$samples$race[which(xpr$samples$race == "")] = NA
groupLabel <- as.factor(xpr$samples$race)
plotData <- data.frame(pc1=PCA_data$x[,1], pc2=PCA_data$x[,2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of RNA Data \n colored by Race")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(figure_path, "processing/exploratory/pca_race.png"), width=10, height=6, dpi=300, units="in", device="png")

# PCA colored by cohort #
groupLabel <- as.factor(xpr$samples$cohort)
plotData <- data.frame(pc1=PCA_data$x[,1], pc2=PCA_data$x[,2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of RNA Data \n colored by Cohort")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(figure_path, "processing/exploratory/pca_cohort.png"), width=10, height=6, dpi=300, units="in", device="png")

# PCA colored by pathogen #
groupLabel <- as.factor(xpr$samples$pathogen)
plotData <- data.frame(pc1=PCA_data$x[,1], pc2=PCA_data$x[,2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of RNA Data \n colored by Pathogen")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(figure_path, "processing/exploratory/pca_pathogen.png"), width=10, height=6, dpi=300, units="in", device="png")

# PCA colored by study #
xpr$samples$studyname[which(xpr$samples$studyname == 'ACESO/ARLG (DU09-03 continuation)')] = 'ACESO/ARLG'
groupLabel <- as.factor(xpr$samples$studyname)
plotData <- data.frame(pc1=PCA_data$x[,1], pc2=PCA_data$x[,2], group=groupLabel)
p <- ggplot(plotData, aes(x=pc1, y=pc2))
p <- p + geom_point(aes(colour=group), size=4)
p <- p + labs(title="PCA of RNA Data \n colored by Study")
p <- p + theme_bw(18)
p <- p + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
p
ggsave(paste0(figure_path, "processing/exploratory/pca_study.png"), width=10, height=6, dpi=300, units="in", device="png")

### Normalization ###

# Apply trimmed mean normalization to data #
```{r}
xpr_tmm <- calcNormFactors(xpr, method = "TMM")
xpr_nlcpm <- cpm(xpr_tmm, log = TRUE)

#write.csv(t(rna_nmf_nlcpm), file=paste0(path, "195/processed_data/nmf/rna_nmf_nlcpm_sex_33.csv"))
```



### Non-Negative Matrix Factorization ###
```{r fig.height = 6, fig.width = 15}
# # NMF on TMM data #
rna_nmf_tmm <- nmf(xpr_tmm$counts, rank=33, seed=123456, .options='t')
#Convert nmf data to log-cpm
rna_nmf_nlcpm_sex <- cpm(rna_nmf_tmm@fit@H, log=TRUE)
# rank estimation
estim.r <- nmf(xpr_tmm$counts, 31:40, nrun = 10, seed = 123456)
estim.r <- nmf(xpr_tmm, 31:40, nrun = 10, seed = 123456)

plot(estim.r)
consensusmap(estim.r, labCol = NA, labRow = NA)
# # Track loss #
figure_path <- "D:/Research/Summer-2020-+DS/output"
png(paste0(figure_path, "rna_nmf_tmm_residuals_20.png"))
plot(rna_nmf_tmm)
dev.off()
# # Save NMF object
saveRDS(rna_nmf_tmm, paste0(path, '196/processed_data/rna_nmf_tmm_33.rds'))
saveRDS(estim.r, paste0(path, '/processed_data/rank_estimation_10-30_10run.rds'))
#plot visualizations
layout(cbind(1,3))
basismap(rna_nmf_tmm, subsetRow = TRUE)
coefmap(rna_nmf_tmm)
consensusmap(rna_nmf_tmm, labCol = NA, labRow = NA)
```



# Save data objects #
```{r}
dim(xpr_nlcpm)
dim(key)
write.csv(t(xpr_nlcpm), file=paste0(path, "195/processed_data/processing/xpr_nlcpm.csv"))
write.csv(key, file=paste0(path, "196/processed_data/processing/key.csv"), row.names=F)
```



# Samples should be in same order
```{r}
identical(colnames(xpr), colnames(xpr_tmm)) # should be TRUE
identical(colnames(xpr_nlcpm), colnames(xpr_tmm)) # should be TRUE
identical(colnames(xpr), key$rna_id) # should be TRUE


```


# Write to file
```{r}
write.csv(t(xpr), paste0(path, "195/processed_data/processing/xpr.csv"))
write.csv(t(xpr_tmm$counts), paste0(path, "195/processed_data/processing/sex/xpr_tmm_sex.csv"))
write.csv(t(rna_nmf_tmm@fit@H), paste0(path, "195/processed_data/nmf/sex/rna_nmf_tmm_sex_33.csv"))
key <- read.csv(file = "D:/Research/Summer-2020-+DS/data/196/processed_data/processing/key.csv")
key_subset <- read.csv(file = "D:/Research/Summer-2020-+DS/data/196/processed_data/processing/key_subset.csv")
xpr_nlcpm <- xpr_nlcpm %>%
   filter(X != "DU18-02S0011611")
```
