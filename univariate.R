# Libraries #
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(xtable))
suppressPackageStartupMessages(library(Vennerable))

# Set seed #
```{r}
set.seed(1337)
```



# Batch 1 or both batches
batch_type = 'all' # 'batch1', 'all'
# Number of top genes to plot #
```{r}
num_genes = 4
pan_num_genes <- 2
```


# Paths #
```{r}
path <- "D:/Research/Summer-2020-+DS/data/processed_data/"
code_path <- paste0(path, "code/")
readable_data_path_readin <- paste0(path, "data/processed/", batch_type, "/readable/")
readable_data_path <- "D:/Research/Summer-2020-+DS/data/"
binary_data_path <- paste0(path, "data/processed/", batch_type, "/binary/")
figure_path <- paste0(path, "figures/", batch_type, "/")
report_path <- paste0(path, "writing/univariate/")
```



# Data #
```{r}
path <- "D:/Research/Summer-2020-+DS/data/195/processed_data/processing/"
xpr_tmm <- read.csv(paste0(path, "xpr_tmm.csv"))
xpr_tmm <- data.frame(xpr_tmm[, -1], row.names = xpr_tmm[, 1])
xpr_tmm <- as.matrix(xpr_tmm)

xpr_nlcpm <- read.csv(paste0(path, "xpr_nlcpm.csv"))
xpr_nlcpm <- data.frame(xpr_nlcpm[, -1], row.names = xpr_nlcpm[, 1])
xpr_nlcpm <- as.matrix(xpr_nlcpm)
```



# Pan Viral Signature #
```{r}
pan_viral_genes <- c('ATF3', 'CCL2', 'DDX58', 'DECR1', 'FARP1', 'FPGS', 'GAPDH', 'GBP1', 'HERC5', 'IFI27', 'IFI44', 'IFI44L', 'IFI6', 'IFIT1', 'IFIT2', 'IFIT3', 'IFIT5', 'ISG15', 'LAMP3', 'LY6E', 'MX1', 'OAS1', 'OAS2', 'OAS3', 'OASL', 'PPIA', 'PPIB', 'RPL30', 'RSAD2', 'RTP4', 'SEPT4', 'SERPING1', 'SIGLEC1', 'TNFAIP6', 'TRAP1', 'XAF1')
pan_viral_genes <- colnames(xpr_nlcpm)[colnames(xpr_nlcpm) %in% pan_viral_genes]
pan_viral_genes <- colnames(rna_nmf_nlcpm)[colnames(rna_nmf_nlcpm) %in% pan_viral_genes]
```



# House keeping genes #
```{r}
hk_genes <- c('TRAP1', 'FPGS', 'GAPDH', 'DECR1', 'PPIB', 'PPIA', 'RPL30', 'FARP1')
```



# Key #
```{r}
key <- read.csv(paste0(path, "key.csv"))
key$subject_id <- as.factor(key$subject_id)
key$cohort <- factor(key$cohort, levels = c('Bacterial', 'COVID-19', 'CoV other', 'Influenza', 'healthy'))




```

# #Timepoint Analysis
# ```{r}
# key_subset <- read.csv(paste0(readable_data_path, "raw_data/key_subset.csv"))
# key_subset$subject_id <- as.factor(key_subset$subject_id)
# key_subset$cohort <- factor(key_subset$cohort, levels = c('COVID-19', 'CoV other', 'Influenza', 'Bacterial', 'healthy'))
# covid_samps <- key_subset %>%
#   dplyr::filter(cohort == "COVID-19")
# covid <- key_subset$rna_id[which(key_subset$cohort == "COVID-19")]
# covid_nmf <- rna_nmf_nlcpm[rownames(rna_nmf_nlcpm) %in% covid,] 
# reorder_idx <- match(covid_samps$rna_id, rownames(covid_nmf))
# covid_nmf <- covid_nmf[rownames(covid_nmf)[reorder_idx],]
# identical(rownames(covid_nmf), covid_samps$rna_id)
# covid_nmf <- cbind(covid_nmf, days_since_onset = covid_samps$days_since_onset)
# covid_nmf <- cbind(covid_nmf, subject_id = covid_time$subject_id)
# covid_nmf <- as.data.frame(covid_nmf)
# covid_nmf$subject_id <- factor(covid_nmf$subject_id)
# 
# 
# a <- ggplot(data = covid_nmf, aes(x = as.numeric(days_since_onset), y = as.numeric(V4))) +
#   geom_point() +
#   geom_smooth(method = "lm", alpha = 0.5)+
#   theme_minimal()+
#   labs(y = "Factor Score", x = "Days Since Onset", title = "Top Ranked Feature")
# b <- ggplot(data = covid_nmf, aes(x = as.numeric(days_since_onset), y = as.numeric(V7))) +
#   geom_point() +
#   geom_smooth(method = "lm", alpha = 0.5)+
#   theme_minimal()+
#   labs(y = "Factor Score", x = "Days Since Onset", title = "Second Ranked Feature")
# c <- ggplot(data = covid_nmf, aes(x = as.numeric(days_since_onset), y = as.numeric(V16))) +
#   geom_point() +
#   geom_smooth(method = "lm", alpha = 0.5)+
#   theme_minimal()+
#   labs(y = "Factor Score", x = "Days Since Onset", title = "Third Ranked Feature")
# d <- ggplot(data = covid_nmf, aes(x = as.numeric(days_since_onset), y = as.numeric(V19))) +
#   geom_point() +
#   geom_smooth(method = "lm", alpha = 0.5)+
#   theme_minimal()+
#   labs(y = "Factor Score", x = "Days Since Onset", title = "Fourth Ranked Feature")
# 
# ggarrange(a, b, c, d,
#                     labels = c("1.", "2.", "3.", "4."),
#                     ncol = 2, nrow = 2)
# 
# ```

##Timepoint Sex Analysis
# 
# ```{r}
# 
# rm <- c("DU18-02S0011654")
# sex_samps <- key_subset %>%
#   dplyr:: filter(cohort == 'COVID-19') %>%
#   dplyr:: filter(rna_id != "DU18-02S0011654")
# 
# keyz_subset <- key_subset %>%
#   dplyr::filter(rna_id != "DU18-02S0011654")
# #identical(rownames(covid_nmf), sex_samps$rna_id)
# #sex_nmf <- cbind(covid_nmf, gender = sex_samps$gender)
# 
# sex <- sex_samps$rna_id
# sex_nmf <- rna_nmf_nlcpm_sex[rownames(rna_nmf_nlcpm_sex) %in% sex,] 
# reorder_idx <- match(sex_samps$rna_id, rownames(sex_nmf))
# sex_nmf <- sex_nmf[rownames(sex_nmf)[reorder_idx],]
# identical(rownames(sex_nmf), sex_samps$rna_id)
# sex_nmf <- cbind(sex_nmf, days_since_onset = sex_samps$days_since_onset)
# #covid_nmf <- cbind(covid_nmf, subject_id = covid_time$subject_id)
# sex_nmf <- as.data.frame(sex_nmf)
# #covid_nmf$subject_id <- factor(covid_nmf$subject_id)
# 
# a <- ggplot(data = sex_nmf, aes(x = as.numeric(days_since_onset), y = as.numeric(V11))) +
#   geom_point() +
#   geom_smooth(method = "lm", alpha = 0.5)+
#   theme_minimal()+
#   labs(y = "Factor Score", x = "Days Since Onset", title = "Top Ranked Feature", subtitle = "Most Differentially Expressed Between Sexes")
# b <- ggplot(data = sex_nmf, aes(x = as.numeric(days_since_onset), y = as.numeric(V13))) +
#   geom_point() +
#   geom_smooth(method = "lm", alpha = 0.5)+
#   theme_minimal()+
#   labs(y = "Factor Score", x = "Days Since Onset", title = "Second Ranked Feature")
# c <- ggplot(data = sex_nmf, aes(x = as.numeric(days_since_onset), y = as.numeric(V29))) +
#   geom_point() +
#   geom_smooth(method = "lm", alpha = 0.5)+
#   theme_minimal()+
#   labs(y = "Factor Score", x = "Days Since Onset", title = "Third Ranked Feature")
# d <- ggplot(data = sex_nmf, aes(x = as.numeric(days_since_onset), y = as.numeric(V26))) +
#   geom_point() +
#   geom_smooth(method = "lm", alpha = 0.5)+
#   theme_minimal()+
#   labs(y = "Factor Score", x = "Days Since Onset", title = "Fourth Ranked Feature")
# 
# ggarrange(a, b, c, d,
#           labels = c("1.", "2.", "3.", "4."),
#           ncol = 2, nrow = 2, title = "Most Differentially Expressed")
# ```
# 
# 
# # Limma
# ```{r}
# key_sex <- key %>%
#   dplyr::filter(gender == "MALE" | gender == "FEMALE")
# 
# 
# key_sex$gender <- factor(key_sex$gender, levels = c("MALE", "FEMALE"))
# design <- model.matrix(~ key_sex$gender, data = key_sex)
# voomObj <- voom(t(rna_nmf_tmm_sex), design, plot=FALSE)
# dupcor <- duplicateCorrelation(voomObj, design, block=key_sex$subject_id)
# voomObj = voom(t(rna_nmf_tmm_sex), design, plot=FALSE, block=key_sex$subject_id, correlation=dupcor$consensus)
# dupcor <- duplicateCorrelation(voomObj, design, block=key_sex$subject_id)
# fitDupCor <- eBayes(lmFit(voomObj, design, block=key_sex$subject_id, correlation=dupcor$consensus))
# topTable(fitDupCor, coef = 2, number = Inf, sort.by = "P")
# ```
# 
# #Sex Boxplot
# ```{r}
# top_sex_boxplot <- function(rna_nmf_sex, key, top_gene, gene_rank){
#   top_box = data.frame(gene = rna_nmf_sex[, top_gene], sex = key$gender)
#   p = ggplot(top_box, aes(x=sex, y=gene, group=sex, fill=sex))
#   p = p + geom_boxplot()
#   p = p + scale_fill_manual(labels = levels(key$gender), values=c("lightblue", "pink"), guide=FALSE)
#   p = p + geom_jitter(aes(colour=sex), position = position_dodge(0.75), size=2)
#   p = p + scale_color_manual(labels = levels(key$sex), values = c("blue", "pink"))
#   if(gene_rank == 'NA'){
#     p = p + labs(title=paste0("Factor Score for ", top_gene), y="Factor Score", x="Sex")
#   }else{
#     p = p + labs(title=paste0("Factor Score for ", top_gene, "\n ranked ", gene_rank, " out of ", ncol(rna_nmf_sex)), y="Factor Score", x="Sex")
#   }
#   p = p + theme_bw(18)
#   p = p + theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=90))
#   p = p + theme(axis.title.x = element_text(size = 25, face="bold"), axis.title.y = element_text(size = 25, face="bold"))
#   p = p + theme(plot.title = element_text(size=30, face="bold", hjust=0.5))
#   p = p+ theme(legend.position = "none")
#   p
#   return(p)
# }
# ```
# 
# 
# #Generate plots
# ```{r}
# figure_path <- "D:/Research/Summer-2020-+DS/output/195/univariate_output/nmf_33/male_vs_female/"
# cohort_filename = c('cov_other', 'influenza', 'bacterial', 'healthy')
# compare_type = c('Other Coronavirus', 'Influenza', 'Bacterial', 'Healthy')
# vol_sig_genes = list()
# for(i in 1:length(cohort_filename)){
#   print(i)
#   topper <- topTable(fitDupCor, coef = 2, number = Inf, sort.by = "P")
#   write.csv(topper, paste0(figure_path, 'topper_nmf_sex.csv'))
#   topper_15 <- topTable(fitDupCor, coef = i+1, number = 15, sort.by = "P")
#   #print.xtable(xtable(topper_15), type = "latex", file = paste0(report_path, 'univariate_covid_other.tex'))
#   # volcano_plot = volcano_results(fitDupCor, paste0(': Male vs Female'), 2)$vol_plot
#   # vol_sig_genes[[1]] = volcano_results(fitDupCor, paste0(': Male vs Female'), 2)$sig_genes
#   # ggsave(paste0(figure_path, 'volcano_sex.png'), width=10, height=6, dpi=300, units="in", device="png")
#   for(k in 1:4){
#     top_gene_boxplot(rna_nmf_nlcpm_sex, key_sex, rownames(topper)[k], k)
#     ggsave(paste0(figure_path, 'gene', k, '_other.png'), width=10, height=6, dpi=300, units="in", device="png")
#   }
# }
# 
# 
# top_box = data.frame(gene = rna_nmf[, "V21"], sex = key_sex$gender)
# ggplot(top_box, aes(x=sex, y=gene, group=sex, fill=sex))+
#   geom_boxplot()+
#   scale_fill_manual(labels = levels(key$sex), values=c("lightblue", "pink"), guide=FALSE) +
#   theme_minimal() 
# 
# ```


# Volcano Plot #
```{r}
volcano_results <- function(fit, compare_type, cf)
{
  ## Volcano Plot ##
  # Results from Limma #
  top_all <- topTable(fit, coef = cf, adjust.method="BH", n=Inf)
  # Color for non-significant genes #
  color <- rep("gray", nrow(top_all))
  # Color for significant genes #
  color[top_all$adj.P.Val < 0.05] <- "red"
  # Volcano plot creation #
  vol_plot <- ggplot(top_all, aes(x=logFC, y=-log10(adj.P.Val)))
  vol_plot <- vol_plot + geom_point(color=color)
  vol_plot <- vol_plot + labs(title=paste0("Volcano plot", compare_type, "\nNumber of significant genes: ", sum(top_all$adj.P.Val < 0.05), " of ", nrow(top_all)), y="-log10(p-value)")
  vol_plot <- vol_plot + theme_bw(18)
  vol_plot <- vol_plot + theme(plot.title = element_text(size = 18, face="bold", hjust=0.5))
  return(list(vol_plot=vol_plot, sig_genes=rownames(top_all)[which(top_all$adj.P.Val < 0.05)]))
}
```



# Gene Expression by Cohort Box plot #
```{r}
top_gene_boxplot <- function(xpr_nlcpm, key, top_gene, gene_rank){
  top_box = data.frame(gene = xpr_nlcpm[, top_gene], cohort = key$cohort)
  p = ggplot(top_box, aes(x=cohort, y=gene, group=cohort, fill=cohort))
  p = p + geom_boxplot()
  p = p + scale_fill_manual(labels = levels(key$cohort), values=c("pink", "peachpuff", "lightgoldenrodyellow", "lightgreen", "lightblue"), guide=FALSE)
  p = p + geom_jitter(aes(colour=cohort), position = position_dodge(0.75), size=2)
  p = p + scale_color_manual(labels = levels(key$cohort), values = c("red", "orange", "gold", "darkgreen", "blue"))
  if(gene_rank == 'NA'){
    p = p + labs(title=paste0("Gene Expression for ", top_gene), y="Expression", x="Cohort")
  }else{
    p = p + labs(title=paste0("Gene Expression for ", top_gene, "\n ranked ", gene_rank, " out of ", ncol(xpr_nlcpm)), y="Expression", x="Cohort")
  }
  p = p + theme_bw(18)
  p = p + theme(axis.text.y = element_text(size=20), axis.text.x = element_text(size=20, angle=90))
  p = p + theme(axis.title.x = element_text(size = 25, face="bold"), axis.title.y = element_text(size = 25, face="bold"))
  p = p + theme(plot.title = element_text(size=30, face="bold", hjust=0.5))
  p
  return(p)
}
```

```{r}
save.image(file = "nmf_26_boxplots.RData")
```


# Limma #
```{r}
design <- model.matrix(~ key$cohort, data = key)
voomObj <- voom(t(xpr_tmm), design, plot=FALSE)
dupcor <- duplicateCorrelation(voomObj, design, block=key$subject_id)
voomObj = voom(t(xpr_tmm), design, plot=TRUE, block=key$subject_id, correlation=dupcor$consensus)
dupcor <- duplicateCorrelation(voomObj, design, block=key$subject_id)
fitDupCor <- eBayes(lmFit(voomObj, design, block=key$subject_id, correlation=dupcor$consensus))
saveRDS(fitDupCor, paste0(path, 'limma_fit.rds'))
```

# Limma NMF #
```{r}
design <- model.matrix(~ key$cohort, data = key)

path <- "D:/Research/Summer-2020-+DS/data/"
rna_nmf_tmm <- read.csv(paste0(path, "195/processed_data/nmf/rna_nmf_tmm_33.csv"))
rna_nmf_tmm <- data.frame(rna_nmf_tmm[, -1], row.names = rna_nmf_tmm[, 1])
rna_nmf_tmm <- as.matrix(rna_nmf_tmm)

rna_nmf_sex <- read.csv(file = paste0(readable_data_path, "195/processed_data/nmf/sex/rna_nmf_tmm_sex_33.csv"))
rna_nmf_sex <- data.frame(rna_nmf_sex[, -1], row.names = rna_nmf_sex[, 1])
rna_nmf_sex <- as.matrix(rna_nmf_sex)
rna_nmf_tmm_sex <- rna_nmf_sex



voomObj <- voom(t(rna_nmf_tmm), design, plot=FALSE)
dupcor <- duplicateCorrelation(voomObj, design, block=key$subject_id)
voomObj = voom(t(rna_nmf_tmm), design, plot=FALSE, block=key$subject_id, correlation=dupcor$consensus)
dupcor <- duplicateCorrelation(voomObj, design, block=key$subject_id)
fitDupCor <- eBayes(lmFit(voomObj, design, block=key$subject_id, correlation=dupcor$consensus))
saveRDS(fitDupCor, paste0(path, 'processed_data/limma_nmf_fit.rds'))
```



```{r}

figure_path <- "D:/Research/Summer-2020-+DS/output/195/univariate_output/non_nmf/cov_vs_all_other/"
#fitDupCor <- readRDS(paste0(path, 'processed_data/limma_nmf_fit.rds'))
cohort_filename = c('cov_other', 'influenza', 'bacterial', 'healthy')
compare_type = c('Other Coronavirus', 'Influenza', 'Bacterial', 'Healthy')
vol_sig_genes = list()
 for(i in 1:length(cohort_filename)){
   print(i)
    topper <- topTable(fitDupCor, coef = i+1, number = Inf, sort.by = "P")
   write.csv(topper, paste0(path, 'univariate_covid_', cohort_filename[i], '.csv'))
 }
  topper_15 <- topTable(fitDupCor, coef = i+1, number = 15, sort.by = "P")
  # print.xtable(xtable(topper_15), type = "latex", file = paste0(path, 'univariate_covid_', cohort_filename[i], '.tex'))
  volcano_plot = volcano_results(fitDupCor, paste0(': COVID-19 vs ', compare_type[i]), i+1)$vol_plot
  vol_sig_genes[[i]] = volcano_results(fitDupCor, paste0(': COVID-19 vs ', compare_type[i]), i+1)$sig_genes
  #
  ggsave(paste0(figure_path, 'nmf_volcano_', cohort_filename[i], '.png'), width=10, height=6, dpi=300, units="in", device="png")
 

 for(k in 1:4){
   top_gene_boxplot(rna_nmf_nlcpm, key, rownames(topper)[k], k)
   ggsave(paste0(figure_path, 'box_', k, '_', cohort_filename[i], '.png'), width=10, height=6, dpi=300, units="in", device="png")
 }
 # for(m in 1:pan_num_genes){
 # top_gene_boxplot(rna_nmf_nlcpm, key, rownames(topper)[which(rownames(topper) %in% pan_viral_genes)[m]], which(rownames(topper) %in% pan_viral_genes)[m])
 # ggsave(paste0(figure_path, '/figures/boxplots/', m, '_', cohort_filename[i], '.png'), width=10, height=6, dpi=300, units="in", device="png")
 # }
}
```

```{r}
save.image(file = "univariate_nmf_boxplots.RData")
```

```{r}
top_gene_boxplot(xpr_nlcpm, key, 'RP11.34P13.7', 1)
```



# COVID19 vs All Others #
```{r}
key_tmp <- key
key_tmp$healthy_other <- as.character(key_tmp$cohort)
key_tmp$healthy_other[which(key_tmp$cohort != 'healthy')] <- 'Other'
key_tmp$healthy_other <- factor(key_tmp$healthy_other, levels=c('healthy', 'Other'))

```



# Limma #
```{r}
design <- model.matrix(~ key_tmp$healthy_other, data = key_tmp)
voomObj <- voom(t(xpr_tmm), design, plot=FALSE)
dupcor <- duplicateCorrelation(voomObj, design, block=key_tmp$subject_id)
voomObj = voom(t(xpr_tmm), design, plot=FALSE, block=key_tmp$subject_id, correlation=dupcor$consensus)
dupcor <- duplicateCorrelation(voomObj, design, block=key_tmp$subject_id)
fitDupCor <- eBayes(lmFit(voomObj, design, block=key_tmp$subject_id, correlation=dupcor$consensus))
```


# saveRDS(fitDupCor, paste0(binary_data_path, 'limma_fit_covid_other.rds'))
fitDupCor <- readRDS(paste0(binary_data_path, 'limma_fit_covid_other.rds'))

```{r}
figure_path <- "D:/Research/Summer-2020-+DS/output/195/univariate_output/nmf_33/cov_vs_other/"
cohort_filename = c('cov_other', 'influenza', 'bacterial', 'healthy')
compare_type = c('Other Coronavirus', 'Influenza', 'Bacterial', 'Healthy')
vol_sig_genes = list()
topper <- topTable(fitDupCor, coef = 2, number = Inf, sort.by = "P")
write.csv(topper, paste0(path, 'topper_healthy_other.csv'))
topper_15 <- topTable(fitDupCor, coef = 2, number = 15, sort.by = "P")
#print.xtable(xtable(topper_15), type = "latex", file = paste0(report_path, 'univariate_covid_other.tex'))
volcano_plot = volcano_results(fitDupCor, paste0(': COVID-19 vs All Others'), 2)$vol_plot
vol_sig_genes[[1]] = volcano_results(fitDupCor, paste0(': COVID-19 vs All Others'), 2)$sig_genes
ggsave(paste0(figure_path, 'volcano_other.png'), width=10, height=6, dpi=300, units="in", device="png")
for(k in 1:4){
  top_gene_boxplot(rna_nmf_nlcpm, key_tmp, rownames(topper)[k], k)
  ggsave(paste0(figure_path, 'gene', k, '_other.png'), width=10, height=6, dpi=300, units="in", device="png")
}
```


for(m in 1:pan_num_genes){
  top_gene_boxplot(xpr_nlcpm, key_tmp, rownames(topper)[which(rownames(topper) %in% pan_viral_genes)[m]], which(rownames(topper) %in% pan_viral_genes)[m])
  ggsave(paste0(figure_path, 'processing/univariate/top_boxplots/pan_gene', m, '_other.png'), width=10, height=6, dpi=300, units="in", device="png")
}

# Boxplots of pan viral genes #
```{r}
for(v in 1:length(pan_viral_genes)){
  top_gene_boxplot(xpr_nlcpm, key, pan_viral_genes[v], 'NA')
  
}
```

for(v in 1:length(pan_viral_genes)){
  top_gene_boxplot(xpr_nlcpm, key, pan_viral_genes[v], 'NA')
  ggsave(paste0(figure_path, 'processing/univariate/pan_viral_genes/', pan_viral_genes[v], '.png'), width=10, height=6, dpi=300, units="in", device="png")
}
# Tables #
```{r}

```

table(key$pathogen, key$cohort)
print.xtable(xtable(table(key$pathogen, key$cohort)), type = "latex", file = paste0(report_path, 'table_pathogen_by_cohort.tex'))

t(table(key$cohort))
print.xtable(xtable(t(table(key$cohort))), type = "latex", file = paste0(report_path, 'table_cohort.tex'))

num_table = matrix(NA, nrow=3, ncol=(length(levels(key$cohort))+1))
colnames(num_table) = c(levels(key$cohort), "Total")
row.names(num_table) = c("Subjects", "Samples", "Total")

for(i in 1:length(levels(key$cohort))){
  num_table[1, i] = length(unique(key$subject_id[which(key$cohort == levels(key$cohort)[i])]))
  num_table[2, i] = length(key$rna_id[which(key$cohort == levels(key$cohort)[i])])
}
num_table[3, ] = colSums(num_table, na.rm=T)
num_table[, length(levels(key$cohort))+1] = rowSums(num_table, na.rm=T)
num_table[3, length(levels(key$cohort))+1] = NA
num_table
print.xtable(xtable(num_table), type = "latex", file = paste0(report_path, "num_table.tex"))

# Venn Diagram #
```{r}
venn_list <- list(Other_Coronavirus=vol_sig_genes[[1]],
                  Influenza=vol_sig_genes[[2]],
                  Bacterial=vol_sig_genes[[3]])
Vstem <- Venn(venn_list)
C <- compute.Venn(Vstem, doWeights = FALSE, type = "circles")
gp <- VennThemes(C)
```

venn_list <- list(Other_Coronavirus=vol_sig_genes[[1]],
                  Influenza=vol_sig_genes[[2]],
                  Bacterial=vol_sig_genes[[3]])
Vstem <- Venn(venn_list)
C <- compute.Venn(Vstem, doWeights = FALSE, type = "circles")
gp <- VennThemes(C)
png(paste0(figure_path, "processing/univariate/venn_diagrams/venn_corona_influenza_bacterial.png"), width=9, height=8, units="in", res=300)
plot(C, gpList = gp)
dev.off()

# Venn Diagram #
venn_list <- list(Other_Coronavirus=vol_sig_genes[[1]],
                  Influenza=vol_sig_genes[[2]],
                  Healthy=vol_sig_genes[[4]])
Vstem <- Venn(venn_list)
C <- compute.Venn(Vstem, doWeights = FALSE, type = "circles")
gp <- VennThemes(C)
png(paste0(figure_path, "processing/univariate/venn_diagrams/venn_corona_influenza_healthy.png"), width=9, height=8, units="in", res=300)
plot(C, gpList = gp)
dev.off()

# Venn Diagram #
venn_list <- list(Other_Coronavirus=vol_sig_genes[[1]],
                  Influenza=vol_sig_genes[[2]],
                  All_Others=vol_sig_genes[[5]])
Vstem <- Venn(venn_list)
C <- compute.Venn(Vstem, doWeights = FALSE, type = "circles")
gp <- VennThemes(C)
png(paste0(figure_path, "processing/univariate/venn_diagrams/venn_corona_influenza_other.png"), width=9, height=8, units="in", res=300)
plot(C, gpList = gp)
dev.off()
