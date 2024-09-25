#########################################################
# Gene expression profiling of ovarian cancer 
#- TCGA-OV project - RNAseq data  

# Loading packages --------------------------------------------------------
library(readr)
library(org.Hs.eg.db)
library(DESeq2)

# Reading Data ------------------------------------------------------------
data.path="C:/Users/marie/Downloads/ovary_data/ovarian_data_TCGA"
pheno.path="C:/Users/marie/Downloads/ovary_data/ovarian_metadata.tsv"

#loading RNA-seq data
files = list.files(data.path, recursive = T, pattern = "tsv")
length(files)

# read the first file for the first time
file = files[1]
filepath = file.path(data.path, files[1])
colnames = unlist(strsplit(readLines(filepath, n=2)[2], "\t"))

#reading the table:
temp = read.table(filepath, header = F, skip = 6)
names(temp) = colnames

genetype_mapper = temp[, c(2,3)]

#creating a storing object in which we store the TPMs of each file
exp = temp[temp$gene_type == "protein_coding", c(1,2)]

for (i in 1:length(files)) {
  file = files[i]
  file.id = strsplit(file, "/")[[1]][1]
  filepath = file.path(data.path, files[i])
  
  temp = read.table(filepath, header = F, skip=6)
  temp = temp[temp[, 3]=="protein_coding", ]
  
  exp = cbind(exp, temp[, 4])
  colnames(exp)[dim(exp)[2]] = file.id
}


# Quality Control ---------------------------------------------------------
sum(duplicated(exp$gene_name)) #there are 24 duplicates

#so let's remove these duplicates by aggregation
exp.data = exp[, 3:dim(exp)[2]]
exp.data = apply(exp.data, 2, as.numeric)
dim(exp.data)

exp.data.agg = aggregate(exp.data, list(exp$gene_name), FUN = mean)
dim(exp.data.agg)
genes =exp.data.agg$Group.1
exp.data.agg = exp.data.agg[, -1]
exp.data.agg = apply(exp.data.agg, 2, as.numeric)
exp.data.agg = round(exp.data.agg)
rownames(exp.data.agg) = genes



# Loading the sample sheet ------------------------------------------------
pheno <- read_delim(pheno.path,"\t", escape_double = FALSE, trim_ws = TRUE)
table(pheno$`Sample Type`) #there are 421 primar tumor samples and 8 recurrent tumor samples

#adjusting colnames of phenotable:
names(pheno) = sub(" ", "_", names(pheno))

########changing colnames of exp.data.agg to sample names from the phenotable:
file.ids.pheno = pheno$File_ID
file.ids.data = colnames(exp.data.agg)

index.files = match(file.ids.data, file.ids.pheno)
colnames(exp.data.agg) = pheno$Sample_ID[index.files]


dim(exp.data.agg)
# Exploratory Analysis ----------------------------------------------------

#checking shape of my data: features must be the rownames and samples are the colnames


#checking mode of data:
mode(exp.data.agg) #it is numeric


# Filtration -------------------------------------------------------------
#removing constant genes:
dim(exp.data.agg)   #19938 rows and 429 columns
rowvar = apply(exp.data.agg, 1, var, na.rm = T)
constvar = (rowvar == 0 | is.na(rowvar))
sum(constvar) #there are 480 constant genes
exp.data.agg = exp.data.agg[!constvar, ]
dim(exp.data.agg)  #19520 rows

#Handling NAs:
sum(apply(exp.data.agg, 1, is.na)) #zero nulls

mat_clean = exp.data.agg
#save(mat_clean, file = "C:/Users/marie/Downloads/ovary_data/mat_clean.RDATA")
load("C:/Users/marie/Downloads/ovary_data/mat_clean.RDATA")

#Removing the least 25% vaariant genes 
0.25*19962  #number of genes thaat should be removed (4990.5)
sds = rowSds(mat_clean)
mat_clean = as.data.frame(cbind(mat_clean, sds))
mat_clean = mat_clean[order(mat_clean$sds), ]
mat_clean = mat_clean[-c(1:4990), ]
dim(mat_clean)

exp = mat_clean[, -430]
write.csv(exp, file = "C:/Users/marie/Downloads/ovary_data/exp.csv")




# Differential Expression Analysis ----------------------------------------
#install.packages(c("ggplot2", "dplyr", "pheatmap", "RColorBrewer", "ggrepel"))

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

#setwd("D:\\bioinformatics\\hackbio")

sample_info = pheno
count_data = exp

dds2 <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = ~Sample_Type)

dds2$Sample_Type <- factor(dds$Sample_Type, levels = c("Primary Tumor", "Recurrent Tumor"))
dds2 <- DESeq(dds2)
deseq_results <- results(dds2)

deseq_results <- as.data.frame(deseq_results)
class(deseq_results)

filtered <- deseq_results %>% filter(deseq_results$padj < 0.05)
filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1)

upregulated <- deseq_results[deseq_results$log2FoldChange > 1 & deseq_results$padj < 0.05, ]
downregulated <- deseq_results[deseq_results$log2FoldChange < -1 & deseq_results$padj < 0.05, ]

write.csv(downregulated, file = "downregulated_genes.csv", row.names = TRUE)
# Save to CSV files
write.csv(upregulated, file = "upregulated_genes.csv", row.names = TRUE)

# Extract normalized count data from DESeq2 object
ntd <- normTransform(dds2)
normalized_counts = assay(ntd)
# Subset normalized counts for upregulated genes
upregulated_genes <- rownames(upregulated)
upregulated_counts <- normalized_counts[upregulated_genes, ]

# Subset normalized counts for downregulated genes
downregulated_genes <- rownames(downregulated)
downregulated_counts <- normalized_counts[downregulated_genes, ]


# Create heatmap for upregulated genes
pheatmap(upregulated_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         main = "Heatmap of Upregulated Genes")


pheatmap(downregulated_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = TRUE, 
         show_colnames = FALSE, 
         fontsize_row = 4,
         main = "Heatmap of Downregulated Genes")



# Load the library
library(EnhancedVolcano)

# Volcano plot using your deseq_results
EnhancedVolcano(deseq_results,
                lab = rownames(deseq_results),
                x = 'log2FoldChange',  # x-axis for fold change
                y = 'pvalue',          # y-axis for p-values
                pCutoff = 0.05,        # Adjust as needed
                FCcutoff = 1,          # Fold-change threshold
                title = 'Volcano Plot',
                xlab = 'Log2 Fold Change',
                ylab = '-Log10 p-value',
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendLabels = c('NS', 'Log2 FC', 'p-value', 'p-value & Log2 FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                labSize = 3.0)


# preparing data for machine learning -------------------------------------
#merging pheno with t.degs.exp:

res.degs = deseq_results[deseq_results$padj<0.05 & abs(deseq_results$log2FoldChange) > log2(2), ]
degs = rownames(res.degs)
degs.exp = normalized_counts[degs, ]
t.degs.exp = t(degs.exp)
t.degs.exp <- read.csv("C:/Users/marie/Downloads/ovary_data/t.degs.exp.csv")
mode(t.degs.exp)
colnames(t.degs.exp)[colnames(t.degs.exp) == "X"] = "Sample_ID"
t.degs.exp$Sample_ID = gsub("\\.", "-", t.degs.exp$Sample_ID)

merged.exp = merge(t.degs.exp, pheno, by = "Sample_ID", all.x = T)


#merged.exp = merged.exp[, -c(139:144)]
rownames(merged.exp) = merged.exp$Sample_ID
#merged.exp = merged.exp[, -1]
primary_tumor = merged.exp[merged.exp$Sample_Type == "Primary Tumor", ]
recurrent_tumor = merged.exp[merged.exp$Sample_Type == "Recurrent Tumor", ]

selected_primary = primary_tumor[sample(1:nrow(primary_tumor), 8), ]

min_matrix = rbind(selected_primary, recurrent_tumor)
dim(min_matrix)
write.csv(min_matrix, file = "C:/Users/marie/Downloads/ovary_data/filter_matrix.csv")
