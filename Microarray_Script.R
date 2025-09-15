# =============================================
# Microarray Data Analysis – R Script
# =============================================

### Step 1: Install Required Packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor packages
BiocManager::install(c("affy", "limma", "makecdfenv", "affxparser"))

# CRAN packages
install.packages(c("curl", "RCurl", "cpp11", "RSQLite", "writexl"))

### Step 2: Load Libraries
library(gdata)
library(gplots)
library(affy)
library(limma)
library(makecdfenv)
library(plogr)
library(BH)
library(affxparser)
library(writexl)

### Step 3: Read Microarray Data
data <- ReadAffy()
print(data)
View(data)

data_expr <- exprs(data)
View(data_expr)

# Initial Boxplots (Raw Data)
boxplot(data)
boxplot(data, las = 2)
boxplot(data, las = 2, col = c("red", "red", "pink", "pink", "pink", "pink"))

### Step 4: Normalize Data (RMA Normalization)
norm_data <- rma(data)
View(norm_data)

norm_data_count <- exprs(norm_data)
View(norm_data_count)

# Boxplot after normalization
boxplot(norm_data_count, las = 2, col = c("red", "red", "pink", "pink", "pink", "pink"))

### Step 5: Statistical Analysis – t-test per Gene
p_val <- apply(norm_data_count, 1, function(x) {
  t.test(x[1:2], x[3:6])$p.value
})
View(p_val)

# Combine expression data with p-values
combine_pval <- cbind(norm_data_count, p_val)
View(combine_pval)

### Step 6: Differential Gene Expression Analysis
# Calculate mean expression of control samples
mean_ctrl <- rowMeans(combine_pval[, c("Ctrl_15.CEL", "Ctrl_16.CEL")], na.rm = TRUE)

# Combine data
all_table <- cbind(combine_pval, mean_ctrl)
View(all_table)

# Calculate fold change: Log2(Treated) - Log2(Control)
fold_change <- all_table[, 3:6] - all_table[, 8]
colnames(fold_change) <- c("T1 vs C", "T2 vs C", "T3 vs C", "T4 vs C")
View(fold_change)

# Final table combining all results
fold_change_table <- cbind(all_table, fold_change)
View(fold_change_table)

### Step 7: Export Results
# Save results as Excel file
fold_change_df <- as.data.frame(fold_change_table)
writexl::write_xlsx(fold_change_df, "microarray.xlsx")

# Save as text and CSV files
write.table(fold_change_table, file = "microarray_TvsC.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(fold_change_table, file = "microarray_TvsC.csv", sep = ",", quote = FALSE, row.names = TRUE, col.names = NA)

# Example filtering
filtered_genes <- fold_change_table[abs(fold_change_table$`T1 vs C`) >= 1 & fold_change_table$p_val <= 0.05, ]
View(filtered_genes)
