library(reshape2)
library(magrittr)
library(dplyr)
library(tibble)
library(pheatmap)


#Process  Star matrix
df_star <- read.table(file = "./df_gene_counts_starhtseq.tsv",header = TRUE,sep = "\t") %>%
  dplyr::mutate(.data = .,File = File %>% gsub(pattern = ".htseq_counts.tsv",replacement = "") )

Tab_STAR_Raw <- acast(data = df_star,formula = File~gene_id+gene_symbol,
                            value.var = "countHTSeq",fill = 0) %>% t %>%
                            as.data.frame %>%
                            rownames_to_column %>%
                            dplyr::rename(.data = .,gene_id_gene_symbol = rowname)


write.table(x = Tab_STAR_Raw,file = "./matrix_gene_counts_starhtseq.txt",append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

# Housekeeping genes analysis
hk.genes <- c("GAPDH","ACTB","RPL36","HINT1","TBP","PPIA","B2M","HPRT1","UBC","RPL13A",
              "PGK1","EIF3K","RPLP0","GUSB","CLTC","HMBS")

# Use Tab_STAR_Raw directly and set row names
x.ge <- Tab_STAR_Raw %>% column_to_rownames(var = "gene_id_gene_symbol")

# Extract gene symbols from the gene_id_gene_symbol column (row names)
gene.symbols <- sapply(strsplit(rownames(x.ge), "_"), function(x) x[length(x)])

# Filter for housekeeping genes
hk.idx <- gene.symbols %in% hk.genes
x.hkg <- x.ge[hk.idx, , drop = FALSE]

# Write housekeeping gene counts per sample
x.hkg.counts <- as.data.frame(x.hkg) %>% 
  rownames_to_column(var = "Gene")

write.table(x.hkg.counts, "HouseKeepingGenes_Counts_mqc.tsv", quote=F, na="", sep="\t", row.names=FALSE, col.names=TRUE)

# Check if we have multiple samples for CV and clustering analysis
n_samples <- ncol(x.hkg)

if (n_samples > 1) {
  # Create CV rates for genes (Coefficient of Variation)
  x.hkg.CV <- data.frame(
    Gene = rownames(x.hkg),
    CV = as.numeric(apply(x.hkg, 1, sd, na.rm=T) / apply(x.hkg, 1, mean, na.rm=T)),
    stringsAsFactors = FALSE
  )
  
  # Write CV table
  write.table(x.hkg.CV, "HouseKeepingGenes_CV_mqc.tsv", quote=F, na="", sep="\t", row.names=FALSE, col.names=TRUE)
  
  # Create clustergram for housekeeping genes
  # Generate PDF
  pheatmap(x.hkg, 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           scale = "row",
           show_rownames = TRUE, 
           show_colnames = TRUE,
           main = "HouseKeeping Genes Expression Clustergram",
           fontsize_row = 10,
           fontsize_col = 10,
           filename = "HouseKeepingGenes_Expression.pdf",
           width = 10,
           height = 8)
  
  # Generate PNG
  pheatmap(x.hkg, 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           scale = "row",
           show_rownames = TRUE, 
           show_colnames = TRUE,
           main = "Housekeeping Genes Expression Clustergram",
           fontsize_row = 10,
           fontsize_col = 10,
           filename = "HKGenes_Expression__mqc.png",
           width = 10,
           height = 8)
} else {
  # Single sample - create empty/placeholder files
  cat("# CV analysis requires multiple samples\n", file = "HouseKeepingGenes_CV_mqc.tsv")
  cat("# Clustergram requires multiple samples\n", file = "HKGenes_Expression__mqc.png")
  cat("# Clustergram requires multiple samples\n", file = "HouseKeepingGenes_Expression.pdf")
}


