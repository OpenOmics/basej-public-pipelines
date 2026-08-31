library(ohchibi)
library(pals)
library(dplyr)
library(ggrepel)
library(ggtree)
library(ggpmisc)


### Functions 
source('/usr/local/bin/rnaseq_custom_functions.R')

################ Read the star structure ####
file <- "df_gene_counts_starhtseq.tsv"

df <- read.table(file,sep="\t",header=T,comment.char="",quote = "") %>%
  dplyr::mutate(.data =.,File = File  %>% gsub(pattern = "_L00.*.htseq_counts.tsv",replacement = ""))

df$NumberName <- df$File %>% factor %>% as.numeric  %>%
  paste0("S",.) %>% factor

df_trans <- df[,c("NumberName","File")] %>% unique
colnames(df_trans)[1] <- "Sample"
max_length <- max(nchar(as.character(df_trans$File)))
p_names <- ggplot() + annotate(geom = "table",x = 1,y = 1,label = list(df_trans)) + theme_void()


#Create matrix
Tab <- acast(data = df,formula=gene_symbol_gene_id~NumberName,value.var = "countHTSeq",fill = 0)  %>% t   

Map <- data.frame(SampleId = rownames(Tab)) 
rownames(Map) <- Map$SampleId

Tax <- data.frame(ID = colnames(Tab),
                  gene_id = colnames(Tab)%>% gsub(pattern = ".*_ENSG",replacement = "ENSG")) %>%
      dplyr::mutate(.data = .,primerid = ID)
rownames(Tax) <- Tax$ID



Map <- Tab %>% apply(X = .,MARGIN = 1,function(x)which(x >0) %>% length) %>%
  data.frame(SampleId = names(.),n_genes_per_cell = .) %>%
  merge(Map, ., by = "SampleId")
rownames(Map) <- Map$SampleId

Dat <- create_dataset(Tab = Tab %>% t,Map = Map,Tax = Tax)


#Apply normalization
Res_wrap <- wrapper_filter_norm(min_cells = 5,Dat = Dat,top_genes = 500,normalize = TRUE,num_genes_per_cell = 1000)

#Perform PCA projection using top 500 hvg 

mpca <- oh.pca(Tab = Res_wrap$Dat_hvg$Tab  %>% t,
       Map = Res_wrap$Dat_hvg$Map,retx = T,center = T,scale = F,id_var = "SampleId")
p_pca <- chibi.pca(list_ohpca = mpca,size =  4,stroke = 0.3) +
  geom_text_repel(aes(label = SampleId),max.overlaps = 1000)+
  theme_ohchibi(size_panel_border = 0.3) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ggtitle("PC: Top 500 genes")

#Perform heatmap using the top 500 hvg genes
Tab <- Res_wrap$Dat_hvg$Tab %>% t %>% scale(center = TRUE,scale = FALSE)  %>% t

res_heatmap <- chibi.heatmap(Tab = Tab,k_rows = 1,k_cols = 1,dist_method_rows = "euclidean",
                             dist_method_cols = "euclidean",hclust_method_rows = "ward.D",
                             hclust_method_cols = "ward.D",range_fill_heatmap = c(-2,2),axis_ticks_row = TRUE,
                             size_axis_text_row  = 4,
                             mtheme = theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust  = 1) ),
                             egg_heights = c(0.05, 1),
                             egg_widths = c(0.05, 1)
                             )

p_heatmap <- res_heatmap$heatmap

# Plot the information  about the cells
df_mt <- read.table(file = "df_mt_gene_counts_starhtseq.tsv",header = T,sep = "\t") %>%
  dplyr::mutate(.data = .,SampleId = File %>% gsub(pattern = "_L00.*",replacement = ""))

p_mt <- ggplot(data = df_mt,aes(SampleId,PropMT)) +
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = 0.1,color = "red",size = 0.5)+
  theme_ohchibi(size_panel_border = 0.3) +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1)) +
  theme(
    axis.text.x = element_text(size = 8,angle = 90,vjust = 0.5,hjust = 1)
  ) +
  ylab(label = "Proportion of counts mapping to mitochondrial genes") +
  xlab(label = "Sample id")

#Plot number of protein genes 
df_gt <- read.table(file = "df_gene_types_detected_summary_starhtseq.tsv",header = T,sep = "\t") %>%
  dplyr::mutate(.data = .,SampleId = File %>% gsub(pattern = "_L00.*",replacement = "")) %>%
  subset(gene_biotype == "protein_coding") %>%
  droplevels
  
p_gt <- ggplot(data = df_gt,aes(SampleId,NumFeatures)) +
  geom_bar(stat = "identity") + 
  theme_ohchibi(size_panel_border = 0.3) +
  scale_y_continuous(breaks = seq(0,10000,500),limits = c(0,10000),oob = rescale_none) +
  theme(
    axis.text.x = element_text(size = 8,angle = 90,vjust = 0.5,hjust = 1)
  ) +
  ylab(label = "Number of protein coding features detected (Count > 0)") +
  xlab(label = "Sample id")



ggsave("pca_mqc.png",plot = p_pca,width = 9 ,height = 9, units = "in")

ggsave("heatmap_mqc.png",plot = p_heatmap,width = 16 ,height = 28, units = "in")
ggsave("mt_mqc.png",plot = p_mt,width = 12 ,height = 11, units = "in")
ggsave("genesdetected_mqc.png",plot = p_gt,width = 12 ,height = 11, units = "in")

#Establish a dynamic size based on the number of samples
minches <- ceiling(nrow(df_trans)/5)
#Adjusting the width based on the length of the filename
width_adjusted <- ifelse(max_length < 50, 5, ifelse(max_length <= 80, 8, 12))
ggsave("samplenames_mqc.png",plot = p_names,width = width_adjusted ,height = minches, units = "in")

