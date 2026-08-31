library(magrittr)
library(dplyr)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
file_name <- args[1]
tx2gene2type_file <- args[2]

tx2gene2type <- read.table(file = tx2gene2type_file,header = F,quote = "",comment.char = "") %>%
	dplyr::rename(.data = .,transcript_id = V1,transcript_biotype = V2,gene_id_raw = V3,gene_biotype=V4,gene_symbol = V5) %>%
	dplyr::mutate(.data =.,gene_id = gene_id_raw %>% gsub(pattern = "\\..*",replacement = "")) %>%
	dplyr::select(.data =.,c(gene_id,gene_id_raw,gene_biotype,gene_symbol)) %>% unique


file <- file_name

df_file <- read.table(file = file,header = F,sep = "\t")  %>%
    dplyr::rename(.data =.,gene_id = V1,gene_symbol = V2,countHTSeq = V3) %>%
    dplyr::select(.data =.,-c("gene_symbol")) %>% 
    dplyr::mutate(.data=.,File = file) %>%
    dplyr::filter(.data=.,! gene_id %>% grepl(pattern = "^__")) %>% 
    dplyr::mutate(.data=.,gene_id = gene_id %>% gsub(pattern = "\\..*",replacement = "")) %>%
    merge(.,tx2gene2type,by = "gene_id",all.x = TRUE) %>%
    dplyr::relocate(.data =.,c(File,gene_id_raw,countHTSeq,gene_biotype,gene_symbol)) %>%
    dplyr::select(.data = .,-c("gene_id")) %>%
    dplyr::rename(.data =.,gene_id = gene_id_raw) %>% 
    subset(countHTSeq != 0) %>% droplevels %>%
    dplyr::mutate(.data =.,gene_symbol_gene_id = paste0(gene_symbol,"_",gene_id))

#Here determine the size is bigger than 0 if not we create empty data frames
if(nrow(df_file) != 0){

#Create summary quantification
df_num_gene_biotype <- df_file %>% 
	dplyr::mutate(.data = .,Dummy = 1) %>%
	aggregate(Dummy~gene_biotype,.,sum) %>%
	dplyr::mutate(.data = ., PropFeatures = Dummy/sum(Dummy),File = file) %>%
	dplyr::rename(.data = .,NumFeatures = Dummy) %>%
	relocate(File)
	
	
df_count_gene_biotype <- df_file %>% 
	aggregate(countHTSeq~gene_biotype,.,sum) %>%
	dplyr::mutate(.data = ., PropcountHTSeq = countHTSeq/sum(countHTSeq),File = file) %>%
	relocate(File)
	
### Merge structure
df_num_gene_biotype <- merge(df_num_gene_biotype,df_count_gene_biotype, by = c("File","gene_biotype"),all = TRUE)

#Add Mitochondrial quantification based on MT tag of gene
num_mt_features_detected <- df_file$countHTSeq[df_file$gene_symbol %>% grep(pattern="^MT-")]  %>% length
sum_mit <- df_file$countHTSeq[df_file$gene_symbol %>% grep(pattern="^MT-")] %>% sum                                                                                                                                                                                                                 
prop_mit <- sum_mit/(df_file$countHTSeq %>% sum)    

df_mt <- data.frame(File = file,TotalFeatures = df_file$gene_id %>% as.character %>% unique %>% length,
    MT_NumFeatures = num_mt_features_detected,
    MT_Counts = sum_mit,Total_Counts = df_file$countHTSeq %>% sum,PropMT = prop_mit )
    
}else{

    df_mt <- data.frame(File = file,TotalFeatures = 0,
    MT_NumFeatures = 0,
    MT_Counts = 0,Total_Counts = 0,PropMT = 0 )
    
    df_num_gene_biotype <- data.frame( File = file,
    gene_biotype = "protein_coding",
    NumFeatures = 0,
    PropFeatures = 0,
    countHTSeq = 0,
    PropcountHTSeq = 0
    )

    
}



write.table(df_file,file = paste0("df_gene_star_htseq_",file_name),
	append = F,sep = "\t",col.names = T,row.names = F,quote = F)


write.table(df_num_gene_biotype,file = paste0("df_sum_detected_gene_star_htseq_",file_name),
	append = F,sep = "\t",col.names = T,row.names = F,quote = F)


write.table(df_mt,file = paste0("df_mtcounts_star_htseq_",file_name),
	append = F,sep = "\t",col.names = T,row.names = F,quote = F)
	
