#Libraries
library("dplyr")
library("tidyr")

args <- commandArgs(trailingOnly=TRUE)
qual.files <- list.files(".",
                         pattern="rnaseq_qc_results.txt", 
                         recursive=TRUE, full.names=TRUE)

#can iterate over folders what we need
##fields stored in Qualimap
qual.headings <- c("bam file", "reads aligned","total alignments","secondary alignments",
                   "non-unique alignments","aligned to genes","ambiguous alignments","no feature assigned",
                   "not aligned","SSP estimation","exonic","intronic","intergenic","overlapping exon","rRNA",
                   "5' bias","3' bias","5'-3' bias","reads at junctions")
out.mat <- matrix(nrow=length(qual.files),ncol=length(qual.headings))
dimnames(out.mat)[[2]] <- qual.headings

for(x.int in 1:length(qual.files)){
  tmp2 <- scan(qual.files[x.int], sep="\n",what="character",comment.char = "",na.strings = "",quiet=T)
  for(i in 1:length(qual.headings)){
    a <- tmp2[grep(qual.headings[i],tmp2)]
    b <- unlist(strsplit(a,"\\="))
    #these 2 lines will get rid of most of the spacial 'fill'
    b <- gsub(" ","",b)
    b <- gsub(",","",b)
    #remove auto-percentages
    if(length(grep("\\(",b[2]))>0){ #don't die if no "(" exists
      d <- unlist(strsplit(b[2],"\\("))
      out.mat[x.int,i] <- d[1]
    }
    else{
      out.mat[x.int,i] <- b[2]
    }
  }
}

##Table generating all metrics no QC
write.csv(out.mat,"qualimap_stats_mqc.csv", row.names = FALSE, quote = FALSE)

##############################
## Generate summary metrics on qualimap
##############################
df<-read.csv("qualimap_stats_mqc.csv", sep=",", header = TRUE)

#df$experiment_name<-gsub("_L001.*", "", df$bam.file)
#df$bulk_sc<-ifelse(grepl("100pg|NTC", df$bam.file, ignore.case = T), "bulk", "sc")
#df$group<-ifelse(grepl("OldW", df$`bam file`, ignore.case = T), "Old Wash", "NewWash")
#df$molecule<-ifelse(grepl("ROR", df$bam.file, ignore.case = T), "RNA", "DNA")

##summarize percent of total 
df<-df %>%
  mutate_at(vars(c("reads.aligned", "exonic", "intronic", "intergenic", "overlapping.exon")) , funs(P_Total = ./df$total.alignments * 100))

##summarize percent of all genomic features
df<-df %>%
  mutate(Total = dplyr::select(., exonic:overlapping.exon) %>% rowSums(na.rm = TRUE))

df<-df%>%
  mutate_at(vars(exonic:overlapping.exon) , funs(P_gen = ./df$Total * 100))


df$reads.aligned_P_Total = df$reads.aligned / (df$reads.aligned + df$not.aligned)
df$reads.aligned_P_Total <- df$reads.aligned_P_Total * 100


#Table generating summarized metrics with benchmark QC

df<-df %>%
  select(bam.file, total.alignments, reads.aligned_P_Total, exonic_P_gen, intergenic_P_gen)

write.csv(df, "pipeline_metrics_summary_percents.csv", row.names = FALSE, quote = FALSE)





