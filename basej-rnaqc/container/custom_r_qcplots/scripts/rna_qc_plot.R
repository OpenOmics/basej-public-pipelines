library(ohchibi)
library(dplyr)
library(ggtree)
library(optparse)

set.seed(130816)

bskb_col<-c("#12284C", "#1082A2","#A0CC2C" , "#DD14D3","#F45D34","#777776", "firebrick","orange","darkgreen","pink")


filter_expresion_matrix <- function(expr_mat, min_num_cells_per_gene, sample_names)  {
  filtered_expr_mat <-expr_mat[sample_names,]
  logical_vec <- apply(filtered_expr_mat, 2, function(.col) {
    return(sum(as.logical(.col)) >= min_num_cells_per_gene)
  })
  filtered_expr_mat <- filtered_expr_mat[, logical_vec] 
  return(filtered_expr_mat)
}

normalize_matrix <- function(expr_mat) {
  #Rows should be samples 
  expr_normed_mat <- t(apply(expr_mat, 1, function(.row) {log((.row/sum(.row)) * 10000  + 1)}))
  return(expr_normed_mat)
}

get_dispersion <- function(log_normed_mat, topn = 500)  {
  temp_mean_and_var <- apply(log_normed_mat, 2, function(.col)  {
    return(c(mean(.col), var(.col)))
  })
  temp_disp_mean_df <- data.frame(ensembl = colnames(temp_mean_and_var), mean_normed_gene_expr = temp_mean_and_var[1,], var_normed_gene_exp = temp_mean_and_var[2,], exp_bin = cut(temp_mean_and_var[1,], 20), stringsAsFactors = FALSE) %>% dplyr::mutate(gene_normed_exp_disp = var_normed_gene_exp/mean_normed_gene_expr)
  temp_disp_mean_df <- dplyr::left_join(temp_disp_mean_df, {dplyr::group_by(temp_disp_mean_df, exp_bin) %>% dplyr::summarize(mean_bin_expression = mean(mean_normed_gene_expr), mean_bin_dispersion = mean(gene_normed_exp_disp), sd_bin_dispersion = sd(gene_normed_exp_disp))}) %>% dplyr:: mutate(abs_normalized_bin_dispersion_deviation = abs((gene_normed_exp_disp - mean_bin_dispersion) / sd_bin_dispersion))
  temp_disp_mean_df <- dplyr::arrange(temp_disp_mean_df, desc(abs_normalized_bin_dispersion_deviation))
  topn_dispersed_mat <- log_normed_mat[,temp_disp_mean_df$ensembl[1:topn]]
  return(topn_dispersed_mat)
}


wrapper_filter_norm <- function(Dat = NA,min_cells = 1,num_genes_per_cell = 1000,top_genes = 500,normalize = TRUE){
  Tab <- Dat$Tab %>% t
  cell_stats_df <- Dat$Map
  Tax <- Dat$Tax
  min_cells <- min_cells
  subset_sample_names <- dplyr::filter(cell_stats_df,  n_genes_per_cell > num_genes_per_cell) %>% dplyr::pull(SampleId)
  subset_cell_stats_df <- cell_stats_df %>% dplyr::filter(SampleId %in% subset_sample_names)
  cell_and_gene_filtered_mat <-  filter_expresion_matrix(Tab, min_cells, subset_sample_names)
  if(normalize == TRUE){
    cell_and_gene_filtered_normed_mat <- normalize_matrix(cell_and_gene_filtered_mat)
    cell_and_gene_filtered_normed_hvg_mat <- get_dispersion(cell_and_gene_filtered_normed_mat, top = top_genes)
    
  }else{
    cell_and_gene_filtered_normed_mat <- cell_and_gene_filtered_mat
    cell_and_gene_filtered_normed_hvg_mat <- get_dispersion(cell_and_gene_filtered_mat, top = top_genes)
  }
  
  Tax_sub <- match(colnames(cell_and_gene_filtered_normed_mat),Tax$ID) %>%
    Tax[.,]
  
  #Create dat objects to return 
  Dat_norm_temp <- create_dataset(Tab =cell_and_gene_filtered_normed_mat %>% t,Map = subset_cell_stats_df,Tax_sub )
  
  Tax_sub <- match(colnames(cell_and_gene_filtered_normed_hvg_mat),Tax$ID) %>%
    Tax[.,]
  Dat_norm_hvg_temp <- create_dataset(Tab =cell_and_gene_filtered_normed_hvg_mat %>% t,Map = subset_cell_stats_df )
  
  return(list(Dat = Dat_norm_temp,
              Dat_hvg = Dat_norm_hvg_temp))
  
}

plot_qc_rna <- function(matrix_file = NULL,
                        metrics_file = NULL,
                        metadata_file = NULL,
                        min_cells = 5,
                        top_genes = 500,
                        cutoff_PropMappability = 0.7,
                        cutoff_PropExonic = 0.7,
                        cutoff_PropIntergenic = 0.1,
                        cutoff_ProportionCountsMitochondrialGenes = 0.1,
                        cutoff_ProteinCodingGenesDetected = 500,
                        size_labels_samples = 1,
                        relative_size_heatmap = 1,
                        relative_size_mappability = 0.2,
                        relative_size_exon = 0.2,
                        relative_size_intergenic = 0.2,
                        relative_size_chrM = 0.2,
                        relative_size_numprotgenes = 0.2,
                        relative_size_dendrogram = 0.2,
                        relative_size_bar_plot_group = 0.03
){
  
  #Read tables and determine usable samples
  Tab <- read.table(file = matrix_file,header = TRUE,row.names = 1,check.names = FALSE) %>% t
  
  
  #Read information about 
  df <- read.table(file = metrics_file,
                   header = TRUE,sep = ",",check.names = TRUE,row.names = 1)  %>%
    tibble::rownames_to_column() %>%
    dplyr::rename(.data =.,SampleId = rowname)
  
  
  ids_tab <- rownames(Tab)
  ids_mets <- df$SampleId
  
  #Check if metadata file was passed
  if(is.null(metadata_file)){

    Map <- data.frame(SampleId = rownames(Tab))
    rownames(Map) <- Map$SampleId
    # Add default group column even when metadata is not provided
    Map$group <- "Group1"

    usable_ids <- intersect(ids_tab,ids_mets)

  }else{
    
    Map <- read.table(file = metadata_file, header = TRUE, sep = ",", check.names = FALSE)
    # Find index for sample column and rename to SampleId
    sample_col_idx <- which(colnames(Map) %in% c("SampleId", "biosampleName"))[1]
    colnames(Map)[sample_col_idx] <- "SampleId"

    # Handle group/groups column
    if ("group" %in% colnames(Map) && "groups" %in% colnames(Map)) {
      # Both exist - remove 'groups' and keep 'group'
      Map <- Map[, colnames(Map) != "groups", drop = FALSE]
    } else if ("groups" %in% colnames(Map) && !("group" %in% colnames(Map))) {
      # Only 'groups' exists - rename it to 'group'
      colnames(Map)[colnames(Map) == "groups"] <- "group"
    } else if (!("group" %in% colnames(Map)) && !("groups" %in% colnames(Map))) {
      # Neither exists - create default group for all samples
      Map$group <- "Group1"
    }
    # If only 'group' exists, no action needed

    # Set rownames as SampleId
    rownames(Map) <- Map$SampleId
    
    ids_meta <- Map$SampleId
    usable_ids <- intersect(ids_tab, ids_mets) %>% intersect(ids_meta)
    
    Map <- match(usable_ids,rownames(Map)) %>% Map[.,]
    
  }

  # Error if less than 3 samples
  if (length(usable_ids) < 3) {
      stop(paste0("Error: Only ", length(usable_ids), " sample(s) found after intersection. At least 3 samples are required for QC plot. Please check your input files."))
  }
  
  #Match usable is to Tab df_mets
  Tab <- match(usable_ids,rownames(Tab)) %>% Tab[.,]
  
  df <- match(usable_ids,df$SampleId) %>% df[.,]
  
  
  Tax <- data.frame(ID = colnames(Tab),
                    gene_id = colnames(Tab)%>% gsub(pattern = ".*_ENSG",replacement = "ENSG")) %>%
    dplyr::mutate(.data = .,primerid = ID)
  rownames(Tax) <- Tax$ID
  
  
  
  Map <- Tab %>% apply(X = .,MARGIN = 1,function(x)which(x >0) %>% length) %>%
    data.frame(SampleId = names(.),n_genes_per_cell = .) %>%
    merge(Map, ., by = "SampleId")
  
  Map <- merge(Map,df,by = "SampleId",all.x = TRUE)
  rownames(Map) <- Map$SampleId
  
  Map <- match(rownames(Tab),rownames(Map)) %>% Map[.,]
  
  Dat <- create_dataset(Tab = Tab %>% t,Map = Map,Tax = Tax)
  
  #Project matrix 
  Res_wrap <- wrapper_filter_norm(min_cells = min_cells,
                                  Dat = Dat,
                                  top_genes = top_genes ,normalize = TRUE,num_genes_per_cell = 0)
  
  
  ### Construct the heatmap using the top 500 genes
  Tab_c <- Res_wrap$Dat_hvg$Tab  %>% t %>% scale(center = TRUE,scale = FALSE) 
  
  clust_samples <- hclust(d = dist(x = Tab_c,method = "euclidean"),method = "ward.D")
  clust_genes <- hclust(d = dist(x = Tab_c %>% t,method = "euclidean"),method = "ward.D")
  
  
  #Construct adjacent panels based on populations
  melted <- Tab_c %>% reshape2::melt() %>% 
    dplyr::rename(.data =.,SampleId = Var1,Gene = Var2)
  
  melted$Gene <- melted$Gene %>% factor(levels = clust_genes$order %>% clust_genes$labels[.])
  melted$SampleId <- melted$SampleId %>% factor(levels = clust_samples$order %>% clust_samples$labels[.])
  
  Map$SampleId <- Map$SampleId %>% 
    factor(levels = clust_samples$order %>% clust_samples$labels[.])
  
  #### Construct the clustered version of the figures with composite score logic
  Mat <- matrix(data = 0,nrow = nrow(Map),ncol = 5)
  rownames(Mat) <- Map$SampleId
  
  Mat[,1][which(Map$PropMappability > cutoff_PropMappability)] <- 1
  Mat[,2][which(Map$PropExonic > cutoff_PropExonic)] <- 1
  Mat[,3][which(Map$PropIntergenic < cutoff_PropIntergenic)] <- 1
  Mat[,4][which(Map$ProportionCountsMitochondrialGenes < cutoff_ProportionCountsMitochondrialGenes)] <- 1
  Mat[,5][which(Map$ProteinCodingGenesDetected > cutoff_ProteinCodingGenesDetected)] <- 1
  
  #Create sum and that will become the cluster
  df_clust <- rowSums(Mat) %>% data.frame %>%
    tibble::rownames_to_column() %>% 
    dplyr::rename(.data =.,SampleId = rowname, Cluster = ".")
  
  df_clust$Cluster <- df_clust$Cluster %>% factor(levels = c("5","4","3","2","1","0"))
  
  #Create the representation with  quality 
  melted$Cluster <- match(melted$SampleId,df_clust$SampleId) %>% df_clust$Cluster[.]
  melted$Cluster <- melted$Cluster %>% factor(levels = c("5","4","3","2","1","0"))
  
  
  p2 <- ggplot(data = melted,aes(Gene,SampleId)) +
    geom_raster(aes(fill = value)) +
    facet_grid(Cluster~.,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3) +
    scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                           oob = squish,limits = c(-0.5,0.5),name = "Standardized\nexpression") +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = size_labels_samples),
      axis.ticks.y = element_line(size = 0.1),
      strip.text.x = element_text(size = 9,color = "grey30"),
      strip.text.y = element_blank(),
      legend.text = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 9),
      panel.spacing.y = unit(0.1, "lines"),
      
    ) +
    xlab(label = "Genes")  +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0))
  
  
  #### Append other figures
  Map$Cluster <- match(Map$SampleId,df_clust$SampleId) %>% df_clust$Cluster[.]
  
  
  df_ag <- Map[,c("SampleId","Cluster","PropMappability")] %>% unique %>%
    aggregate(PropMappability~Cluster,.,quantile)
  
  colnames(df_ag$PropMappability) <- colnames(df_ag$PropMappability) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercPropMappability")
  
  Map <- match(Map$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(Map,.)
  
  
  #
  p_map <- ggplot(data = Map,aes(PropMappability,SampleId)) +
    geom_rect(aes(xmin =PercPropMappability25 ,xmax =PercPropMappability75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercPropMappability50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_PropMappability,color = "red",linetype = "longdash")+
    #geom_bar(stat = "identity",,width = 0,fill = "grey30",color = NA) +
    geom_line(group = 1)+
    facet_grid(Cluster~.,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3)+
    theme(
      legend.position = "none",
      panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(size = unit(0.1,"line")),
      axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Proportion\nmappability") +
    scale_x_continuous(expand = c(0,0),limits = c(0,1),oob = squish,breaks = seq(0,1,0.1))
  
  ####
  df_ag <- Map[,c("SampleId","Cluster","PropExonic")] %>% unique %>%
    aggregate(PropExonic~Cluster,.,quantile)
  
  colnames(df_ag$PropExonic) <- colnames(df_ag$PropExonic) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercPropExonic")
  
  Map <- match(Map$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(Map,.)
  
  
  p_exon <- ggplot(data = Map,aes(PropExonic,SampleId)) +
    geom_rect(aes(xmin =PercPropExonic25 ,xmax =PercPropExonic75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercPropExonic50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_PropExonic,color = "red",linetype = "longdash")+
    geom_line(group = 1)+
    #geom_bar(stat = "identity",,width = 1,fill = "grey30",color = NA) +
    facet_grid(Cluster~.,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3)+
    theme(
      legend.position = "top",
      panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(size = unit(0.1,"line")),
      axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Proportion\nexonic") +
    scale_x_continuous(expand = c(0,0),limits = c(0,1),oob = squish,breaks = seq(0,1,0.1))
  
  ####
  df_ag <- Map[,c("SampleId","Cluster","PropIntergenic")] %>% unique %>%
    aggregate(PropIntergenic~Cluster,.,quantile)
  
  colnames(df_ag$PropIntergenic) <- colnames(df_ag$PropIntergenic) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercPropIntergenic")
  
  Map <- match(Map$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(Map,.)
  
  
  p_inter <- ggplot(data = Map,aes(PropIntergenic,SampleId)) +
    geom_rect(aes(xmin =PercPropIntergenic25 ,xmax =PercPropIntergenic75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercPropIntergenic50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_PropIntergenic,color = "red",linetype = "longdash")+
    geom_line(group = 1)+
    #geom_bar(stat = "identity",,width = 1,fill = "grey30",color = NA) +
    facet_grid(Cluster~.,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3)+
    theme(
      legend.position = "top",
      panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(size = unit(0.1,"line")),
      axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Proportion\nintergenic") +
    scale_x_continuous(expand = c(0,0),limits = c(0,1),oob = squish,breaks = seq(0,1,0.1))
  
  ####
  df_ag <- Map[,c("SampleId","Cluster","ProportionCountsMitochondrialGenes")] %>% unique %>%
    aggregate(ProportionCountsMitochondrialGenes~Cluster,.,quantile)
  
  colnames(df_ag$ProportionCountsMitochondrialGenes) <- colnames(df_ag$ProportionCountsMitochondrialGenes) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercProportionCountsMitochondrialGenes")
  
  Map <- match(Map$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(Map,.)
  
  
  p_mit <- ggplot(data = Map,aes(ProportionCountsMitochondrialGenes,SampleId)) +
    geom_rect(aes(xmin =PercProportionCountsMitochondrialGenes25 ,xmax =PercProportionCountsMitochondrialGenes75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercProportionCountsMitochondrialGenes50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_ProportionCountsMitochondrialGenes,color = "red",linetype = "longdash")+
    geom_line(group = 1)+
    #geom_bar(stat = "identity",,width = 1,fill = "grey30",color = NA) +
    facet_grid(Cluster~.,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3)+
    theme(
      legend.position = "top",
      panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(size = unit(0.1,"line")),
      axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      strip.text.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9)
    ) +
    xlab(label = "Proportion\nmitochondria") +
    scale_x_continuous(expand = c(0,0),limits = c(0,1),oob = squish,breaks = seq(0,1,0.1))
  
  ####
  df_ag <- Map[,c("SampleId","Cluster","ProteinCodingGenesDetected")] %>% unique %>%
    aggregate(ProteinCodingGenesDetected~Cluster,.,quantile)
  
  colnames(df_ag$ProteinCodingGenesDetected) <- colnames(df_ag$ProteinCodingGenesDetected) %>%
    gsub(pattern = "%",replacement = "") %>%
    gsub(pattern = "^",replacement = "PercProteinCodingGenesDetected")
  
  Map <- match(Map$Cluster,df_ag$Cluster) %>% df_ag[.,-1] %>% cbind(Map,.)
  
  
  p_prot <- ggplot(data = Map,aes(ProteinCodingGenesDetected,SampleId)) +
    geom_rect(aes(xmin =PercProteinCodingGenesDetected25 ,xmax =PercProteinCodingGenesDetected75,ymin = -Inf,ymax = Inf,group = Cluster),
              fill = "#A0CC2C",color = NA,alpha = 0.1)+
    geom_vline(aes(xintercept = PercProteinCodingGenesDetected50, group = Cluster), colour = '#DD14D3',
               size = 1.3) +
    geom_vline(xintercept = cutoff_ProteinCodingGenesDetected,color = "red",linetype = "longdash")+
    geom_line(group = 1)+
    #geom_bar(stat = "identity",,width = 1,fill = "grey30",color = NA) +
    facet_grid(Cluster~.,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 0.3)+
    theme(
      legend.position = "top",
      panel.grid.major.x = element_line(linetype = "dotted",color= "grey"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(size = unit(0.1,"line")),
      axis.text.x = element_text(size = 9,color = "grey30",angle = 90,vjust = 0.5,hjust = 1),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(color = 'black',size = 0.3),
      panel.spacing.y = unit(0.1, "lines"),
      axis.title.x = element_text(size = 9),
      strip.text.y = element_text(size = 9,angle = 0)
      
    ) +
    xlab(label = "Number\nprotein-coding genes") 
  
  #Append the dendrogram
  p_tree_cols <- ggtree(clust_genes, ladderize = F, size = 0.3) +
    coord_flip() +scale_x_reverse(expand =  c(0.001, 0.001)) 
  p_blank <- ggplot() + theme_void()
  
  #Bar###
  ########
  
  
  ### Check if metadata is passed and append it as a bar graph
  if(is.null(metadata_file)){
    
    composition <- egg::ggarrange(p_tree_cols,p_blank,p_blank,p_blank,p_blank,p_blank,
                                  p2,p_map,p_exon,p_inter,p_mit,p_prot,
                                  nrow = 2,ncol = 6,
                                  heights = c(relative_size_dendrogram,1),
                                  widths = c(relative_size_heatmap,relative_size_mappability,
                                             relative_size_exon,relative_size_intergenic,
                                             relative_size_chrM,relative_size_numprotgenes))
    dev.off()
    
    
    
  }else{
    
    #Create bar graph plotting group
    Map$Bar <- "Bar"
    
    p_bar <- ggplot(data = Map,aes(Bar,SampleId))+
      geom_raster(aes(fill = group))+
      facet_grid(Cluster~.,space = "free",scales = "free") +
      theme_ohchibi(size_panel_border = 0.3) +
      theme(
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 9,color = "grey30"),
        strip.text.y = element_text(size = 9,angle = 0),
        legend.text = element_text(size = 9,color = "grey30"),
        legend.title = element_text(size = 9),
        plot.title = element_text(size = 9),
        panel.spacing.y = unit(0.1, "lines"),
        
        
      )  +
      scale_x_discrete(expand = c(0,0)) +
      scale_fill_paletteer_d("ggthemes::Tableau_20",name = "Group") 
    
    
    #### Create composition
    composition <- egg::ggarrange(p_tree_cols,p_blank,p_blank,p_blank,p_blank,p_blank,p_blank,
                                  p2,p_map,p_exon,p_inter,p_mit,p_prot +theme(strip.text.y = element_blank()),p_bar,
                                  nrow = 2,ncol = 7,
                                  heights = c(relative_size_dendrogram,1),
                                  widths = c(relative_size_heatmap,relative_size_mappability,
                                             relative_size_exon,relative_size_intergenic,
                                             relative_size_chrM,relative_size_numprotgenes,relative_size_bar_plot_group))
    dev.off()
    
  }
  
  

  merged <- Mat %>% as.data.frame   %>% 
    tibble::rownames_to_column() %>% 
    dplyr::rename(.data =.,SampleId = rowname,
                  Verdict_PropMappability = V1,
                  Verdict_PropExonic = V2,
                  Verdict_PropIntergenic = V3,
                  Verdict_ProportionCountsMitochondrialGenes  = V4,
                  Verdict_ProteinCodingGenesDetected  = V5) %>%
    merge(df_clust,.,by = "SampleId") %>%
    dplyr::rename(.data =.,CompositeScore = Cluster)
  
  
  #Create a table summarizing how many fall per group 
  df_sum_cat <- merged$CompositeScore %>% table %>% data.frame %>%
    dplyr::rename(.data =.,Category = ".",NumberCells = Freq) %>%
    dplyr::mutate(.data =.,ProportionCells = ((NumberCells/sum(NumberCells))*100) %>% round(2))
  
  # Always add Group column and generate group summary since we now always have a group column
  merged$Group <- match(merged$SampleId,Map$SampleId) %>% Map$group[.]
  merged <- merged %>%
    dplyr::relocate(.data =.,c("SampleId","CompositeScore","Group"))

  df_sum_cat_group <- merged %>%
    dplyr::mutate(.data =.,Dummy =1) %>%
    reshape2::dcast(formula = Group~CompositeScore,fun.aggregate = sum,fill = 0,value.var = "Dummy")

  return(list(df_verdict = merged,
              df_sum_verdict = df_sum_cat,
              df_sum_verdict_group = df_sum_cat_group,
              composition = composition))
  
}

option_list <- list(
  make_option("--matrix_file", type = "character", help = "Expression matrix file"),
  make_option("--metrics_file", type = "character", help = "Metrics file"),
  make_option("--metadata_file", type = "character", default = NULL, help = "Metadata file (optional)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$matrix_file) || is.null(opt$metrics_file)) {
  stop("Please provide --matrix_file and --metrics_file")
}

matrix_file <- opt$matrix_file
metrics_file <- opt$metrics_file
metadata_file <- opt$metadata_file

res <- plot_qc_rna(matrix_file = matrix_file,
                   metrics_file = metrics_file,
                   metadata_file = metadata_file,
                   min_cells = 5,
                   top_genes = 500,
                   cutoff_PropMappability = 0.7,
                   cutoff_PropExonic = 0.7,
                   cutoff_PropIntergenic = 0.1,
                   cutoff_ProportionCountsMitochondrialGenes = 0.1,
                   cutoff_ProteinCodingGenesDetected = 500,
                   size_labels_samples = 1,
                   relative_size_heatmap = 1,
                   relative_size_mappability = 0.2,
                   relative_size_exon = 0.2,
                   relative_size_intergenic = 0.2,
                   relative_size_chrM = 0.2,
                   relative_size_numprotgenes = 0.2,
                   #Setting relative_size_dendrogram to 0 disables dendrogram on top
                   relative_size_dendrogram = 0,
                   relative_size_bar_plot_group = 0.03)

#We can save composition using oh.save.pdf function that enable shigh quality PDF saving
oh.save.pdf(p = res$composition,outname = "composition_rnaqc.pdf",outdir = ".",
            width = 16,height = 8)

# Save as JPEG using ggsave
ggplot2::ggsave(filename = "RNAQC_composition_mqc.jpg", plot = res$composition, width = 16, height = 8, units = "in", dpi = 300)


#Produces summarizing tables
res$df_sum_verdict %>% gt::gt()

res$df_sum_verdict_group %>% gt::gt()

write.table(res$df_verdict, "RNA-QC_ConsensusScores.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(res$df_sum_verdict,
            file = "summary_verdict.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

write.table(res$df_sum_verdict_group,
            file = "summary_verdict_group.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
