library(ohchibi)

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