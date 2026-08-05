suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(parallel))
suppressMessages(library(parallelly))

#Create function to create segments
range_vectors <- function(x,y) {
return(seq(x,y,1))
}

#Function to decompose the read base information and extract correspondances 
dissect_read_bases <- function(string = NA){

    #Change \ for | for better parsing. \ represents quality and it is always preceded by ^
    x <- string %>% gsub(pattern = "\\\\",replacement = "|")%>% strsplit(split = "") %>% unlist  
    
    #Locate positions of $. $ represents the last position of a read 
    pos_dol <- x %>% grep(pattern = "\\$")

    #Locate ^ positions. ^ Represents first base and is followed by quality
    hat_loc <- x %>% grep(pattern = "\\^")    
    res_hat <- mapply(range_vectors,hat_loc, hat_loc+2)
    num_hat <- ncol(res_hat)
    pos_hat <- res_hat %>% unlist %>%  as.vector

    #Pos pure dol are end of read matches and the symbol we can omit
    #Here we differentiate these from the $ that is quality
    pos_pure_dol <- pos_dol[which(!(pos_dol %in% pos_hat))]

    #Locate here the positions that are indels
    minus_plus_loc <- x %>% grep(pattern = "[-+]")    

    #We differentiate them from the +- that represent quality
    minus_plus_loc <- minus_plus_loc[which(!(minus_plus_loc %in% pos_hat))]

    #Get the position of *. * represents pure deletions of the base
    star_loc <- x %>% grep(pattern = "\\*")    
    
    #Differentiate the  * from the * that represent quality , quality ones preceded by with ^
    pos_star <- star_loc[which(!(star_loc %in% pos_hat))] %>% unlist %>% as.vector 

    #Here check if we have any pure +- (insertions/deletions) events
    if(length(minus_plus_loc) > 0){

        #Here evaluate which positions +2 and +3 are numeric and which ones only have one symbol
        minus_plus_loc_three <- minus_plus_loc[which(sapply(x[minus_plus_loc+3],function(x)x %in% 0:9)) %>% as.vector]

        minus_plus_loc_two <- minus_plus_loc[which(sapply(x[minus_plus_loc+2],function(x)x %in% 0:9)) %>% as.vector]

        minus_plus_loc_one <- minus_plus_loc[which(!(minus_plus_loc %in% c(minus_plus_loc_two,minus_plus_loc_three)))] 

        minus_plus_num_bases_add_three <- paste0(x[minus_plus_loc_three+1],x[minus_plus_loc_three+2],x[minus_plus_loc_three+3]) %>% as.numeric
        
        minus_plus_num_bases_add_two <- paste0(x[minus_plus_loc_two+1],x[minus_plus_loc_two+2]) %>% as.numeric

        minus_plus_num_bases_add_one <- x[minus_plus_loc_one+1] %>% as.numeric

        minus_plus_loc <- c(minus_plus_loc_one,minus_plus_loc_two,minus_plus_loc_three)

        minus_plus_num_bases_add <- c(minus_plus_num_bases_add_one+1,minus_plus_num_bases_add_two+2,minus_plus_num_bases_add_three+3) %>% na.omit %>% as.vector

        res_mp <- mapply(range_vectors,minus_plus_loc-1, minus_plus_loc+minus_plus_num_bases_add)
    

        if(class(res_mp)[1] == "list"){

            num_mp <- length(res_mp)
            pos_mp <- res_mp %>% unlist %>% as.vector

        }else{

            num_mp <- ncol(res_mp)
            pos_mp <- res_mp %>% as.vector


        }

    }else{
        
        #Case in which we dont have any +/- events
        num_mp <- 0 
        pos_mp <- NULL

    }

    #There is a case in which ^ can be before deletion or insertion so if there is that overlap that corresponds to the deletion or insertion
    #Here check the overlaps and remove the procedent two characters from the vector it is always ^QOverlap
    if(is.null(ncol(res_hat)) == FALSE){

        todrop_hats <- NULL
        for(co in 1:ncol(res_hat)){

                vec_pos_hat <- res_hat[,co]
                res_inter <- which(vec_pos_hat %in% pos_mp)
                if(length(res_inter) > 0){
                    todrop_hats <- c(todrop_hats,co)

                }
        }
        
        if(length(todrop_hats) > 0){
            
            res_hat <- res_hat[,-todrop_hats]



        }

        #Clean hat loc should be first row of remaining
        if(is.null(nrow(res_hat))){
            
            res_hat <- res_hat %>% matrix(nrow = 3)

        }
        hat_loc <- res_hat[1,]
    }

    num_star <- length(pos_star)

    #Add unique sort to deal with the overlap of pos hat with insertions and deletions 
    pos_pec <- c(pos_mp,pos_hat,pos_star,pos_pure_dol) %>% unique %>% sort

    num_norm <- length(x)-length(pos_pec)

    num_tot <- c(num_norm,num_mp,num_hat,num_star)  %>% sum

    #Indices norm represent the indices for snps
    indices_snps <- which(!(seq(1,length(x),1) %in% pos_pec)) 

    #Append the indices of the hat
    indices_snps <- c(indices_snps,hat_loc+2)

    #Get the reference positions 
    indices_ref_forward <- indices_snps[ toupper(x[indices_snps]) %>% grep(pattern = "[\\.]")]

    indices_ref_reverse <- indices_snps[ toupper(x[indices_snps]) %>% grep(pattern = "[,]")]
    
    #Store positions in a 4 column data frame (Allele, Type,Indices and Unique Id)
    df_res <- NULL
    if(length(indices_ref_forward) > 0){

        df_res <- data.frame(Allele = "REF",Type = "Forward",Indices = indices_ref_forward,UId = paste0("REF_Forward_",1:length(indices_ref_forward))) %>%
            rbind(df_res,.)

    }

    if(length(indices_ref_reverse) > 0){

        df_res <- data.frame(Allele = "REF",Type = "Reverse",Indices = indices_ref_reverse,UId = paste0("REF_Reverse",1:length(indices_ref_reverse))) %>%
            rbind(df_res,.)

    }

    #Now get all SNPs
    sorted_alt <- x[indices_snps] %>% grep(pattern ="[\\.,]",invert = TRUE,value = TRUE)  %>% toupper %>% table %>% sort(decreasing = TRUE)

    if(length(sorted_alt) > 0){

        #Loop over the alleles and get the information
        for(ntalt in names(sorted_alt)){

            indices_alt_reverse <- indices_snps[ (x[indices_snps]) %>% grep(pattern = tolower(ntalt))]
            indices_alt_forward <- indices_snps[ (x[indices_snps]) %>% grep(pattern = toupper(ntalt))]

            if(length(indices_alt_forward) >0){

                df_res <- data.frame(Allele = ntalt,Type = "Forward",Indices = indices_alt_forward,UId = paste0(ntalt,"_Forward_",1:length(indices_alt_forward))) %>%
                    rbind(df_res,.)

            }

            if(length(indices_alt_reverse) >0){

                df_res <- data.frame(Allele = ntalt,Type = "Reverse",Indices = indices_alt_reverse,UId = paste0(ntalt,"_Reverse_",1:length(indices_alt_reverse))) %>%
                    rbind(df_res,.)

            }
            

        }
    }

    #Process deletions and insertions 
    if(num_mp > 0){

        if(class(res_mp)[1] == "list"){

            #Populate df_res
            for(z in 1:length(res_mp)){

                mindices <- res_mp[[z]]
                mchars <- x[mindices]

                nom <- mchars[-1] %>% paste0(collapse = "") %>% toupper()
                type_alt <- "Forward"
                if(mchars[1] == ","){

                    type_alt <- "Reverse"

                }
                df_res <- data.frame(Allele = nom,Type = type_alt,Indices = mindices,UId = paste0(nom,"_",type_alt,"_",z)) %>%
                    rbind(df_res,.)

            }


        }else{

            
            #Create a matrix of characters fro res_mp
            res_mp_char <- apply(res_mp,MARGIN=2,FUN=function(y)x[y] %>% toupper)
            
            #Append to df_res structure
            for(z in 1:ncol(res_mp_char)){
                
                nom <- res_mp_char[-1,z] %>% paste0(collapse = "")
                type_alt <- "Forward"
                if(res_mp_char[1,z] == ","){

                    type_alt <- "Reverse"

                }
                df_res <- data.frame(Allele = nom,Type = type_alt,Indices = res_mp[,z],UId = paste0(nom,"_",type_alt,"_",z)) %>%
                    rbind(df_res,.)

            }
        }
        
    }

    #Process deletions coded as *
    if(length(pos_star) > 0){

        df_res <- data.frame(Allele = "DELREFCIGAR",Type = "Forward",Indices = pos_star,UId = paste0("DELREFCIGAR_FWD_",1:length(pos_star))) %>%
            rbind(df_res,.)
            

    }

    ### Perform translation to original positions
    df_res <- with(df_res,order(Indices)) %>% df_res[.,]

    df_res$UId <- df_res$UId %>% factor(levels = df_res$UId %>% as.character %>% unique)

    df_res$NumPos <- df_res$UId %>% as.numeric

    df_ag <- df_res[,c("Allele","Type","NumPos")] %>% unique %>%
        aggregate(NumPos~Allele+Type,.,FUN = function(x)sort(x) )

    #Return list with allele name, indices forward, indices reverse
    list_res <- list()
    for(y in 1:nrow(df_ag)){

        mline <- df_ag[y,]
        mnom <- mline[1]%>% unlist %>% as.character
        mtype <- mline[2]%>% unlist %>% as.character
        mvec <- mline[3]%>% unlist %>% as.vector

        list_res[[mnom]][[mtype]] <- mvec

    }

    return(list_res)

}

dissect_AS <- function(string = NA){

    return(strsplit(x = string,split = ",") %>% unlist %>% as.numeric)

}

dissect_MQ_BQ <- function(string = NA){

    return(string %>% strsplit(split = "") %>% unlist %>% sapply(FUN = function(x)utf8ToInt(x)-33)  %>% as.vector)

}

dissect_BPSTART <- function(string = NA){

    return(strsplit(x = string,split = ",") %>% unlist %>% as.numeric)

}

dissect_READ_IDS <- function(string = NA){

    return(strsplit(x = string,split = ",") %>% unlist)

}

pileup_by_row_to_parsed_dataframe <- function(row_to_parse = NA,sample_name = NA,threshold_mq = 30,threshold_bq = 30,threshold_bp = 15,read_length = 150){
    #To process ZA
    if(is.na(row_to_parse$READ_BASES) || row_to_parse$READ_BASES == ""){
        row_to_parse$READ_BASES <- "ZA"
    }
    list_rd_temp <- dissect_read_bases(row_to_parse$READ_BASES)
    vec_as_temp <- dissect_AS(row_to_parse$AS) %>% as.numeric
    vec_mq_temp <- dissect_MQ_BQ(row_to_parse$MQ) %>% as.numeric
    vec_bq_temp <- dissect_MQ_BQ(row_to_parse$BQ) %>% as.numeric
    vec_bstart_temp <- dissect_BPSTART(row_to_parse$BPSTART) %>% as.numeric
    vec_cigars_temp <- dissect_READ_IDS(row_to_parse$CIGAR)

    #Logic to decompose the list of alleles found
    malleles <- names(list_rd_temp)

    tot_depth_all <- 0
    tot_depth_hq <- 0

    vec_num_forward_all <- NULL
    vec_num_reverse_all <- NULL

    vec_median_focal_as <- NULL
    vec_prop_clipped <- NULL

    vec_num_forward_hq<- NULL
    vec_sd_forward_hq <- NULL
    vec_mad_forward_hq <- NULL
    vec_prop_forward_hq <- NULL
    vec_num_reverse_hq <- NULL
    vec_sd_reverse_hq <- NULL
    vec_mad_reverse_hq <- NULL
    vec_prop_reverse_hq <- NULL
    vec_prop_upper_forward_hq <- NULL
    vec_prop_upper_reverse_hq <- NULL

    vec_depth_all_position <- NULL
    vec_depth_hq_position <- NULL

    for(mallele in malleles){

        #print(mallele)

        indices_forward <- list_rd_temp[[mallele]]$Forward
        indices_reverse <- list_rd_temp[[mallele]]$Reverse

        both_indices <- c(indices_forward,indices_reverse)

        #Calculate AS median
        median_as <- vec_as_temp[both_indices] %>% median

        #Proportion of clipped reads
        prop_bases_clipped <- length(vec_cigars_temp[both_indices] %>% grep(pattern = "H|S"))/(length(both_indices))

        #Number of MQ>=30 BQ>=30 bases supporting read        
        forward_indices_hq_mq_bq <- intersect(indices_forward[which(vec_mq_temp[indices_forward] >= threshold_mq)],indices_forward[which(vec_bq_temp[indices_forward] >= threshold_bq)])
        reverse_indices_hq_mq_bq <- intersect(indices_reverse[which(vec_mq_temp[indices_reverse] >= threshold_mq)],indices_reverse[which(vec_bq_temp[indices_reverse] >= threshold_bq)])


        ### Calculate the base SD and MAD using all MQ30 BQ30 reads supporting a variant         
        bpstart_hq_mq_bq_forward <- vec_bstart_temp[forward_indices_hq_mq_bq]
    
        bpstart_hq_mq_bq_reverse <- vec_bstart_temp[reverse_indices_hq_mq_bq]
    
        prop_under_forward <- length(which(bpstart_hq_mq_bq_forward < threshold_bp))/length(bpstart_hq_mq_bq_forward)

        prop_under_reverse <- length(which(bpstart_hq_mq_bq_reverse < threshold_bp))/length(bpstart_hq_mq_bq_reverse)

        up_limit <- read_length - threshold_bp

        prop_upper_forward <- length(which(bpstart_hq_mq_bq_forward > up_limit))/length(bpstart_hq_mq_bq_forward)

        #print(bpstart_hq_mq_bq_forward)
        
        prop_upper_reverse <- length(which(bpstart_hq_mq_bq_reverse > up_limit))/length(bpstart_hq_mq_bq_reverse)

        #print(bpstart_hq_mq_bq_reverse)

    
        sd_forward_hq_mq_bq<- sd(bpstart_hq_mq_bq_forward,na.rm = TRUE)
    
        mad_forward_hq_mq_bq <- mad(bpstart_hq_mq_bq_forward,na.rm = TRUE)
    
        sd_reverse_hq_mq_bq <- sd(bpstart_hq_mq_bq_reverse,na.rm = TRUE)
        
        mad_reverse_hq_mq_bq <- mad(bpstart_hq_mq_bq_reverse,na.rm = TRUE)

        #Append to vectors holding all metrics
        vec_num_forward_all <- c(vec_num_forward_all,length(indices_forward))
        vec_num_reverse_all <- c(vec_num_reverse_all,length(indices_reverse))

        vec_median_focal_as <- c(vec_median_focal_as,median_as)
        vec_prop_clipped <- c(vec_prop_clipped,prop_bases_clipped)

        vec_num_forward_hq <- c(vec_num_forward_hq,length(forward_indices_hq_mq_bq))
        vec_prop_forward_hq <- c(vec_prop_forward_hq,prop_under_forward)
        vec_prop_upper_forward_hq <- c(vec_prop_upper_forward_hq,prop_upper_forward)
        vec_sd_forward_hq <- c(vec_sd_forward_hq,sd_forward_hq_mq_bq)
        vec_mad_forward_hq <- c(vec_mad_forward_hq,mad_forward_hq_mq_bq)

        vec_num_reverse_hq <- c(vec_num_reverse_hq,length(reverse_indices_hq_mq_bq))
        vec_prop_reverse_hq <- c(vec_prop_reverse_hq,prop_under_reverse)
        vec_prop_upper_reverse_hq <- c(vec_prop_upper_reverse_hq,prop_upper_reverse)
        vec_sd_reverse_hq <- c(vec_sd_reverse_hq,sd_reverse_hq_mq_bq)
        vec_mad_reverse_hq <- c(vec_mad_reverse_hq,mad_reverse_hq_mq_bq)

        tot_depth_all <- tot_depth_all + length(both_indices)
        tot_depth_hq <- tot_depth_hq + (length(forward_indices_hq_mq_bq) + length(reverse_indices_hq_mq_bq))


    }

    df_temp <- data.frame(
                        CHROM = rep(row_to_parse$CHROM,length(malleles)),
                        POS = rep(row_to_parse$POS,length(malleles)),
                        REF = rep(row_to_parse$REF,length(malleles)),
                        ALT = malleles,
                        NUM_FRAGMENTS_ALLQ_MQ_BQ_F = vec_num_forward_all,
                        NUM_FRAGMENTS_ALLQ_MQ_BQ_R = vec_num_reverse_all,
                        NUM_FRAGMENTS_HQ_MQ_BQ_F = vec_num_forward_hq,
                        NUM_FRAGMENTS_HQ_MQ_BQ_R = vec_num_reverse_hq,
                        MEDIAN_AS_VARIANT_READS = vec_median_focal_as,
                        PROP_BASES_CLIPPED = vec_prop_clipped,
                        PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_F = vec_prop_forward_hq,
                        PROP_FRAGMENTS_BPSTART_UNDER_HQ_MQ_BQ_R = vec_prop_reverse_hq,
                        SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_F = vec_sd_forward_hq,
                        SD_BPSTART_FRAGMENTS_HQ_MQ_BQ_R = vec_sd_reverse_hq,
                        MAD_BPSTART_RAGMENTS_HQ_MQ_BQ_F = vec_mad_forward_hq,
                        MAD_BPSTART_RAGMENTS_HQ_MQ_BQ_R = vec_mad_reverse_hq,
                        NUM_FRAGMENTS_ALLQ_POSITION = rep(tot_depth_all,length(malleles)),
                        NUM_FRAGMENTS_HQ_POSITION = rep(tot_depth_hq,length(malleles)),
                        PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_F = vec_prop_upper_forward_hq,
                        PROP_FRAGMENTS_BPSTART_UPPER_HQ_MQ_BQ_R = vec_prop_upper_reverse_hq
                )

    indices_ins <- df_temp$ALT %>% grep(pattern = "\\+") 
    indices_del <- df_temp$ALT %>% grep(pattern = "\\-") 
    indices_ref <-  df_temp$ALT %>% grep(pattern = "^REF$") 
    indices_delrefcigar <-  df_temp$ALT %>% grep(pattern = "^DelREFCIGAR$") 
    indices_snps <- (1:nrow(df_temp))[which(!(1:nrow(df_temp) %in% c(indices_ins,indices_del,indices_ref,indices_delrefcigar)))]

    df_snps <- NULL
    df_ref <- NULL
    df_ref_dc <- NULL
    df_ins <- NULL
    df_del <- NULL

    if(length(indices_snps) >0){
        df_snps <- df_temp[indices_snps,]
        df_snps$VariantId <- paste0(df_snps$CHROM,"_",df_snps$POS,"_",df_snps$REF,"_",df_snps$ALT)
    }

    if(length(indices_ref) >0){
        df_ref <- df_temp[indices_ref,]
        df_ref$VariantId <- paste0(df_ref$CHROM,"_",df_ref$POS,"_",df_ref$REF,"_",df_ref$REF)
    }

    if(length(indices_delrefcigar) >0){
        df_ref_dc <- df_temp[indices_delrefcigar,]
        df_ref_dc$VariantId <- paste0(df_ref_dc$CHROM,"_",df_ref_dc$POS,"_",df_ref_dc$REF,"_DELREFCIGAR")
    }

    if(length(indices_ins) >0){
        df_ins <- df_temp[indices_ins,]
        df_ins$VariantId <- paste0(df_ins$CHROM,"_",df_ins$POS,"_",df_ins$REF,"_",paste0(df_ins$REF,"",gsub(x =df_ins$ALT,pattern = "[\\+][0-9]+",replacement = "")))
    }

    if(length(indices_del) >0){
        df_del <- df_temp[indices_del,]
        df_del$VariantId <- paste0(df_del$CHROM,"_",df_del$POS,"_",paste0(df_del$REF,"",gsub(x =df_del$ALT,pattern = "[\\-][0-9]+",replacement = "")),"_",df_del$REF)
    }

    df_temp <- rbind(df_ref,df_snps,df_ref_dc,df_ins,df_del) %>% dplyr::arrange(.data =.,POS)

    df_temp$SampleId <- sample_name
    
    df_temp <- df_temp %>% dplyr::relocate(.data =.,"SampleId","VariantId") 

    return(df_temp)
}

#Read arguments
args = commandArgs(trailingOnly=TRUE)
pileup_cigars_file <- args[1]
sample_name <- args[2]
threshold_mq <- as.numeric(args[3])
threshold_bq <- as.numeric(args[4])
threshold_bp <- as.numeric(args[5])
chunk_size <- as.numeric(args[6])
read_length <- as.numeric(args[7])

### Parallelization ###
numCores <- max(1, parallelly::availableCores())
cl <- makeCluster(numCores,type = "FORK")
registerDoParallel(numCores) 

cat("Number of cores to use: ",numCores,"\n")

output_file <- "res.tsv"

#Process each chunk
chunk = 1
done = FALSE
while(!done){

    cat("Working on chunk ",chunk," of size ",chunk_size,"\n")

  #Read chunk — force all columns to character to prevent fread type inference
    #from reading MQ/AS/BPSTART as integers (happens when a sample has only 1-2 pileup rows)
    chunk_data <- tryCatch(fread(pileup_cigars_file,sep = "\t", header = FALSE, data.table = FALSE,quote = "",skip=(chunk-1)*chunk_size,nrows=chunk_size-1,
            colClasses = list(character = 1:12)),
            error = function(e)return(NA))

    if (class(chunk_data) != "data.frame"){
        done = TRUE
        break
    } 

    chunk = chunk + 1

    colnames(chunk_data) <- c("CHROM","POS","REF","NREADS","READ_BASES","BQ","MQ","READS","FLAGS","BPSTART","AS","CIGAR")

    # Parallel processing
    chunk_results <- foreach(i = 1:nrow(chunk_data), .combine = rbind, .packages = 'dplyr') %dopar% {
        pileup_by_row_to_parsed_dataframe(row_to_parse = chunk_data[i,],sample_name =sample_name,threshold_mq = threshold_mq,threshold_bq = threshold_bq,threshold_bp = threshold_bp,read_length = read_length )
    }
    # Append results to file
    fwrite(chunk_results, output_file, append = TRUE, sep = "\t", col.names = FALSE)
}

stopCluster(cl)
