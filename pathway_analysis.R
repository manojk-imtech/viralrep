# Loading required libraries
library(pheatmap)
library(RColorBrewer)

# Input location of RData file
cat("Enter location of working directory: ")
wd <- readLines(file("stdin"),1)
setwd(wd)
load("./req_data.RData")

# Mapping uniprot IDs to gene entrez IDs for all viruses
kegg_count_tab <- list()
for(i in 1:length(gene2kegg)){  
  map <- gene2kegg[[i]]
  genes <- as.numeric(unlist(sapply(map$Uniprot_id, function(x){
    as.numeric(uniprot_mapping[(which(uniprot_mapping[,1]==as.character(x) & uniprot_mapping[,2]=="GeneID")),3])
  })))
  uni_genes <- unique(genes)
  df <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = F)
  
  for(j in uni_genes){
    for(k in 1:length(genes.by.pathway)){
      if(j%in%as.numeric(genes.by.pathway[[k]])){
        df <- rbind(df, c(j,names(genes.by.pathway[k])), stringsAsFactors = FALSE, factor.exclude	= NA)
      }
      
    }
  }
  
  res <- keggList("pathway", "hsa")
  pathway.names <- unname(res[paste0("path:",df[,2])])
  df <- cbind(df, pathway.names)
  df <- cbind(df, as.numeric(table(df[,3])[df[,3]]))
  colnames(df) <- c("gene_entrez_id", "kegg_pathway_code", "kegg_pathway_name", "kegg_pathway_count")
  
  count_tab <- data.frame(kegg_pathways=names(table(df[,3])), count=as.numeric(table(df[,3])))
  count_tab <- count_tab[order(count_tab$count, decreasing = T),]
  
  kegg_count_tab[[i]] <- count_tab 
  
}

# Generating a collection of unique pathways which are present in the top 10 for any virus
all_pathways <- sapply(kegg_count_tab, function(x){
  x[1:10,"kegg_pathways"]
})

sel_pathways <- unique(as.vector(all_pathways))[!is.na(unique(as.vector(all_pathways)))]

# Generating the data matrix for the heatmap
hm_data <- data.frame(matrix(ncol=13,nrow = length(sel_pathways)))
rownames(hm_data) <- sel_pathways

for(i in 1:ncol(hm_data)){
  x <- kegg_count_tab[[i]]
  if(dim(x)[1]>=10){
    
    hm_data[x[which(x[,1]%in%rownames(hm_data)),1],i] <-x[which(x[,1]%in%rownames(hm_data)),2]
  }else{
    hm_data[x[1:dim(x)[1],1],i] <-x[1:dim(x)[1],2]
  }
}

rownames(hm_data) <- sapply(strsplit(sel_pathways, split =" - Homo"), "[[", 1)

colnames(hm_data) <- c("CHIKV", "CCHF", "EBV", "HeV", "IAV/IBV", "LASV", "MARV", "MERS", "nCov-19", "NiV", "RVFV", "SARS", "ZIKV")

p_codes <- sapply(as.character(rownames(hm_data)), function(x){
  names(pathways.list[grep(paste0(x, " - Homo sapiens (human)"), pathways.list, fixed = T, ignore.case = F)])
})
clean_p_codes <- unname(sapply(strsplit(p_codes, split = ":", fixed=T), "[[", 2))
gene_in_pathway <- sapply(clean_p_codes, function(x){
  length(genes.by.pathway[[x]])
})

rownames(hm_data) <- paste0(rownames(hm_data)," [",gene_in_pathway,"]")

# Code for plotting the heatmap
#tiff("./pathway_heatmap_all.tiff", height = 900, width = 800)
pheatmap(hm_data[names(sort(rowSums(hm_data, na.rm = T),decreasing =F)), ], cluster_rows = F, cluster_cols = F, na_col = "lightgrey", color = brewer.pal(9, "Blues"), angle_col = 90, fontsize_col = 18, fontsize_row = 12)
#dev.off()

