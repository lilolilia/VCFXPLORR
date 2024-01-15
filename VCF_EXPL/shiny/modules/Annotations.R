library("biomaRt")
library('ReactomePA')
library('org.Hs.eg.db')

#This file gather the functions to annotate the vcf

#FILTER THE SNPS DATABASE BASED ON THE SELECTED CHROMOSOME

GET_SNPS<- function(data){
  tar_chr <- as.vector(seqnames(data)@values)
  tar_chr <- gsub("chr", "", tar_chr)
  tar_chr[grepl(tar_chr, pattern = "M")] <- "MT"
  my_snps <- snpsBySeqname(all_snps, c(tar_chr))
  seqlevelsStyle(my_snps) <- "UCSC"
  genome(my_snps) <- "hg38"
  return(my_snps)
} #end of function

#MATS 

GET_MATS <- function(data,snp_id){
  matV1<-as.data.frame(geno(data)$AF)
  colnames(matV1)<-"AF"
  matV1$Variant<-rownames(matV1)
  # matV1 <- data.frame(Variant = names(data), stringsAsFactors = FALSE)
  matV1$chromosome <- gsub("(.*):(.*)_(.*)/(.*)", "\\1", matV1$Variant)
  matV1$start <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
  matV1$end <- gsub("(.*):(.*)_(.*)/(.*)", "\\2", matV1$Variant)
  matV1$ref_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\3", matV1$Variant)
  matV1$alt_allele <- gsub("(.*):(.*)_(.*)/(.*)", "\\4", matV1$Variant)
  matV1$posIDX <- gsub("(.*)_(.*)", "\\1", matV1$Variant)
  matS <- merge(matV1,snp_id,all.x=TRUE,by="posIDX")
  matS <- dplyr::select(matS,-posIDX)
  return(matS)
} #end of function

#annotations

GET_ANNOTATIONS <- function(data){
  ensinfo<- getBM(attributes=c("refsnp_id", "ensembl_gene_stable_id",
                               "sift_score","associated_gene","clinical_significance",
                               "minor_allele_freq"),
                filters="snp_filter", values=(data$rsID),
                mart=ensembl, uniqueRows=TRUE)
  theBM <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol", "entrezgene_id",
                              "entrezgene_description", "p_value"), 
              filters = 'ensembl_gene_id', 
              values =  ensinfo$ensembl_gene_stable_id, 
              mart = ensembl1)
  names(ensinfo)[names(ensinfo) == 'ensembl_gene_stable_id'] <- "ensembl_gene_id"
  ensinfo<-as.data.frame(ensinfo)
  theBM<-as.data.frame(theBM)
  df<-merge(theBM,ensinfo,all.x=TRUE,by="ensembl_gene_id")
  names(df)[names(df) == 'refsnp_id'] <- "rsID"
  df2<-merge(data,df,all.x=TRUE,by="rsID")
  return(df2)
}

#reactome pathways

GET_PATHWAYS <- function(data){
  # Extract unique Ensembl Gene IDs
  gene_ids <- data$entrezgene_id[!is.na(data$entrezgene_id)]
  # Run pathway enrichment analysis
  enrich_result <- enrichPathway(gene = gene_ids)
  # Print the enrichment result
  results<- as.data.frame(enrich_result@result)
  return(results)
}

#MATA

GET_MATA <- function(data){
  subset_vcf <- data
  seqlevels(subset_vcf, pruning.mode="coarse") <- target_chromosomes
  coding <- predictCoding(subset_vcf, txdb, seqSource = Hsapiens)
  matA <- data.frame(Variant=names(coding),
                     chromosome=seqnames(coding),
                     start=start(coding),end=end(coding),
                     ref_allele=as.character(coding$REF),
                     alt_allele=unlist(lapply(lapply(
                       coding$ALT,`[[`,1),as.character)),
                     GeneID=coding$GENEID,
                     TxID=coding$TXID,
                     Protein_posi=unlist(lapply(lapply(
                       coding$PROTEINLOC,`[[`,1),as.integer)),
                     ref_AA=as.character(coding$REFAA),
                     alt_AA=as.character(coding$VARAA),
                     Type=coding$CONSEQUENCE)
  matA$aaChange <- paste0("p.",matA$ref_AA,matA$Protein_posi,matA$alt_AA)
  matA <- dplyr::select(matA,-Protein_posi,-ref_AA,-alt_AA)
  return(matA)
} #end of function

#TACDAT

GET_TACDT<-function(data, matA){
  var_in_coding <- data.frame(varName = names(data), 
                              in_coding = names(data) %in%
                                matA$Variant, stringsAsFactors = FALSE)
  taC <- table(matA$Type)
  taC_dat <- as.data.frame(taC)
  return(taC_dat)
} #end of function
