library("fastmap")
library('VariantAnnotation')
library('ggplot2')

#This function is designed to have the variants corresponding to the ids.

Get_Variants <- function(sample1DT) {
  
  rd <- rowRanges(sample1DT)
  
  var_1 <- rownames(geno(sample1DT)$GT)[
    geno(sample1DT)$GT %in% c("0/1", "1/1", "0|1", "1|1")]
  
  if (!is.null(var_1)) {
    
    varTab1 <- data.frame(variant=names(rd)[names(rd) %in% var_1],
                          chr=as.vector(seqnames(rd)[names(rd) %in% var_1]),
                          start=start(rd)[names(rd) %in% var_1],
                          end=end(rd)[names(rd) %in% var_1],
                          stringsAsFactors = FALSE)
    
    #get their ref base
    ref_base <- ref(sample1DT)[rownames(sample1DT) %in% var_1]
    varTab1$refBase <- as.character(ref_base)
    
    #alt base 
    alt_base <- lapply(alt(sample1DT)[rownames(sample1DT) %in% var_1],`[[`,1)
    alt_base <- lapply(alt_base,as.character)
    varTab1$altBase <- unlist(alt_base)
    
    #count from AD 
    
    adCount <-geno(sample1DT)$AD[rownames(geno(sample1DT)$AD) %in% var_1]
    varTab1$refCount <- unlist(lapply(adCount,`[[`,1))
    varTab1$altCount <- unlist(lapply(adCount,`[[`,2))
    varTab1$genoType <- geno(sample1DT)$GT[rownames(geno(sample1DT)$GT) %in% var_1]
    varTab1$gtQuality <- geno(sample1DT)$GQ[rownames(geno(sample1DT)$GQ) %in% var_1]
    
  }

  else {

    varTab1 <- NULL
  }

  #get GT for  1/2

  var_2 <- rownames(geno(sample1DT)$GT)[
    geno(sample1DT)$GT %in% c("1/2", "1|2")]

  if(!is.null(var_2)){

    varTab2 <- data.frame(variant=names(rd)[names(rd) %in% var_2],
                        chr=as.vector(seqnames(rd)[names(rd) %in% var_2]),
                        start=start(rd)[names(rd) %in% var_2],
                        end=end(rd)[names(rd) %in% var_2],
                        refBase=unlist(lapply(lapply(
                          alt(sample1DT)[rownames(sample1DT) %in% var_2],`[[`,1),as.character)),
                        altBase=unlist(lapply(lapply(
                          alt(sample1DT)[rownames(sample1DT) %in% var_2],`[[`,2),as.character)),
                        refCount=unlist(lapply(
                          geno(sample1DT)$AD[rownames(geno(sample1DT)$AD) %in% var_2],`[[`,2)),
                        altCount=unlist(
                          lapply(geno(sample1DT)$AD[rownames(geno(sample1DT)$AD) %in% var_2],`[[`,3)),
                        genoType=geno(sample1DT)$GT[rownames(geno(sample1DT)$GT) %in% var_2],
                        gtQuality=geno(sample1DT)$GQ[rownames(geno(sample1DT)$GQ) %in% var_2],
                        stringsAsFactors = FALSE)
  }

  else {
    varTab2 <- NULL
  }

  # Combine the two variant tables
  varTab <- if (!is.null(varTab1) && is.null(varTab2)) {
    varTab1
  } else if (is.null(varTab1) && !is.null(varTab2)) {
    varTab2
  } else {
    rbind(varTab1, varTab2)
  }
  
  for (k in 1:length(varTab$variant)) {
      if (width(varTab$refBase[k]) < width(varTab$altBase[k])) {
        varTab$mutType[k] <- "INS"
      } else if (width(varTab$refBase[k]) > width(varTab$altBase[k])) {
        varTab$mutType[k] <- "DEL"
      } else if (width(varTab$refBase[k]) == 1 & width(varTab$altBase[k]) == 1) {
        varTab$mutType[k] <- "SNP"
      } else {
        varTab$mutType[k] <- "Others"
      }
  }
  
  varTab$nuSub <- paste0(varTab$refBase,">",varTab$altBase)
  varTab$TiTv[varTab$nuSub %in% ti] <- "Ti"
  varTab$TiTv[varTab$nuSub %in% tv] <- "Tv"

  return(varTab)
}