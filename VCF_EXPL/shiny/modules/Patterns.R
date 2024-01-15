#THIS FILE GATHERS OTHER USEFUL FUNCTIONS TO IDENTIFY PATTERNS

#GET THE PATTERNS IN THE VCF FILES

GET_PATTERNS <-function(data){
  rd_idx <- str_split(names(data), "_", simplify = T)
  rd_sub <- data[rd_idx[, 2] == "C/T"]
  #extract sequence beneath mutation
  rd_sub$triNu <- getSeq(Hsapiens,
                         seqnames(rd_sub),
                         start = pmax(1, start(rd_sub) - 1),
                         end = end(rd_sub) + 1)
  tbl <- table(rd_sub$triNu)
  tbl_dat <- as.data.frame(tbl)
  tbl_dat$APOBEC_target <- tbl_dat$Var1 %in% c("TCA", "TCT")
  return(tbl_dat)
  } #end of function
