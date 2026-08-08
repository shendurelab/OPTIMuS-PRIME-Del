library(Seurat)
library(stringr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(BPCells)
library(rtracklayer)

obj <- readRDS("./K562-L1-L2.filtered-normalized-obj.rds")

sites_granges <- function(sites) {
    
    # Split the input string and handle cases with or without the "end" value
    site_df <- data.frame(do.call(rbind, lapply(strsplit(sites, "-"), function(x) {
        if (length(x) == 2) {
            x <- c(x, x[2])  # Add the start value as the end value
        }
        return(x)
    })))

    colnames(site_df) <- c("seqnames", "start", "end")

    site_df$start <- as.numeric(site_df$start)
    site_df$end <- as.numeric(site_df$end)
    
    site_gr <- makeGRangesFromDataFrame(site_df)
    site_gr
    
}


CRE.gr <- sites_granges(unique(obj@meta.data[obj@meta.data$CRE != "" &
                               grepl("_",obj@meta.data$CRE) == FALSE,]$CRE))

CRE.gr.expanded <- resize(CRE.gr,
                          width = width(CRE.gr) + 2000,
                          fix = "center")

CRE.df <- as.data.frame(CRE.gr.expanded)[,c("seqnames","start","end")]

write.table(
  CRE.df,
  file = "blacklist_t7_ranges.bed",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

