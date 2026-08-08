library(Seurat)
library(stringr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(rtracklayer)
library(GenomicRanges)
library(Signac)
library(AnnotationHub)
library(reshape2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
library(patchwork)
library(ggrastr)
library(cowplot)

#############################
# reading in the processed 
# Seurat object
#############################
obj <- readRDS("./t7-filtered.K562-L1-L2.normalized-obj.rds")

#############################
# getting the number of 
# cells assigned per
# deletion
#############################
n_cells <- as.data.frame(matrix(ncol=0,nrow=0))

CREs <- unique(obj@meta.data[obj@meta.data$CRE_id != "" &
              grepl("_",obj@meta.data$CRE_id) == FALSE,]$CRE_id)
for (CRE in CREs) {
        n_cells <- rbind(n_cells,c(CRE,nrow(obj@meta.data[grepl(CRE,obj@meta.data$CRE_id) == TRUE,])))
}

colnames(n_cells) <- c("CRE","n.cells")
n_cells$n.cells <- as.numeric(n_cells$n.cells)
n_cells <- n_cells[order(-n_cells$n.cells),]
n_cells$CRE <- factor(n_cells$CRE, levels = n_cells$CRE)

#############################
# function(s) to get genes
# close to each targeted
# enhancer
#############################
# Build transcript-level TSS once (from your same GENCODE GTF)
make_transcript_tss <- function(gtf) {
  tx <- gtf[gtf$type == "transcript"]
  biotype_col <- if ("gene_type" %in% names(mcols(tx))) "gene_type" else "gene_biotype"
  tx <- tx[mcols(tx)[[biotype_col]] == "protein_coding"]
  tx <- keepStandardChromosomes(tx, pruning.mode = "coarse")
  seqlevelsStyle(tx) <- "UCSC"
  tss_tx <- promoters(tx, upstream = 1, downstream = 1)
  # unpack tags into a simple character vector when present
  if ("tag" %in% names(mcols(tss_tx))) {
    tss_tx$tag_str <- vapply(mcols(tss_tx)$tag, function(x) paste0(x, collapse=","), "")
  } else {
    tss_tx$tag_str <- ""
  }
  tss_tx
}

# Helper: pick one transcript TSS for a gene (MANE > APPRIS principal > closest-to-CRE)
select_one_tss <- function(tss_tx_gene, CRE_mid) {
  if (length(tss_tx_gene) == 0) return(NULL)
  # 1) MANE Select
  mane_idx <- grepl("MANE_Select", tss_tx_gene$tag_str, fixed = TRUE)
  if (any(mane_idx)) return(tss_tx_gene[which(mane_idx)[1]])
  # 2) APPRIS principal
  appr_idx <- grepl("appris_principal", tss_tx_gene$tag_str, fixed = TRUE)
  if (any(appr_idx)) return(tss_tx_gene[which(appr_idx)[1]])
  # 3) Closest transcript TSS to the CRE midpoint
  d <- abs(start(tss_tx_gene) - start(CRE_mid)[1])
  return(tss_tx_gene[which.min(d)])
}


# FINAL: signed distance (transcript-aware)
get_signed_distance_tx <- function(gene, CRE_range, tss_tx) {
    
  CRE_gr  <- Signac::StringToGRanges(CRE_range, sep = c("-", "-"))
  CRE_mid <- resize(CRE_gr, width = 1, fix = "center")

  # transcripts for this gene (GENCODE uses 'gene_name')
  gene_hits <- which(mcols(tss_tx)$gene_name == gene)
  tss_tx_gene <- tss_tx[gene_hits]

  # same chromosome only
  tss_tx_gene <- tss_tx_gene[as.character(seqnames(tss_tx_gene)) == as.character(seqnames(CRE_mid))[1]]
  if (length(tss_tx_gene) == 0) return(NA_integer_)

  tss_one <- select_one_tss(tss_tx_gene, CRE_mid)
  if (is.null(tss_one)) return(NA_integer_)

  as.integer(start(tss_one) - start(CRE_mid)[1])
}
                             
get_nearby_genes <- function(CRE,tss_gr) {
    CRE_range <- unique(obj@meta.data[obj@meta.data$CRE_id == CRE,]$CRE)
    CRE_range <- Signac::StringToGRanges(CRE_range, sep = c("-","-"))

    ### A) Distance-based (robust near chrom ends)
    hits <- distanceToNearest(tss_gr, CRE_range, ignore.strand = TRUE)
    near_idx <- queryHits(hits)[mcols(hits)$distance <= 1e6]   # 1,000,000 bp
    tss_near <- tss_gr[near_idx]

    ### B) Window-based (±1 Mb around the CRE)
    win <- resize(CRE_range, width = width(CRE_range) + 2e6, fix = "center")
    # If your seqinfo has seqlengths, you can also do: win <- trim(win)
    tss_near_win <- subsetByOverlaps(tss_gr, win, ignore.strand = TRUE)

    genes <- unique(as.data.frame(tss_near_win)$gene_name)
    genes <- genes[genes %in% rownames(obj)]

    genes 
    }
                             
gtf_file <- "../top-100_K562-3N-3d-T7/gencode.v43.basic.annotation.gtf"
gtf <- import(gtf_file) 
                             
tss_tx_gr <- make_transcript_tss(gtf)
                             
                             
#############################
# performing the differential
# expression testing for 
# each cCRE with at least
# 100 cells
#############################                           
DefaultAssay(obj) <- "RNA"
CRE_de <- as.data.frame(matrix(ncol=0,nrow=0))

for (CRE in n_cells[n_cells$n.cells >= 100,]$CRE) {
    try({
    message(CRE)



    genes <- get_nearby_genes(CRE,tss_tx_gr)

    CRE_cells <- rownames(obj@meta.data[grepl(CRE,obj@meta.data$CRE_id) == TRUE,])

    # Get the total number of cells with a deletion for 
    # this CRE
    total.cells <- length(CRE_cells)

    # Get all cells assigned a deletion
    all_cells <- rownames(obj@meta.data[obj@meta.data$CRE_id != "",])
    all_cells <- all_cells[!(all_cells %in% CRE_cells)]

    # Get the gene that the deletion overlaps
    CRE_gene <- unique(obj@meta.data[obj@meta.data$CRE_id == CRE,]$CRE_gene)

    obj@meta.data$group <- ifelse(rownames(obj@meta.data) %in% CRE_cells, "deletion",
                           ifelse(obj@meta.data$CRE_id == "", "no deletion",""))
    Idents(obj) <- obj@meta.data$group

    de_ <- FindMarkers(obj,
                      ident.1 = "deletion",
                      ident.2 = "no deletion",
                      features = genes,
                      min.pct  = 0, 
                      logfc.threshold = 0, 
                      #test.use="t",
                      min.diff.pct    = -Inf)
    de_$CRE <- CRE
    de_$gene <- rownames(de_)
    rownames(de_) <- paste0(CRE,":",de_$gene)
    de_$CRE_gene <- CRE_gene
        
    de_$total_cells <- total.cells
        
    CRE_range <- unique(obj@meta.data[obj@meta.data$CRE_id == CRE,]$CRE)
    de_$distance <- vapply(de_$gene,
                           get_signed_distance_tx,
                           CRE_range = CRE_range,
                           tss_tx    = tss_tx_gr,
                           FUN.VALUE = integer(1))

    CRE_de <- rbind(CRE_de,de_)
    })
    }
                             
#############################
# subsetting to candidate
# enhancer-gene pairs
#############################                         
CRE_de <- CRE_de[(CRE_de$pct.1 >= .05 & CRE_de$pct.2 >= .05),]

CRE_de$type <- ifelse(CRE_de$avg_log2FC <= -.3 &
       (CRE_de$pct.1 >= .05 | CRE_de$pct.2 >= .05) &
        CRE_de$p_val < .1, "p-value < .1",
                       "p-value > .1") 
                             
CRE_de <- CRE_de[CRE_de$type == "p-value < .1", ]
                             
                             
                        
                             
#############################
# FDR
#############################                              
empirical_fdr <- function(CRE,gene)  {                           
    CRE_cells <- rownames(obj@meta.data[grepl(CRE,obj@meta.data$CRE_id) == TRUE,])

    # Get the total number of cells with a deletion for 
    # this CRE
    total.cells <- length(CRE_cells)

    # Get all cells assigned a deletion
    all_cells <- rownames(obj@meta.data[obj@meta.data$CRE_id != "",])
    all_cells <- all_cells[!(all_cells %in% CRE_cells)]

    obj@meta.data$group <- ifelse(rownames(obj@meta.data) %in% CRE_cells, "deletion",
                           ifelse(obj@meta.data$CRE_id == "", "no deletion",""))
    Idents(obj) <- obj@meta.data$group

    de_ <- FindMarkers(obj,
                      ident.1 = "deletion",
                      ident.2 = "no deletion",
                      features = gene,
                      min.pct  = 0, 
                      logfc.threshold = 0, 
                      min.diff.pct    = -Inf)
    de_$type <- "true deletion"
    de_$gene <- rownames(de_)
    de_$CRE <- CRE
    rownames(de_) <- paste0(CRE,":",gene," true deletion")

    for (i in 1:1000) {

        set.seed(i)
        message(i)

        CRE_cells <- sample(all_cells,total.cells)

        obj@meta.data$group <- ifelse(rownames(obj@meta.data) %in% CRE_cells, "deletion",
                               ifelse(obj@meta.data$CRE_id == "", "no deletion",""))
        Idents(obj) <- obj@meta.data$group

        de_sc <- FindMarkers(obj,
                          ident.1 = "deletion",
                          ident.2 = "no deletion",
                          features = gene,
                          min.pct  = 0, 
                          logfc.threshold = 0, 
                          min.diff.pct    = -Inf)
        de_sc$type <- "sim deletion"
        de_sc$gene <- rownames(de_sc)
        de_sc$CRE <- CRE
        rownames(de_sc) <- paste0(CRE,":",gene," sim deletion ",i)
        de_ <- rbind(de_,de_sc)  

    }
    de_
    }

#############################
# Finding FDR for all candidate
# enhancer-gene pairs
############################# 
fdr_out <- as.data.frame(matrix(ncol=0,nrow=0))
for (row in rownames(CRE_de)) {
    sub_de <- CRE_de[row,]
    gene <- sub_de$gene
    CRE <- sub_de$CRE
    message(paste0("Finding FDR for ",CRE, " DE of ",gene))
    fdr_out <- rbind(fdr_out,empirical_fdr(CRE,gene))
    message("")
}

saveRDS(fdr_out,file="empirical_fdr_results.wilcox.rds") 
