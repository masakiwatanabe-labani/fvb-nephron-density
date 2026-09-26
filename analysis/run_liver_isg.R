#!/usr/bin/env Rscript
# Is the FVB/N interferon-low signature kidney-specific, or a strain-wide property?
# Public replication set: GSE167328 vehicle-control livers (C57BL/6J n=3, FVB/NJ n=3, pregnant females,
# 6-8 wk). That study also aligned to strain-specific pseudogenomes (UNC CCstatus + Lapels), so it
# carries the same class of reference-bias correction as ours.
suppressPackageStartupMessages({library(DESeq2); library(fgsea); library(msigdbr); library(data.table)})
a <- commandArgs(trailingOnly = TRUE); g <- function(k, d=NULL){i<-which(a==paste0("--",k)); if(length(i)) a[i+1] else d}
counts_f <- g("counts"); map_f <- g("id2name"); outdir <- g("outdir"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

m <- as.matrix(read.delim(counts_f, row.names = 1, check.names = FALSE))
rownames(m) <- sub("\\..*$", "", rownames(m))
m <- m[rowSums(m) >= 10, ]
cond <- factor(ifelse(grepl("^FVB", colnames(m)), "FVB", "B6"), levels = c("B6", "FVB"))
dds <- DESeqDataSetFromMatrix(m, DataFrame(strain = cond), ~ strain)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, alpha = 0.1))
res$gene_id <- rownames(res)

map <- fread(map_f, header = FALSE, col.names = c("gid", "sym"))
map[, gid := sub("\\..*$", "", gid)]
map <- map[!duplicated(gid)]
res$sym <- map$sym[match(res$gene_id, map$gid)]
write.csv(res, file.path(outdir, "liver_FVB_vs_B6_DESeq2.csv"), row.names = FALSE)
cat(sprintf("genes tested: %d   padj<0.1: %d\n", nrow(res), sum(res$padj < 0.1, na.rm = TRUE)))

r <- res[!is.na(res$stat) & !is.na(res$sym) & res$sym != "", ]
r <- r[order(r$sym, r$gene_id), ]; r <- r[!duplicated(r$sym), ]
stats <- setNames(r$stat, r$sym)

sets <- list()
for (cc in list(c("H", NA), c("C2", "CP:REACTOME"))) {
  mm <- if (is.na(cc[2])) msigdbr(species="Mus musculus", collection=cc[1])
        else msigdbr(species="Mus musculus", collection=cc[1], subcollection=cc[2])
  sets <- c(sets, split(mm$gene_symbol, mm$gs_name))
}
sets <- sets[!duplicated(names(sets))]
keep <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "REACTOME_INTERFERON_SIGNALING","REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
          "HALLMARK_MYC_TARGETS_V1","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_FATTY_ACID_METABOLISM")
set.seed(42)
fg <- fgsea(sets, stats, minSize = 15, maxSize = 500, eps = 0, sampleSize = 1001)
fg <- fg[pathway %in% keep][order(padj)]
fg[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
fwrite(fg, file.path(outdir, "liver_gsea_targeted.tsv"), sep = "\t")
cat("\n--- targeted sets, liver, FVB vs B6 (n=3 each) ---\n")
print(fg[, .(pathway, NES, pval, padj, size)])

core <- c("Ifit3","Gbp3","Irf7","Oas2","Oas3","Oasl1","Rnasel","Bst2","Xaf1","Ifih1",
          "Usp18","Mx2","Isg15","Rtp4","Ddx60","Stat1","Stat2","Irf9")
cg <- res[res$sym %in% core, c("sym","baseMean","log2FoldChange","padj")]
cg <- cg[order(cg$log2FoldChange), ]
cat("\n--- core ISGs in liver ---\n"); print(cg, row.names = FALSE)
cat(sprintf("\nliver core-ISG median log2FC = %+.2f   padj<0.1: %d/%d\n",
            median(cg$log2FoldChange, na.rm=TRUE), sum(cg$padj < 0.1, na.rm=TRUE), nrow(cg)))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("\nLIVER_ISG_DONE\n")
