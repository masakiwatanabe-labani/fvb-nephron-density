## E13.5 vs P1 comparison of the interferon response using PRE-SPECIFIED MSigDB gene sets.
## The 18-gene panel used for Figure 7B was assembled from the leading edge of the P1 result and is
## therefore post-hoc. Here the comparison is repeated on HALLMARK_INTERFERON_GAMMA_RESPONSE and
## HALLMARK_INTERFERON_ALPHA_RESPONSE taken whole, which are defined independently of these data.
suppressPackageStartupMessages({library(msigdbr); library(data.table)})
a <- commandArgs(TRUE); g <- function(k,d=NULL){i<-which(a==paste0("--",k)); if(length(i)) a[i+1] else d}
deseq <- g("deseq-dir"); id2 <- g("id2name"); out <- g("out")
m <- msigdbr(species="Mus musculus", collection="H")
sets <- split(m$gene_symbol, m$gs_name)
keep <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE")
map <- fread(id2, header=FALSE, col.names=c("gene_id","sym")); map <- map[!duplicated(gene_id)]
res <- list()
for (tp in c("E13.5","P1")) {
  d <- fread(file.path(deseq, paste0("DESeq2_", tp, ".tsv")))
  d <- merge(d, map, by="gene_id", all.x=TRUE)
  setorder(d, sym, gene_id); d <- d[!duplicated(sym) & !is.na(sym) & sym != ""]
  for (s in keep) {
    x <- d[sym %in% sets[[s]]]
    res[[length(res)+1]] <- data.table(set=s, timepoint=tp, n_in_set=length(unique(sets[[s]])),
      n_tested=nrow(x), median_log2FC=median(x$log2FoldChange, na.rm=TRUE),
      median_baseMean=median(x$baseMean, na.rm=TRUE),
      n_padj_lt_0.1=sum(x$padj < 0.1, na.rm=TRUE),
      frac_padj_lt_0.1=round(sum(x$padj < 0.1, na.rm=TRUE)/sum(!is.na(x$padj)), 3))
  }
}
r <- rbindlist(res); fwrite(r, out, sep="\t"); print(r)
