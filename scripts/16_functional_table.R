nc <- readRDS("data/processed/seqtab_nochim.rds")
taxa <- readRDS("data/processed/taxa.rds")
md <- read.delim("data/processed/digester_metadata.tsv")
res <- read.delim("results/14_diffabund.tsv")
rd <- read.table("data/genus_roles.txt", sep="|", quote="", comment.char="", stringsAsFactors=FALSE)
names(rd) <- c("genus", "role", "source")
stopifnot(identical(rownames(nc), md$sample))
g <- unname(taxa[, "Genus"])
g[is.na(g)] <- "no_genus"
gt <- t(rowsum(t(nc), group=g))
rel <- 100 * gt / rowSums(gt)
top <- head(res$genus, 12)
all <- unique(c(top, rd$genus))
cs <- colSums(gt > 0)
mc <- colMeans(rel[md$group == "cattle", , drop=FALSE])
ms <- colMeans(rel[md$group == "swine", , drop=FALSE])
out <- data.frame(genus=all, stringsAsFactors=FALSE)
out$set <- ifelse(all %in% top, "top-12 differential", "reference genus")
out$n_samples <- unname(cs[all])
out$pct_cattle <- round(unname(mc[all]), 3)
out$pct_swine <- round(unname(ms[all]), 3)
out$q <- round(res$q[match(all, res$genus)], 4)
out$role <- rd$role[match(all, rd$genus)]
out$source <- rd$source[match(all, rd$genus)]
none <- is.na(out$role)
out$role[none] <- "no role annotated (no source verified)"
out$source[none] <- "-"
write.table(out, "results/16_functional_table.tsv", sep="\t", quote=FALSE, row.names=FALSE)
o2 <- lapply(out, function(x) ifelse(is.na(x), "n/a", as.character(x)))
hdr <- "| Genus | Set | Samples (of 20) | Mean % cattle | Mean % swine | q | Role from literature | Sources |"
sep <- "|---|---|---|---|---|---|---|---|"
body <- paste0("| ", do.call(paste, c(o2, sep=" | ")), " |")
writeLines(c(hdr, sep, body), "results/16_functional_table.md")
cat("Genera in table:", nrow(out), "\n")
cat("Not found in our data:", paste(out$genus[is.na(out$n_samples)], collapse=", "), "\n")
cat("No role annotated:", sum(none), "\n")
