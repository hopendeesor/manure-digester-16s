nc <- readRDS("data/processed/seqtab_nochim.rds")
taxa <- readRDS("data/processed/taxa.rds")
rel <- 100 * nc / rowSums(nc)
lab <- function(x) ifelse(is.na(x), "unclassified", x)
ph <- lab(taxa[, "Phylum"])
ge <- lab(taxa[, "Genus"])
top <- function(x, tax, n) {
  v <- sort(tapply(x, tax, sum), decreasing=TRUE)[1:n]
  paste0(names(v), " ", round(v, 1), "%", collapse="; ")
}
show <- function(name, x) {
  cat("\n", name, "\n")
  cat("  phyla: ", top(x, ph, 4), "\n")
  cat("  genera:", top(x, ge, 5), "\n")
}
for (s in c("s12", "c7", "s13", "s14")) show(s, rel[s, ])
show("MEAN OF ALL 20", colMeans(rel))
