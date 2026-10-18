library(vegan)
library(ggplot2)
set.seed(1)
nc <- readRDS("data/processed/seqtab_nochim.rds")
md <- read.delim("data/processed/digester_metadata.tsv")
stopifnot(identical(rownames(nc), md$sample))
depth <- min(rowSums(nc))
rr <- rrarefy(nc, depth)
a <- data.frame(sample=md$sample, group=md$group)
a$observed <- rowSums(rr > 0)
a$shannon <- diversity(rr, index="shannon")
a$fa <- md$fa
a$ph <- md$ph
a$nh4n <- md$nh4n
write.table(a, "results/11_alpha.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("Rarefied to", depth, "reads per sample\n")
print(aggregate(cbind(observed, shannon) ~ group, a, median))
cat("Wilcoxon p, observed:", wilcox.test(observed ~ group, a, exact=FALSE)$p.value, "\n")
cat("Wilcoxon p, shannon:", wilcox.test(shannon ~ group, a, exact=FALSE)$p.value, "\n")
for (v in c("fa", "ph", "nh4n")) {
  for (m in c("observed", "shannon")) {
    r <- cor.test(a[[v]], a[[m]], method="spearman", exact=FALSE)
    cat(m, "vs", v, "rho =", round(r$estimate, 2), "p =", round(r$p.value, 4), "\n")
  }
}
l1 <- data.frame(group=a$group, metric="Observed ASVs", value=a$observed)
l2 <- data.frame(group=a$group, metric="Shannon", value=a$shannon)
l <- rbind(l1, l2)
p <- ggplot(l, aes(group, value)) + geom_boxplot(outlier.shape=NA)
p <- p + geom_jitter(aes(colour=group), width=0.1, height=0)
p <- p + facet_wrap(~metric, scales="free_y") + theme_bw()
ggsave("figures/11_alpha_by_group.png", p, width=7, height=4, dpi=150)
