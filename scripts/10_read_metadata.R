f <- "data/raw/meta/DataSheet1.xlsx"
x <- readxl::read_excel(f, "Table S1", col_names=FALSE, .name_repair="minimal")
d <- as.data.frame(x[3:22, ])
n1 <- c("sample","province","date","digester","volume_m3","substrate")
n2 <- c("srt_day","biogas","temp_c","ph","cod","nh4n","fa","phosphate","cn")
names(d) <- c(n1, n2)
for (v in n2) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
d$volume_m3 <- suppressWarnings(as.numeric(d$volume_m3))
meta <- read.delim("data/processed/metadata.tsv")
stopifnot(setequal(d$sample, meta$sample))
d <- d[match(meta$sample, d$sample), ]
d$group <- meta$group
write.table(d, "data/processed/digester_metadata.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("Substrate code vs our group:\n")
print(table(d$substrate, d$group))
cat("\nMissing values per column:\n")
print(colSums(is.na(d)))
cat("\n")
print(d[, c("sample","province","digester","temp_c","ph","fa","nh4n","cn")])
