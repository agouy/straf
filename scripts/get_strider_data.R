# library(ggplot2)
# library(ggrepel)

fname <- "./inst/application/www/STRidER_frequencies_2024-09-24.csv"
strider_frequencies <- freq_to_mds(fname)
save(strider_frequencies, file="./inst/application/www/strider_frequencies.rda")
load("./inst/application/www/strider_frequencies.rda")


strider_frequencies <- strider_frequencies[!rownames(strider_frequencies) %in% c("SAUDI_ARABIA"), ]

X_2 <- strider_frequencies[ 10:11,]
idx_2 <- colSums(is.na(X_2)) == 0

X <- strider_frequencies[-10:-11,]
idx <- colSums(is.na(X)) == 0
idx

common_names <- intersect(names(idx[idx]), names(idx_2[idx_2]))
loci_names <- unique(do.call(rbind, strsplit(common_names, "_"))[, 1])
print(loci_names)

d <- X %*% t(X)
vec <- sqrt(diag(d))
d <- d/vec[col(d)]
d <- d/vec[row(d)]
d <- -log(d)
d <- as.dist(d)

mds <- cmdscale(d)
MDS <- data.frame(ax1 = mds[, 1], ax2 = mds[, 2], pop = rownames(mds))

p <- ggplot(MDS, aes(x=ax1, y=ax2, color = pop, label = pop)) +
  geom_point() +
  geom_text_repel() + 
  labs( x = "MDS Axis 1", y = "MDS Axis 2", title = "MDS based on Nei's distance")  +
  theme_minimal()
plot(p)
