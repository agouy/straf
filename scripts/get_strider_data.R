library(tidyr)
fname <- "./STRidER_frequencies_2019-08-02.csv"

ln <- readLines(fname)
ln2 <- lapply(ln, function(x) strsplit(x, ",")[[1]])
hd <- lengths(ln2)
names_idx <- which(hd == 1)
st_idx <- names_idx + 1
en_idx <- names_idx - 1
en_idx <- c(en_idx[-1], length(ln2))

df <- lapply(seq_along(names_idx), function(i) {
  loc_id <- names_idx[i]
  loc_name <- ln2[[loc_id]]
  
  if(en_idx[i] - st_idx[i] > 1) {
    mat <- do.call(rbind, ln2[st_idx[i]:en_idx[i]])
    colnames(mat) <- mat[1, ]
    mat[mat == ""] <- "0"
    df <- as.data.frame(mat[-1:-2, ])
    colnames(df) <- gsub(pattern = " ", replacement = "_", colnames(df))
    colnames(df) <- gsub(pattern = "\"", replacement = "", colnames(df))
    
    print(i)
    print(loc_id)
    print(df)
    
    df_long <- gather(df, location, frequency, -Allele, factor_key=TRUE)
    df_long$locus <- loc_name
    return(df_long)
    
  } else {
    return(NULL)
  }
})
df_l <- do.call(rbind, df)

df_l$frequency <- as.numeric(df_l$frequency)
df_l$location  <- as.character(df_l$location)

library(reshape2)
tt <- reshape2::acast(df_l, location ~ locus + Allele, value.var = 'frequency', fun.aggregate = mean, fill = -1)

ct <- rownames(tt)
tt <- tt %>% as_tibble()

library(dplyr)
df_f <- tt %>% as_tibble() %>% mutate_all(~ifelse(.x == -1, mean(.x[.x != -1], na.rm = TRUE), .x))          


df_f
matt <- (as.matrix(df_f))
rownames(matt) <- ct

strider_frequencies <- matt
save(strider_frequencies, file="./www/strider_frequencies.rda")


library(ggplot2)
library(ggrepel)
load("./www/strider_frequencies.rda")

X <- strider_frequencies
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
