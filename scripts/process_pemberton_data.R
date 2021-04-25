library(dplyr)
library(tidyr)
d <- read.table(
  "C:/Users/alexa/Desktop/pembertonEtAl2013.MS5795_formatted.stru",
  header=TRUE
)

sc <- seq_len(nrow(d)) %% 2 == 0
info <- d[sc, c(1:5)]
loc1 <- d[sc, -c(1:5)]
loc2 <- d[!sc, -c(1:5)]
df2 <- cbind(loc1, loc2)
df2 <- cbind(info, df2[, order(colnames(df2))])
ind <- paste(df2$popid, df2$id, sep = "_")
df2 <- cbind(ind, df2)
df3 <- df2[, -which(colnames(df2) %in% c("id", "popid", "subpop", "superpop"))]
df3[df3 == "-9"] <- 0

df3 <- rbind(colnames(df3), df3)
df3[1, ] <- gsub("[.]1", "", df3[1, ])

test <- read.table("./scripts/UniSTS.aliases", sep = "\t")

std_loc <- c("CSF1PO", "FGA", "TH01", "TPOX", " VWA", "D3S1358", "D5S818", "D7S820", "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11", " D1S1656", "D2S441", "D2S1338", "D10S1248", "D12S391", "D19S433", "D22S1045", "FGA", "TH01", "VWA", "D1S1656", "D2S441", "D3S1358", "D8S1179", " D10S1248", "D12S391", "D18S51", "D21S11", "D22S1045", "D2S1338", "D16S539", "D19S433", "SE33", "Amelogenin", "FGA", "TH01", "VWA", "D2S1338", " D3S1358", " D8S1179", "D16S539", "D18S51", "D19S433", "D21S11", "FGA", "TH01", "VWA", "D3S1358", "D8S1179", "D18S51", "D21S11")
std_loc <- unique(trimws(std_loc))

req <- paste0(std_loc, collapse = "|")
idx <- grep(req, test[,2])

strr <- paste(test[idx, 2], collapse = ";")

cidx <- colnames(df3)[grep(pattern = gsub(";", "|", strr), colnames(df3))]
df4 <- df3[, c("ind", "pop", cidx)]
df4[1,] <- gsub("TPO[.]D2S", "TPOX", df4[1,])

df4[df4 == "0"] <- NA

df4 <- df4 %>% tidyr::drop_na()
tb <- table(df4$pop)

df4 <- df4[c(1, which(df4$pop %in% names(tb[tb >= 10]))), ]
df4 <- df4[which(!df4$ind %in% "1114_71587"), ]
 
saveRDS(df4, file = "./scripts/pemberton2013.rds")

write.table(df4, sep = "\t",
            quote=FALSE, row.names = FALSE, col.names = FALSE,
            file="./scripts/pemberton2013_straf.txt")



source("./scripts/helpers.R")
Ifile <- list()
Ifile$datapath <- "./scripts/pemberton2013_straf.txt"
library(adegenet)
dat2 <- createGenind(Ifile, 2, 2, "Diploid")
obj <- genind2genpop(dat2, quiet = FALSE)
cn <- colnames(obj@tab)
cn <- gsub("[.]", "_", cn)
cn <- gsub("[-]", ".", cn)
colnames(obj@tab) <- cn

mt <- as.matrix(obj@tab)
saveRDS(mt, "./www/pemberton_frequencies.rds")
