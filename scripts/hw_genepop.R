# install.packages("genepop")

library(genepop)
locinfile <- ('./tests/input_files/exampleSTRAFdiplo.txt')
run_gp_ld
run_gp_ld <- function(input_file) {
  
}

gp_tmp <- tempfile()
on.exit(unlink(gp_tmp))
gp_tmp_out <- tempfile()
on.exit(unlink(gp_tmp_out))

cat(straf2genepop(locinfile), file = gp_tmp)

test_LD(gp_tmp, gp_tmp_out, dememorization = 10000, batches = 100, iterations = 5000)

tx <- readLines(gp_tmp_out)
clean_workdir()

tx <- gsub("Not possible", "NA NA NA", tx)
tx <- gsub("No information", "NA NA NA", tx)

tx <- gsub("& ", "", tx)
# tx <- gsub("Penta E", "PentaE", tx)
spl <- strsplit(tx, "\\s+")
ln <- lengths(spl)
tb <- do.call(rbind, spl[ln == 6])

df <- as.data.frame(tb[-c(1:5), ])
colnames(df) <- c("Population", "Locus_1", "Locus_2", "P_value", "Std_Error", "Switches")

df$P_value <- as.numeric(df$P_value)
df$Std_Error <- as.numeric(df$Std_Error)

df$Locus_1 <- factor(df$Locus_1, levels = unique(df$Locus_1))
df$Locus_2 <- factor(df$Locus_2, levels = unique(df$Locus_2))


library(dplyr)
library(tidyr)

df$fdr <- p.adjust(df$P_value, method = "bonferroni") + 1e-10

df$lab <- NA
df$lab[df$fdr < 0.05] <- "*"
df$lab[df$fdr < 0.01] <- "**"
df$lab[df$fdr < 0.001] <- "***"

library(ggplot2)
f.plot <- "C:/Users/alexa/Desktop/LD_plot_v2.png"
df2 <- df[df$Population=="3k",][1:20,]
plt <- ggplot(data = df2, aes(Locus_1, Locus_2, fill = -log10(P_value))) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#fee6ce", mid = "#fdae6b", high = "#e6550d", 
                       midpoint = 2.5, limit = c(-0.1,5), space = "Lab", 
                       name="Linkage Disequilibrium\nSignificance\n\n-log10(p-value)\n") +
  theme_minimal()+ 
  scale_x_discrete(position = "top") +
  theme(
    axis.text.x = element_text(angle = 90, size = 9), axis.text.y = element_text(angle = 0, size = 9),
    panel.grid.major = element_blank(),
    legend.position=c(0.85, 0.4),
    plot.caption=element_text(vjust = 30, hjust = 0.95, size = 10)
  ) + 
  geom_text(aes(label = sub('^(-)?0[.]', '\\1.', round(P_value, 3))), size = 2) +
  ggtitle("Linkage disequilibrium between loci") +
  xlab("") + ylab("")+
  coord_fixed() +  geom_text(aes(label=lab), size=3, col = ifelse(df2$fdr < 0.001,"white", "black")) +
  labs("\u03b1 = 0.05 / 496 = 0.0001")

plt

# df
# hist(df$fdr)

ggsave(f.plot, plot = plt, device = "png", width = 20, height = 20, units = "cm", dpi = 600)

df_out <- df %>% select(Locus_1, Locus_2, P_value, Std_Error )
df_out$P_value <- round(df_out$P_value, 5)
df_out$Std_Error  <- round(df_out$Std_Error , 5)
# df_out[df_out == 0] <- "< 0.00001"

df_out_wide <- df_out
df_out_wide$val <- paste0(df_out_wide$P_value, " (+/-", df_out_wide$Std_Error , ")")
df_out_wide <- df_out_wide %>% select(Locus_1, Locus_2, val) %>% spread(Locus_2, val, fill = "-")
df_out_wide <- df_out_wide[, c("Locus_1", rev(colnames(df_out_wide)[-1]))]

library("xlsx")
file <- "C:/Users/alexa/Desktop/LD_results_genepop.xlsx"
write.xlsx(df_out, file, sheetName = "LD_analysis", 
             col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx(df_out_wide, file, sheetName = "LD_analysis_alternative_table", 
           col.names = TRUE, row.names = FALSE, append = TRUE)
