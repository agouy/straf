#' Convert input file to genind object
#' @param Ifile Input file object.
#' @param ploidy Ploidy ("Haploid" or "Diploid").
#' @return An object of class genind.
#' @export
#' @noRd
createGenind <- function(Ifile, ploidy) {
  
  if(ploidy == 2) {
    
    mat <- readLines(Ifile$datapath)
    mat <- strsplit(mat, "[\t]")
    
    mat <- matrix(unlist(mat), nrow = length(mat), ncol = length(mat[[1]]),
                  byrow = TRUE)
    mat[mat == "0"] <- NA
    colnames(mat) <- mat[1,]
    rownames(mat) <- mat[,1]
    
    mat <- mat[-1,]
    loci <- unique(colnames(mat[,-1:-2]))
    freqTAB <- NULL
    mat2 <- mat
    mat2 <- sub("[.]", "-", mat2)
    
    for(i in 1:length(loci)){
      ids <- which(colnames(mat)==loci[i])
      alleles <- unique(c(mat[,ids]))
      alleles <- sub("[.]","-",alleles)
      alleles <- alleles[!is.na(alleles)]
      nameCol <- paste(loci[i],".",alleles,sep = "")
      
      newmat <- matrix(NA, ncol = length(nameCol), nrow = dim(mat)[1])
      for(ii in 1:length(alleles)){
        newmat[, ii] <- apply(mat2[, ids] == alleles[ii], 1, sum)
        colnames(newmat) <- nameCol
      }
      freqTAB <- cbind(freqTAB,newmat)
    }
    rownames(freqTAB) <- mat[, 1]
    colnames(freqTAB) <- sub(" ", "", colnames(freqTAB))
    dat2 <- as.genind(tab = freqTAB)
    pop(dat2) <- mat[, "pop"]
    
  } else {
    dat <- read.table(Ifile$datapath, header = TRUE,
                      sep = "\t", colClasses = "character")
    rownames(dat) <- dat$ind
    dat_tmp <- dat[, -1:-2]
    if(length(grep("[.]", unlist(dat_tmp))) > 0) {
      new_dat <- apply(dat_tmp, MARGIN = 2, function(x) {
        x <- gsub("[.]", "", x)
        x[nchar(x) == 1] <- paste0("0", x[nchar(x) == 1], "0")
        x[nchar(x) == 2] <- paste0(x[nchar(x) == 2], "0")
        if(any(nchar(x) != 3)) stop("Allele encoding error.")
        return(x)
      })
      dat[, -1:-2] <- new_dat
    }
    
    dat2 <- df2genind(dat[, -1:-2], ncode = 3, ploidy = ploidy)
    pop(dat2) <- dat$pop
  }
  return(dat2)
}

#' Get allele frequencies for each population.
#' @param data a genind object
#' @return an allele frequency table
#' @export
#' @noRd
getFreqAllPop <- function(data) {
  freq <- list()
  freq$all <- getFreqFromGenind(data)
  for(popu in unique(data$pop)) {
    freq <- c(freq, x = NA)
    mat <- getFreqFromGenind(data[data@pop == popu, ])
    freq$x <- mat
    names(freq)[length(freq)] <- popu
  }
  return(freq)
}

#' Get allele frequencies for a given population.
#' @param data a genind object
#' @return an allele frequency table
#' @export
#' @noRd
getFreqFromGenind <- function(data) {
  freq <- apply(data@tab, 2, sum, na.rm = TRUE)
  nam <- strsplit(names(freq), split = "[.]")
  loc <- as.factor(unlist(
    lapply(nam, function(x) x[1])
  ))
  alle <- as.numeric(unlist(
    lapply(nam, function(x) sub("-", ".", x[2]))
  ))
  
  DAT <- data.frame(freq, loc, alle)
  N <- tapply(DAT$freq, DAT$loc, sum)
  DAT$frequency <- DAT$freq / N[DAT$loc]
  
  DAT$alle[is.na(DAT$alle)] <- 0
  matr <- matrix(
    NA,
    ncol = length(unique(DAT$loc)),
    nrow = length(unique(DAT$alle))
  )
  rownames(matr) <- sort(unique(DAT$alle))
  colnames(matr) <- sort(unique(DAT$loc))
  
  for(i in sort(unique(DAT$alle))) {
    matr[
      as.character(i),
      DAT[DAT$alle == i, ]$loc
    ] <- DAT[DAT$alle==i, ]$frequency
  }
  
  return(matr)
}

#' Plot PCA.
#' @export
#' @noRd
plotPCA <- function(pca, popus, coul, axis) {
  
  var1 <- round(100 * (pca$eig / sum(pca$eig))[axis[1]], 2)
  var2 <- round(100 * (pca$eig / sum(pca$eig))[axis[2]], 2)
  
  plot(pca$li[, axis[1]],
       pca$li[, axis[2]],
       col = transp(coul[popus], 0.5),
       pch = 16, cex = 1.7,
       xlab = paste0("PC", axis[1], " (", var1, " %)"),
       ylab = paste0("PC", axis[2], " (", var2, " %)"),
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5,
       bty = "l",
       main = "PCA projection"
  )
  
  sapply(unique(popus), function(x) {
    car::ellipse(
      c(mean(pca$li[popus %in% x, axis[1]]), mean(pca$li[popus %in% x, axis[2]])),
      cov(pca$li[, axis[1:2]]),
      0.95,
      col = transp(coul[x], 1),
      lty = 1,
      lwd = 2,
      center.pch = 0,
      fill = TRUE,
      fill.alpha = 0.2,
      segments = 100
    )
  })
  
  invisible(sapply(unique(popus), function(x) {
    legend(
      mean(pca$li[popus %in% x, axis[1]]),
      mean(pca$li[popus %in% x, axis[2]]),
      x,
      xjust = 0.5,
      yjust = 0.5,
      text.col=coul[x],
      text.font = 2,
      box.col = NA,
      bg = transp("white", 0.5),
      adj = 0.2
    )
  }))
  
}


#' Convert a POPTREE file to a format suitable for use in STRAF
#' @param fname Input file path
#' @return An object of class genind.
#' @export
#' @noRd
poptree2straf <- function(fname) {
  # df <- readxl::read_excel(fname, sheet = 1, skip = 0, col_names = FALSE)
  df <- readlines(fname)
  
  popnames <- df[[1]]
  popnames <- popnames[2:(which(is.na(popnames))[1] - 1)]
  n_pops <- length(popnames)
  
  a <- grep("locus", df[[1]])
  b <- grep("#", df[[1]])
  if(length(a) != length(b)) stop("problem")
  str_out <- c()
  for(i in seq_along(a)) {
    
    df_tmp <- df[a[i]:(b[i]-1), ]
    pop_names <- df_tmp[1, 4:(4+n_pops-1)]
    Locus <- df_tmp[1, 2][[1]]
    freq <- df_tmp[, 4:(4+n_pops-1)]
    all_names <- c("Allele", unname(unlist(df_tmp[-1, 2])))
    df_tmp2 <- cbind(all_names, freq)
    df_tmp3 <- apply(df_tmp2, 1, paste0, collapse = ",")
    str_out <- c(str_out, Locus, df_tmp3)
  }
  str_out <- gsub(",[.]", ",0." , str_out)
  return(str_out)
}


#' Pairwise Linkage Disequilibrium computation.
#' @param data a genind object
#' @return a table containing Genepop's LD p-values
#' @export
#' @noRd
getGenepopLD <- function(data, ploidy, n_iter) {
  
  gp_tmp <- tempfile()
  on.exit(unlink(gp_tmp))
  
  out_tmp <- tempfile()
  on.exit(unlink(out_tmp))
  
  straf2genepop(data, gp_tmp, ploidy)
  genepop::test_LD(
    gp_tmp, 
    out_tmp, 
    dememorization = 5000, 
    batches = 100, 
    iterations = n_iter
  )
  
  df <- get_ld_gp(out_tmp)
  
  lnam <- sort(unique(c(as.character(df$Locus_2), as.character(df$Locus_1))))
  df$Locus_1 <- factor(df$Locus_1, levels=rev(lnam))
  df$Locus_2 <- factor(df$Locus_2, levels=rev(lnam))
  
  df$plab <- sub('^(-)?0[.]', '\\1.', round(df$P_value, 3))
  df$plab[df$plab == "0"] <- "< .0001"
  par(mar = c(4,4,4,8))
  df <- dplyr::distinct(df)
  return(df)
}

#' Plot Pairwise Linkage Disequilibrium.
#' @param df a table containing LD test results
#' @param pop string, population to represent
#' @return a heatmap representing LD test p-values
#' @export
#' @noRd
plot_ld <- function(df, pop) {
  df <- df[df$Population == df$Population[1], ]
  
  plt <- ggplot(data = df, aes(Locus_1, Locus_2, fill = -log10(P_value))) +
    geom_tile(color = "white", na.value='lightgrey') +
    scale_fill_gradient2(low = "#fee6ce", mid = "#fdae6b", high = "#e6550d", 
                         midpoint = 2.5, limit = c(-0.1,5), space = "Lab", 
                         name="LD\nSignificance\n\n-log10(p-value)\n") +
    theme_minimal() + 
    scale_x_discrete(position = "top") +
    theme(
      axis.text.x = element_text(angle = 90, size = 9), axis.text.y = element_text(angle = 0, size = 9),
      panel.grid.major = element_blank(),
      legend.position=c(1.2, 0.6),
      plot.caption=element_text(vjust = 30, hjust = 1.5, size = 10),
      plot.margin = margin(t = 0.2, r =  4, b = 0.2, l = 0.2, "cm")
    ) + 
    geom_text(aes(label = plab), size = 2, col = "black") +
    ggtitle(label = "Linkage disequilibrium between loci") +
    xlab("") + ylab("") +
    coord_fixed()
  
  plotly::ggplotly(plt)
  
}

#' Parse Linkage Disequilibrium Genepop results.
#' @param out.name path to genepop file
#' @return a dataframe containing LD test p-values
#' @export
#' @noRd
get_ld_gp <- function(out.name) {
  tx <- readLines(out.name)
  genepop::clean_workdir()
  
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
  df$P_value <- df$P_value + 1/10000
  df$P_value[df$P_value > 1] <- 1
  
  df$fdr <- p.adjust(df$P_value, method = "bonferroni")
  
  df$lab <- NA
  df$lab[df$fdr < 0.05] <- "*"
  df$lab[df$fdr < 0.01] <- "**"
  df$lab[df$fdr < 0.001] <- "***"
  tmp <- colnames(df)
  tmp[2:3] <- c("Locus_2", "Locus_1")
  df_tmp <- df[,tmp]
  colnames(df_tmp) <- colnames(df)
  df <- rbind(df, df_tmp)
  
  return(df)
}


#' Hardy-Weinberg Equilibrium test..
#' @param data path to straf input file
#' @param ploidy ploidy
#' @return a table containing Genepop's HW test p-values
#' @export
#' @noRd
getGenepopHW <- function(data, ploidy, n_iter = 10000) {
  gp_tmp <- tempfile()
  on.exit(unlink(gp_tmp))
  
  out_tmp <- tempfile()
  on.exit(unlink(out_tmp))
  
  straf2genepop(data, gp_tmp, ploidy)
  genepop::test_HW(
    gp_tmp,
    which = "Proba", 
    out_tmp,
    batches = 100, 
    iterations = n_iter 
  )
  
  df <- get_hw_gp(out_tmp)
  return(df)
}

#' Parse HWE test Genepop results.
#' @param fname path to genepop file
#' @return a dataframe containing HWE test p-values
#' @export
#' @noRd
get_hw_gp <- function(fname) {
  tx <- readLines(fname)
  genepop::clean_workdir()
  
  spl <- strsplit(tx, "\\s+")
  
  pop_idx <- grep("Pop", spl)
  pop_names <- unlist(lapply(spl[pop_idx], "[", 3))
  
  all_idx <- grep("Fisher's", spl)
  all_idx <- all_idx[all_idx < pop_idx[1]]
  all_pvals <- unlist(lapply(spl[all_idx + 3], "[", 4))
  
  ## ln = 7 and between pop_idx, add pop_names
  for(i in seq_along(pop_idx)) {
    if(i == length(pop_idx)) idd <- pop_idx[i]:length(spl)
    else idd <- pop_idx[i]:pop_idx[i+1]
    for(ii in idd) spl[[ii]] <- c(pop_names[i], spl[[ii]])  
  }

  spl <- spl[-seq_len(pop_idx[1])]
  ln <- lengths(spl)

  
  tb <- do.call(rbind, spl[ln == 8])
  
  df <- as.data.frame(tb)

  colnames(df) <- c("Population", "Locus", "P_value", "Std_Error", 
                    "Fis_WC", "Fis_RH", "Steps", "Switched")
  # df <- tidyr::complete(df, Locus, tidyr::nesting(Population))
  
  # df$Population <- rep(pop_names, each = nrow(tb) / length(pop_names))
  
  df <- dplyr::select(df, Population, Locus, P_value)
  
  df$P_value <- as.numeric(df$P_value)
  df$P_value <- df$P_value + 1/10000
  df$P_value[df$P_value > 1] <- 1
  
  df_2 <- data.frame(Population = "all", Locus = unique(df$Locus), P_value = all_pvals)
  df <- rbind(df, df_2)
  return(df)
  
}
