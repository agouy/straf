#' Convert input file to genind object
#' @param Ifile Input file object.
#' @param Imicrovariants Number of microvariants.
#' @param Incode Number of digits encoding allele size.
#' @param Iploidy Ploidy ("Haploid" or "Diploid").
#' @return An object of class genind.
#' @export
#' @noRd
#' @importFrom adegenet as.genind df2genind pop<- locNames
#' @importClassesFrom adegenet genind
createGenind <- function(Ifile, Imicrovariants, Incode, Iploidy) {
  
  if(Imicrovariants == 2) {
    
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
    
    if(Iploidy == "Haploid") {
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
    }
    
    dat2 <- df2genind(dat[, -1:-2],ncode = switch(
      Incode, "2" = 2, "3" = 3), ploidy = switch(
        Iploidy, Diploid = 2, Haploid = 1))
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
