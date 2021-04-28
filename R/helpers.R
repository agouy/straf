### straf helpers
# EXTERNAL FUNCTIONS -----------------------------------------------------------

### convert input to genind object
createGenind <- function(Ifile, Imicrovariants, Incode, Iploidy) {
  
  if(Imicrovariants == 2) {
    
    mat <- readLines(Ifile$datapath)
    mat <- strsplit(mat, "[\t]")
    
    mat <- matrix(unlist(mat), nrow = length(mat), ncol = length(mat[[1]]),
                  byrow = TRUE)
    mat[mat == "0"] <- NA ###
    colnames(mat) <- mat[1,]
    rownames(mat) <- mat[,1]
    
    mat <- mat[-1,]
    loci <- unique(colnames(mat[,-1:-2]))
    freqTAB <- NULL
    mat2 <- mat
    mat2 <- sub("[.]","-",mat2)
    
    for(i in 1:length(loci)){
      ids <- which(colnames(mat)==loci[i])
      alleles <- unique(c(mat[,ids]))
      alleles <- sub("[.]","-",alleles)
      alleles <- alleles[!is.na(alleles)] ###
      nameCol <- paste(loci[i],".",alleles,sep = "")
      
      newmat <- matrix(NA,ncol = length(nameCol),nrow = dim(mat)[1])
      for(ii in 1:length(alleles)){
        newmat[,ii] <- apply(mat2[,ids]==alleles[ii],1,sum)
        colnames(newmat) <- nameCol
      }
      freqTAB <- cbind(freqTAB,newmat)
    }
    rownames(freqTAB) <- mat[,1]
    colnames(freqTAB) <- sub(" ","",colnames(freqTAB))
    dat2 <- as.genind(tab = freqTAB)
    pop(dat2) <- mat[,"pop"]
    
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

## getFreqAllPop
## returns allele frequencies for each population
# input: genind object
# output: a list of allele frequencies tables
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

## getFreqFromGenind
## returns allele frequencies for a given population
# input: genind object
# output: an allele frequencies tables
getFreqFromGenind <- function(data) {
  
  freq <- apply(data@tab, 2, sum, na.rm = TRUE)
  nam <- strsplit(
    names(freq),
    split = "[.]"
  )
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


## getIndicesAllPop
## returns indices for each population
# input: genind object
# output: a list of indices tables
getIndicesAllPop <- function(data, hw = FALSE, hwperm = 1000, ploidy = "Diploid") {
  ind <- list()
  ind$all <- getIndicesFromGenind(data, hw, hwperm, ploidy)
  for(popu in unique(data$pop)) {
    ind <- c(ind, x = NA)
    mat <- getIndicesFromGenind(data[data@pop == popu, ], hw, hwperm, ploidy)
    ind$x <- mat
    names(ind)[length(ind)] <- popu
  }
  return(ind)
}

## getIndicesFromGenind
## returns indices for a given population
# input: genind object
# output: an indices table
getIndicesFromGenind <- function(data,
                                 hw = FALSE,
                                 hwperm = 1000,
                                 ploidy = "Diploid") {
  
  freq <- apply(data@tab, 2, sum, na.rm = TRUE)
  
  nam <- strsplit(
    names(freq),
    split="[.]"
  )
  
  loc <- as.factor(unlist(
    lapply(nam, function(x) x[1])
  ))
  
  alle <- as.numeric(unlist(
    lapply(nam, function(x) sub("-", ".", x[2]))
  ))
  DAT <- data.frame(freq, loc, alle)
  N <- tapply(DAT$freq, DAT$loc, sum)
  DAT$frequency <- DAT$freq / N[DAT$loc]

  PIC <- NULL
  for(i in unique(loc)) {
    
    #FR <- c(DAT$frequency[names(DAT$frequency) == i])
    FR <- c(DAT$frequency[DAT$loc == i])
    
    xu <- outer(FR, FR, "fu")
    som <- sum(xu[lower.tri(xu)])
    PIC[i] <-  1 - sum(FR ^ 2) - som
    
  } 
  Nall <- tapply(
    DAT[DAT$freq>0, ]$freq,
    DAT[DAT$freq>0, ]$loc,
    length
  )
  
  GD <- tapply(
    DAT$frequency,
    DAT$loc,
    function(x) 1 - sum(x ^ 2)
  )
  
  GD <- GD * N / (N - 1) 
  PIC <- PIC[names(GD)]
  D2 <- genind2loci(data)
  sumloc <- summary(D2)[names(GD)]
  PM1 <- lapply(sumloc, function(x) {
    sum((x$genotype / sum(x$genotype)) ^ 2)
  })
  
  PM <- unlist(PM1)
  
  DF <- data.frame(
    locus = names(GD),
    N = N,
    Nall = Nall,
    GD = GD,
    PIC = PIC,
    PM = PM,
    PD = 1 - PM
  )

  if(ploidy == "Diploid") {
    DF$Hobs <- adegenet::summary(data)$Hobs[names(GD)]
    DF$PE <- (DF$Hobs ^ 2) * (1 - 2 * (DF$Hobs) * ((1 - DF$Hobs) ^ 2))
    DF$TPI <- 1 / (2 * (1 - DF$Hobs))
  }
  
  
  if(length(unique(data@pop)) > 1 & length(locNames(data)) > 1) {
    basicstat <- basic.stats(
      data,
      diploid = switch(ploidy, Diploid = TRUE, Haploid = FALSE),
      digits = 4
    )$perloc
    rownames(basicstat) <- as.character(unique(data@loc.fac))
    Fst <- wc(
      data,
      diploid = switch(ploidy, Diploid = TRUE, Haploid = FALSE)
    )$per.loc$FST
    names(Fst) <- as.character(unique(data@loc.fac))

    DF$Fst <- Fst[names(GD)]
    # DF$Fst <- basicstat[names(GD), "Fst"]
    DF$Ht <- basicstat[names(GD), "Ht"]
    DF$Fis <- basicstat[names(GD), "Fis"]
    
  }
  
  if(ploidy == "Diploid" & hw) {
    withProgress(message = 'Performing HW test...', value = 0, {
      DF$pHW <- hw.test(data, B = hwperm)[names(GD), 4]
    })
  } 
  
  return(DF)
}


fu <- function(a, b){
  2 * (a ^ 2) * (b ^ 2)
}


plotPCA <- function(pca, popus, coul, axis) {
  
  var1 <- round(100*(pca$eig/sum(pca$eig))[axis[1]], 2)
  var2 <- round(100*(pca$eig/sum(pca$eig))[axis[2]], 2)
  
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
    ellipse(
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

straf2genepop <- function(f.name, ploidy = 2) {
  df <- readLines(f.name)
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  df_tmp <- lapply(df[, -1:-2], function(x) gsub("[.]", "", x))
  
  # add leading zeros
  df_tmp2 <- lapply(df_tmp, function(x) {
    x[nchar(x) == 1] <- paste0(x[nchar(x) == 1], "00")
    x[nchar(x) == 2] <- paste0(x[nchar(x) == 2], "0")
    if(any(nchar(x) != 3)) stop("Error while converting allele labels.")
    return(x)
  })
  
  # concatenate
  idx <- seq_len(length(df_tmp2))
  
  if(ploidy == 2) {
    ids <- idx %% 2
    
    df_out <- list()
    for(i in idx[as.logical(ids)]) {
      nm <- names(df_tmp2[i])
      nm <- gsub(" ", "", nm)
      df_out[[nm]] <- paste0(df_tmp2[[i]], df_tmp2[[i + 1]])
    }
  } else if (ploidy == 1) {
    df_out <- list()
    for(i in idx) {
      nm <- names(df_tmp2[i])
      nm <- gsub(" ", "", nm)
      df_out[[nm]] <- df_tmp2[[i]]
    }
  }
  
  df_out <- as.data.frame(df_out)
  
  first.line <- "STRAF-generated GENEPOP input file."
  loci <- colnames(df_out)
  
  str_out <- apply(df_out, 1, paste0, collapse = "\t")
  
  ## get pops
  populations <- unique(df$pop)
  vec_out <- c()
  for(i in populations) {
    idx <- df$pop %in% i
    vec_out <- c(
      vec_out,
      "Pop",
      paste(df[idx, ]$ind, str_out[idx], sep = "\t,\t")
    )
  }
  
  ## write file
  output <- c(first.line, loci, vec_out)
  output <- paste(output, "\n", collapse = "")
  
  return(output)
}


straf2familias <- function(f.name) {
  
  df <- readLines(f.name)
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  df_tmp <- df[, -1:-2]
  
  # add leading zeros
  # concatenate
  idx <- seq_len(length(df_tmp))
  ids <- as.logical(idx %% 2)
  
  df_out <- list()
  for(i in idx[ids]) {
    nm <- names(df_tmp[i])
    nm <- gsub(" ", "", nm)
    df_out[[nm]] <- c(df_tmp[[i]], df_tmp[[i + 1]])
  }
  
  
  tbs <- lapply(df_out, table)
  prop.tb <- lapply(tbs, prop.table)
  
  str.list <- lapply(prop.tb, function(x) {
    vec <- x
    str_loc <- paste0(names(vec), "\t", unname(vec), collapse = "\n")
    return(str_loc)
  })
  
  output <- paste(names(prop.tb), str.list, sep = "\n")
  out <- paste(output, collapse = "\n\n")
  out <- paste0(out, "\n")
  return(out)
}

straf2arlequin <- function(f.name) {
  df <- readLines(f.name)
  
  spt <- do.call("rbind", strsplit(df, "\t"))
  colnames(spt) <- spt[1, ]
  df <- as.data.frame(spt[-1, ])
  
  df_tmp2 <- list()
  
  for(i in seq_len(ncol(df[, -1:-2]))) {
    x <- df[, -1:-2][, i]
    if(i %% 2 != 0) {
      x2 <- df[, -1:-2][, i + 1]
      dot_idx <- grep("[.]", x)
      dot_idx2 <- grep("[.]", x2)
      if(length(dot_idx) > 0 | length(dot_idx2) > 0) {
        x <- gsub("[.]", "", x)
        x[-dot_idx] <- paste0(x[-dot_idx], "0")
        x2 <- gsub("[.]", "", x2)
        x2[-dot_idx2] <- paste0(x2[-dot_idx2], "0")
        if(length(dot_idx) == 0) x <- paste0(x, "0")
        if(length(dot_idx2) == 0) x2 <- paste0(x2, "0")
      }
      df_tmp2[[i]] <- x
      df_tmp2[[i+1]] <- x2
      
    } 
    
  }
  df_tmp2 <- as.data.frame(df_tmp2)
  
  # concatenate
  idx <- seq_len(length(df_tmp2))
  ids <- as.logical(idx %% 2)
  
  out_str <- c()
  for(pp in unique(df$pop)) {
    df_pop <- df[df$pop == pp, ]
    out_str <- c(out_str, paste0('SampleName="', pp, '"\nSampleSize=',nrow(df_pop),'\nSampleData={\n'))
    
    for(i in which(df$pop == pp)) {
      samp_nm <- df[i, 1]
      l1 <- paste0(c(samp_nm, "1", unname(unlist(df_tmp2[i, c(idx[ids])]))), collapse = "\t")
      l2 <- paste0(c("", "", unname(unlist(df_tmp2[i, c(idx[!ids])]))), collapse = "\t")
      out_str <- c(out_str, l1, l2)
    }
    out_str <- c(out_str, "}\n\n")
  }
  npop <- length(unique(df$pop))
  header <- paste0('[Profile]\nTitle="STRAF-generated Arlequin file."\nNbSamples=',npop,'\nDataType=MICROSAT\n
GenotypicData=1\nGameticPhase=0\nMissingData="?"\nLocusSeparator=WHITESPACE\n\n[Data]\n\n[[Samples]]\n\n')
  
  output <- c(header, out_str)
  ## write file
  output <- paste(output, "\n", collapse = "")
  
  return(output)
}

freq_to_mds <- function(fname) {
  ln <- readLines(fname)
  ln2 <- lapply(ln, function(x) strsplit(x, ",")[[1]])
  ln3 <- lapply(ln2, function(x) {
    if(sum(nchar(x[-1]) == 0) == length(x[-1])) return(x[1])
    else return(x)
  })
  hd <- lengths(ln3)
  names_idx <- which(hd == 1)
  st_idx <- names_idx + 1
  en_idx <- names_idx - 1
  en_idx <- c(en_idx[-1], length(ln2))
  df <- lapply(seq_along(names_idx), function(i) {
    loc_id <- names_idx[i]
    loc_name <- ln3[[loc_id]]
    if(en_idx[i] - st_idx[i] > 1) {
      mat <- do.call(rbind, ln2[st_idx[i]:en_idx[i]])
      colnames(mat) <- mat[1, ]
      mat[mat == ""] <- "0"
      df <- as.data.frame(mat[-1:-2, ])
      colnames(df) <- gsub(pattern = " ", replacement = "_", colnames(df))
      colnames(df) <- gsub(pattern = "\"", replacement = "", colnames(df))
      
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
  tt <- reshape2::acast(df_l, location ~ locus + Allele, value.var = 'frequency', fun.aggregate = mean, fill = -1)
  ct <- rownames(tt)
  tt <- tt %>% as_tibble()
  df_f <- tt %>% as_tibble() %>% mutate_all(~ifelse(.x == -1, NA, .x)) #mean(.x[.x != -1], na.rm = TRUE)
  matt <- (as.matrix(df_f))
  rownames(matt) <- ct
  return(matt)
}
