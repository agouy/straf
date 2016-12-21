# SERVER -----------------------------------------------------------------------

library(shiny);library(ade4);library(adegenet)
library(pegas);library(hierfstat);library(DT)

shinyServer(function(input, output) {

# GET INPUT FILE -------------------------------------------------------------
    
  getData <- reactive({
    if(is.null(input$file1)) return(NULL)
    inFile <- input$file1
  })
  
  output$fileUploaded <- reactive({
    return(!is.null(getData()))
  })
  
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  
# DISPLAY DATASET --------------------------------------------------------------
  
  output$contents <- renderDataTable({
    if (is.null(getData()) | !input$displayTable)
      return(NULL)
    X<-read.table(getData()$datapath, header = TRUE,
                  sep = "\t",colClasses="character")
    DT::datatable(X)
  })
  
# ALLELE FREQUENCIES -----------------------------------------------------------

  output$alleleFreq <- renderPlot({
    if (is.null(getData()) | !input$displayAlleleFreq)
      return(NULL)
    
    dat2<-createGenind(Ifile=getData(),Imicrovariants=input$microvariants,
                       Incode=input$ncode,Iploidy=input$ploidy)
    
    freq<-apply(dat2@tab,2,sum,na.rm=TRUE)
    
    
    nam<-strsplit(names(freq),split="[.]")
    loc<-as.factor(unlist(lapply(nam,function(x)x[1])))
    alle<-as.numeric(unlist(lapply(nam,function(x)sub("-",".",x[2]))))
    DAT<-data.frame(freq,loc,alle)
    DAT<-DAT[order(DAT$alle),]
    
    if(length(unique(DAT$loc))<=5) ncolo<-2
    if(length(unique(DAT$loc))>=6) ncolo<-3
    if(length(unique(DAT$loc))>=10) ncolo<-4
    par(mfrow=c(ceiling(length(unique(DAT$loc))/ncolo),ncolo))
    par(mar=c(2,2,2,2))
    for(i in unique(DAT$loc)){
      barplot(DAT$freq[DAT$loc==i],names.arg=DAT$alle[DAT$loc==i],
              main=i,col=transp(input$barplotcolor,input$transparency),
              border=as.numeric(input$borderbarplot))
    }
  }
  )
  
  output$plotAF <- renderUI({
    plotOutput('alleleFreq', width=paste(input$width,"%",sep=""), height=input$height)
  })
  
# ALLELE FREQUENCIES TABLE -----------------------------------------------------
  
  output$tableFreq <- renderDataTable({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayAlleleTable)
      return(NULL)
    
    dat2<-createGenind(Ifile=inFile,Imicrovariants=input$microvariants,
                       Incode=input$ncode,Iploidy=input$ploidy)
    
    freq<-apply(dat2@tab,2,sum,na.rm=TRUE)
    nam<-strsplit(names(freq),split="[.]")
    loc<-as.factor(unlist(lapply(nam,function(x)x[1])))
    alle<-as.numeric(unlist(lapply(nam,function(x)sub("-",".",x[2]))))
    
    DAT<-data.frame(freq,loc,alle)
    N<-tapply(DAT$freq,DAT$loc,sum)
    DAT$frequency<-DAT$freq/N[DAT$loc]
    
    DAT$alle[is.na(DAT$alle)]<-0
    matr<-matrix(NA,ncol=length(unique(DAT$loc)),nrow=length(unique(DAT$alle)))
    rownames(matr)<-sort(unique(DAT$alle))
    colnames(matr)<-sort(unique(DAT$loc))
    
    for(i in sort(unique(DAT$alle))){
      matr[as.character(i),DAT[DAT$alle==i,]$loc]<-DAT[DAT$alle==i,]$frequency  
    }
    DT::datatable(matr) %>% formatRound(columns=colnames(matr),digits=3)
  }
  )
  
  tabfreqInput <- reactive({
    inFile <- input$file1
      
    if (is.null(inFile) | !input$displayAlleleTable)
      return(NULL)
    
    dat2<-createGenind(Ifile=inFile,Imicrovariants=input$microvariants,
                       Incode=input$ncode,Iploidy=input$ploidy)
    
    freq<-apply(dat2@tab,2,sum,na.rm=TRUE)
    nam<-strsplit(names(freq),split="[.]")
    loc<-as.factor(unlist(lapply(nam,function(x)x[1])))
    alle<-as.numeric(unlist(lapply(nam,function(x)x[2])))
    
    DAT<-data.frame(freq,loc,alle)
    N<-tapply(DAT$freq,DAT$loc,sum)
    DAT$frequency<-DAT$freq/N[DAT$loc]
    
    DAT$alle[is.na(DAT$alle)]<-0
    matr<-matrix(NA,ncol=length(unique(DAT$loc)),nrow=length(unique(DAT$alle)))
    rownames(matr)<-sort(unique(DAT$alle))
    colnames(matr)<-sort(unique(DAT$loc))
    
    for(i in sort(unique(DAT$alle))){
      matr[as.character(i),DAT[DAT$alle==i,]$loc]<-DAT[DAT$alle==i,]$frequency  
    }
    tabF<-cbind(allele=rownames(matr),matr)
  }
  )
  
  output$dlTabfreq <- downloadHandler(
    filename = function() { 
      paste('allele_frequencies.tsv', sep='') 
    },
    content = function(file) {
      write.table(tabfreqInput(), file, sep="\t", na = "", row.names = FALSE)
    }
  )
  
# STATISTICS TABLE -------------------------------------------------------------
  
  reacIndices <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    dat2<-genind4LD(Ifile=inFile,Imicrovariants=input$microvariants,
                    Incode=input$ncode,Iploidy=input$ploidy)
    
    freq<-apply(dat2@tab,2,sum,na.rm=TRUE)
    nam<-strsplit(names(freq),split="[.]")
    loc<-as.factor(unlist(lapply(nam,function(x)x[1])))
    alle<-as.numeric(unlist(lapply(nam,function(x)sub("-",".",x[2]))))
    
    DAT<-data.frame(freq,loc,alle)
    N<-tapply(DAT$freq,DAT$loc,sum)
    DAT$frequency<-DAT$freq/N[DAT$loc]
    
    
    ##### maj
    PIC<-NULL
    for(i in unique(loc)){
      FR<-c(DAT$frequency[names(DAT$frequency)==i])
      som<-0
      for(ii in 1:(length(FR)-1)){
        for(iii in (ii+1):(length(FR))){
          som<-som+2*(FR[ii]^2)*(FR[iii]^2)
        }
      }
      PIC[i]<- 1-sum(FR^2)-som
    } 
    
    D2<-genind2loci(dat2)
    sumloc<-summary(D2)
    PM1<-lapply(sumloc,function(x){
      sum((x$genotype/sum(x$genotype))^2)
    })
    PM<-unlist(PM1)
    PD<-1-PM
    ##### maj
    
    
    Nall<-tapply(DAT$freq,DAT$loc,length)
    GD<-tapply(DAT$frequency,DAT$loc,function(x)1-sum(x^2))
    GD<-GD*N/(N-1) ##### maj
    
    # PIC<-tapply(DAT$frequency,DAT$loc,function(x){1-sum(x^2)-(sum(x^2))^2+sum(x^4)})
    
    if(input$ploidy=="Diploid") {
      Hobs<-summary(dat2)$Hobs[names(GD)]
      PE<-(Hobs^2)*(1-2*(Hobs)*((1-Hobs)^2))
      TPI<-1/(2*(1-Hobs))
    }
    
    
    
    Nall<-tapply(DAT$freq,DAT$loc,length)
#     GD<-tapply(DAT$frequency,DAT$loc,function(x)1-sum(x^2))
#     GD<-GD*N/(N-1) ##### maj, Nei
#     PIC<-tapply(DAT$frequency,DAT$loc,function(x){1-sum(x^2)-(sum(x^2))^2+sum(x^4)})
#     PD<-tapply(DAT$frequency,DAT$loc,function(x){1-2*((sum(x^2))^2)-sum(x^4)})
#     PE<-(GD^2)*(1-(1-GD)*(GD^2))
#     TPI<-1/(2*(1-GD))
    
    if(length(unique(dat2@pop)) > 1){
    basicstat <- basic.stats(dat2, diploid = switch(
      input$ploidy, Diploid = TRUE, Haploid = FALSE), digits = 4)
    #add fstats
    
    Fst<-basicstat$perloc[names(GD),"Fst"]
    
      
    Ht<-basicstat$perloc[names(GD),"Ht"]
    Fis<-basicstat$perloc[names(GD),"Fis"]
    if(input$ploidy=="Diploid") {
      if(input$computeHW) {
          pHW<-hw.test(dat2,B=1000)[names(GD),4]
          DF<-data.frame(locus=names(GD),N,Nall,Hobs,GD,Ht,PIC,PM,PD,PE,TPI,Fst,Fis,pHW);DF
      } else {
        DF<-data.frame(locus=names(GD),N,Nall,Hobs,GD,Ht,PIC,PM,PD,PE,TPI,Fst,Fis);DF
      }
    } else {
      DF<-data.frame(locus=names(GD),N,Nall,GD,Ht,PIC,PM,PD,Fst);DF
    }
      }else{
        if(input$ploidy=="Diploid") {
          if(input$computeHW) {
            pHW<-hw.test(dat2,B=1000)[names(GD),4]
            DF<-data.frame(locus=names(GD),N,Nall,Hobs,GD,PIC,PM,PD,PE,TPI,pHW);DF
          } else {
            DF<-data.frame(locus=names(GD),N,Nall,Hobs,GD,PIC,PM,PD,PE,TPI);DF
          }
        } else {
          DF<-data.frame(locus=names(GD),N,Nall,GD,Ht,PIC,PM,PD);DF
        }
      }
      indicesNames <- colnames(DF)
      DF
  })
  
  ### forensics display + download

  output$forensics <- renderTable({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    dat2<-genind4LD(Ifile=inFile,Imicrovariants=input$microvariants,
                    Incode=input$ncode,Iploidy=input$ploidy)
    if(length(unique(dat2@pop))>1) reacIndices()[,-which(colnames(reacIndices())%in%c("Ht","Fis","Fst"))]
    else reacIndices()
  }, digits = 4)
  
  output$dlForensics <- downloadHandler(
    filename = function() { 
      paste('forensics_parameters.tsv', sep='') 
    },
    content = function(file) {
      write.table(reacIndices()[,-which(colnames(reacIndices())%in%c("Ht","Fis","Fst"))], 
                  file, sep="\t", row.names = FALSE)
    }
  )
  
  ### popgen display + download
  
  output$diversity <- renderTable({
    reacIndices()[,-which(colnames(reacIndices())%in%c("PIC","PD","PE","TPI"))]
  }, digits = 4)
  
  output$dlPopgen <- downloadHandler(
    filename = function() { 
      paste('popgen_statistics.tsv', sep='') 
    },
    content = function(file) {
      write.table(reacIndices()[,-which(colnames(reacIndices())%in%c("PIC","PD","PE","TPI"))], 
                  file, sep="\t",row.names = FALSE)
    }
  )
  
  ### plots
  
  output$uiFOR <- renderUI({
    selectInput("plotIndicesFOR", "Select the statistic to plot:",
                choices = colnames(reacIndices())[-1:-2],
                selected = "Nall"
    )
  })
  
  output$uiPG <- renderUI({
    selectInput("plotIndicesPG", "Select the statistic to plot:",
                choices = colnames(reacIndices())[-1:-2],
                selected = "Nall"
    )
  })
  
  output$plotIndices<- renderPlot({
    if(is.null(reacIndices())) return(NULL)
    dat<-reacIndices()
    barplot(dat[,input$plotIndicesFOR],names.arg = dat[,"locus"],horiz=TRUE,
            las=1,col=transp(input$barplotcolor,input$transparency),
            border=as.numeric(input$borderbarplot))
  })
  output$plotFOR <- renderUI({
    plotOutput('plotIndices', width=paste(input$width,"%",sep=""), height=input$height)
  })
  
  output$plotIndicesPopgen<- renderPlot({
    if(is.null(reacIndices())) return(NULL)
    datpl<-reacIndices()
    barplot(datpl[,input$plotIndicesPG],names.arg = datpl[,"locus"],horiz=TRUE,
            las=1,col=transp(input$barplotcolor,input$transparency),
            border=as.numeric(input$borderbarplot))
  })
  output$plotPG <- renderUI({
    plotOutput('plotIndicesPopgen', width=paste(input$width,"%",sep=""), height=input$height)
  })
  
  
  
  
  
# PAIRWISE FST -----------------------------------------------------------------
  
  fstMatInput <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayFstMat)
      return(NULL)
    
    dat2<-createGenind(Ifile=inFile,Imicrovariants=input$microvariants,
                       Incode=input$ncode,Iploidy=input$ploidy)
    
    if (length(unique(dat2@pop)) == 1)
      stop("Multiple populations are required to perform this analysis")
    matFST<-pairwise.fst(dat2, res.type = "matrix")
    matFST[lower.tri(matFST)]<-NA
    matFST
  })
  
  output$FstMat <- renderTable({
    matFST<-apply(fstMatInput(),1,rev)
  }, digits = 4,include.rownames =TRUE,na="")
  
  output$dlFstMat <- downloadHandler(
    filename = function() { 
      paste('pairwise_fst.tsv', sep='') 
    },
    content = function(file) {
      matFST<-apply(fstMatInput(),1,rev)
      write.table(matFST, file, sep="\t", na = "",row.names = TRUE)
    }
  )
  
# LD ---------------------------------------------------------------------------
  
  reacLDtable <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    D<-genind4LD(Ifile=inFile,Imicrovariants=input$microvariants,
                 Incode=input$ncode,Iploidy=input$ploidy)
    
    datLD<-genind2loci(D)
    loci<-(unique(D@loc.fac))
    nloc<-length(unique(D@loc.fac))
    LDmat<-matrix(NA,nrow=nloc,ncol=nloc)
    
    npairs<-nloc*(nloc-1)/2
    withProgress(message = 'LD computation', value = 0, {
      for(i in 1:nloc){
        for(ii in 1:nloc){
          if(i<ii){
            if(input$ploidy=="Diploid") lx<-LD2(datLD,
                                                locus=c(loci[i],loci[ii]))$T2
            if(input$ploidy=="Haploid") lx<-LD(datLD,
                                               locus=c(loci[i],loci[ii]))$T2
            LDmat[i,ii]<-lx[3]
            LDmat[ii,i]<-lx[1]
            incProgress(1/npairs,message="Computing LD...")
          } 
        }
      }
    })
    
    rownames(LDmat)<-colnames(LDmat)<-loci
    LDmat
  })
  
  output$LDtable <- renderTable({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayLDtable)
      return(NULL)
    
    M<-reacLDtable()
    M[lower.tri(M)] <- NA
    M<-apply(M,1,rev)
  }, digits = 4,include.rownames =TRUE,na="")
  
  output$plotLD <- renderPlot({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayLDplot)
      return(NULL)
    
    M<- -log10(reacLDtable())
    M[lower.tri(M)] <- NA
    col<-redpal(100)
    image(M,col = col,frame=F,xaxt="n",yaxt="n")
    axis(2,at=seq(0,1,length.out=ncol(M)),labels=colnames(M),las=2,cex.axis=0.8)
    axis(3,at=seq(0,1,length.out=nrow(M)),labels=rownames(M),las=2,cex.axis=0.8)
  })
  output$plotLD2 <- renderUI({
    plotOutput('plotLD', width=paste(input$width,"%",sep=""), height=input$height)
  })
  
  output$plotLDpval <- renderPlot({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayLDplot)
      return(NULL)
    
    M<- (reacLDtable())
    par(mfrow=c(1,2))
    pv<-M[upper.tri(M)]
    hist(pv,xlab="P-value",main = "LD p-values distribution")
    qqplot(pv,qunif(seq(0,1,0.01)),pch=16,xlab="Observed quantiles",
           ylab="Expected quantiles",
           main=paste("KS test p-value:",
                      round(ks.test(pv,qunif(seq(0,1,0.01)))$p.value,
                            digits = 4)))
    abline(0,1,lty=2,lwd=2,col="grey")
  })
  
  output$plotLDpval2 <- renderUI({
    plotOutput('plotLDpval', width=paste(input$width,"%",sep=""), height=input$height)
  })
  
  output$dlLDtable <- downloadHandler(
    filename = function() { 
      paste('LD_pvalues.tsv', sep='') 
    },
    content = function(file) {
      pairLD<-reacLDtable();pairLD[lower.tri(pairLD)] <- NA
      pairLD<-apply(pairLD,1,rev)
      write.table(pairLD, file, sep="\t", na = "",row.names = TRUE)
    }
  )
  
# PCA --------------------------------------------------------------------------
  
  output$runPCA <- renderPlot({
    inFile <- input$file1
    
    if (is.null(inFile) | !input$displayPCA)
      return(NULL)
    
    dat2<-createGenind(Ifile=inFile,Imicrovariants=input$microvariants,
                       Incode=input$ncode,Iploidy=input$ploidy)
    
    FREQ<-makefreq(dat2,missing="mean")
    pcaY<-dudi.pca(FREQ, nf = 3, scannf = FALSE)
    
    
    if(length(input$PCAaxis)==2){
      par(mfrow=c(1,1))
      s.class(pcaY$li, fac=pop(dat2),
              col=transp(funky(12),.6),
              axesel=FALSE, cstar=0, cpoint=2,xax=as.numeric(input$PCAaxis[1]),
              yax=as.numeric(input$PCAaxis[2]))
    }
    if(length(input$PCAaxis)==3){
      
      par(mfrow=c(1,2))
      s.class(pcaY$li, fac=pop(dat2),
              col=transp(funky(12),.6),
              axesel=FALSE, cstar=0, cpoint=2,xax=as.numeric(input$PCAaxis[1]),
              yax=as.numeric(input$PCAaxis[2]))
      s.class(pcaY$li, fac=pop(dat2),
              col=transp(funky(12),.6),
              axesel=FALSE, cstar=0, cpoint=2,xax=as.numeric(input$PCAaxis[1]),
              yax=as.numeric(input$PCAaxis[3]))
    }
  })
  output$plotPCA <- renderUI({
    plotOutput('runPCA', width=paste(input$width,"%",sep=""), height=input$height)
  })
  
  
  output$loadings <- renderPlot({
    inFile <- input$file1
    if(is.null(inFile) | !input$displayPCA) return(NULL)
    
    dat2<-createGenind(Ifile=inFile,Imicrovariants=input$microvariants,
                       Incode=input$ncode,Iploidy=input$ploidy)
    FREQ<-makefreq(dat2,missing="mean")
    pcaY<-dudi.pca(FREQ, nf = 3, scannf = FALSE)
    loadingplot(pcaY$c1^2)
  })
  output$plotLoadings <- renderUI({
    plotOutput('loadings', width=paste(input$width,"%",sep=""), height=input$height)
  })
})

# EXTERNAL FUNCTIONS -----------------------------------------------------------

createGenind <- function(Ifile,Imicrovariants,Incode,Iploidy) {
  if(Imicrovariants==2){
    readLines(Ifile$datapath)->matrix
    matrix<-strsplit(matrix,"[\t]")
    
    mat<-matrix(unlist(matrix),nrow=length(matrix),ncol=length(matrix[[1]]),
                byrow=TRUE)
    mat[mat=="0"]<-NA ###
    colnames(mat)<-mat[1,]
    rownames(mat)<-mat[,1]
    
    mat<-mat[-1,]
    loci<-unique(colnames(mat[,-1:-2]))
    freqTAB<-NULL
    mat2<-mat
    mat2<-sub("[.]","-",mat2)
    for(i in 1:length(loci)){
      ids<-which(colnames(mat)==loci[i])
      alleles<-unique(c(mat[,ids]))
      alleles<-sub("[.]","-",alleles)
      alleles<-alleles[!is.na(alleles)] ###
      nameCol<-paste(loci[i],".",alleles,sep="")
      
      newmat<-matrix(NA,ncol=length(nameCol),nrow=dim(mat)[1])
      for(ii in 1:length(alleles)){
        newmat[,ii]<-apply(mat2[,ids]==alleles[ii],1,sum)
        colnames(newmat)<-nameCol
      }
      freqTAB<-cbind(freqTAB,newmat)
    }
    rownames(freqTAB)<-mat[,1]
    colnames(freqTAB)<-sub(" ","",colnames(freqTAB))
    dat2<-as.genind(tab=freqTAB)
    pop(dat2)<-mat[,"pop"]
    
  } else {
    dat<-read.table(Ifile$datapath, header = TRUE,
                    sep = "\t",colClasses="character")
    rownames(dat)<-dat$ind
    dat2<-df2genind(dat[,-1:-2],ncode=switch(
      Incode, "2" = 2, "3" = 3),ploidy=switch(
        Iploidy, Diploid = 2, Haploid = 1));pop(dat2)<-dat$pop
  }
  return(dat2)
}

genind4LD<-function(Ifile,Imicrovariants,Incode,Iploidy){
  if(Imicrovariants==2){
    readLines(Ifile$datapath)->matrix
    matrix<-strsplit(matrix,"[\t]")

  mat<-matrix(unlist(matrix),nrow=length(matrix),ncol=length(matrix[[1]]),
              byrow=TRUE)
  mat[mat=="0"]<-NA ###
  colnames(mat)<-mat[1,]
  rownames(mat)<-mat[,1]
  
  mat<-mat[-1,]
  mat<-sub("[.]","",mat)
#   mat[nchar(mat)==2]<-paste(mat[nchar(mat)==2],"0",sep="")
#   mat[nchar(mat)==1]<-paste(mat[nchar(mat)==1],"00",sep="")
  mat[nchar(mat)==2&!is.na(mat)]<-paste(mat[nchar(mat)==2&!is.na(mat)],"0",sep="") ###
  mat[nchar(mat)==1&!is.na(mat)]<-paste(mat[nchar(mat)==1&!is.na(mat)],"00",sep="") ###
  loci<-unique(colnames(mat[,-1:-2]))
  freqTAB<-NULL
  for(i in 1:length(loci)){
    ids<-which(colnames(mat)==loci[i])
    alleles<-unique(c(mat[,ids]))
    alleles<-alleles[!is.na(alleles)] ###
    
    nameCol<-paste(loci[i],".",alleles,sep="")
    
    newmat<-matrix(NA,ncol=length(nameCol),nrow=dim(mat)[1])
    for(ii in 1:length(alleles)){
      newmat[,ii]<-apply(mat[,ids]==alleles[ii],1,sum)
      colnames(newmat)<-nameCol
    }
    freqTAB<-cbind(freqTAB,newmat)
  }
  rownames(freqTAB)<-mat[,1]
  colnames(freqTAB)<-sub(" ","",colnames(freqTAB))
  
  D<-genind(tab=freqTAB,pop=mat[,"pop"])
  
  } else {
    dat<-read.table(Ifile$datapath, header = TRUE,
                    sep = "\t",colClasses="character")
    rownames(dat)<-dat$ind
    D<-df2genind(dat[,-1:-2],ncode=switch(
      Incode, "2" = 2, "3" = 3),ploidy=switch(
        Iploidy, Diploid = 2, Haploid = 1));pop(D)<-dat$pop
  }
  return(D)
}