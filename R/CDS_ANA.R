codonList <- c(
  "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG",
  "AAT", "AAC", "GAT", "GAC", "TGT", "TGC", "CAA", "CAG", "GAA", "GAG",
  "GGT", "GGC", "GGA", "GGG", "CAT", "CAC", "ATT", "ATC", "ATA", "TTA",
  "TTG", "CTT", "CTC", "CTA", "CTG", "AAA", "AAG", "ATG", "TTT", "TTC",
  "CCT", "CCC", "CCA", "CCG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
  "ACT", "ACC", "ACA", "ACG", "TGG", "TAT", "TAC", "GTT", "GTC", "GTA",
  "GTG", "TAA", "TGA", "TAG")
CodonPercN <- "CornCodonTable_V1N.txt"

# String Find: pos <- instr(fnm, "_") + 1
instr <- function(str1, str2, startpos=1, n=1){
  aa=unlist(strsplit(substring(str1, startpos), str2))
  if(length(aa) < n+1 ) return(0);
  return(sum(nchar(aa[1:n])) + startpos + (n-1)*nchar(str2) )
}#instr

filestem <- function(filenm){
  len <- nchar(filenm)
  str1 <- substr(filenm, len-3, len)
  if(str1 == ".txt") return(substr(filenm, 1, len-4))
}#filestem

#' Create a file with protein sequence in the third column
#'
#' @param infile "DNA sequence file", "Codon Usage Table file", "Result Excel file and "Input file serial number"
#' @return None (Result file with "*_AddAAS.txt" and Excel sheet "AAseq" will be created)
#' @export
TransProtein <- function(filenm, CodonPerc, xlsxFS, f){
  fnm <- filenm
  D1 <- read.table(fnm, header = TRUE, sep = "\t")
  D1N <- nrow(D1)
  D2 <- read.table(CodonPerc, header = TRUE, sep = "\t")
  D1$AAS <- "M"
  s <- 1
  for(s in c(1:D1N)){
    str1 <- toupper(D1[s, 2])
    str2 <- stri_sub(str1, seq(1, stri_length(str1),by=3), length=3)
    strN <- length(str2)
    RLT <- ""
    t <- 1
    for(t in c(1:strN)){
      str3 <- str2[t]
      value <- D2[D2[[4]]==str3, 1]
      RLT <- paste(RLT, value, sep = "")
    }#t
    D1[s, 3] <- RLT
  }#s
  fnm0 <- filestem(filenm)
  fnm <- paste(fnm0, "_AddAAS.txt", sep = "")
  write.table(D1, file = fnm, quote = FALSE, sep = "\t", row.names = FALSE)
  sheetS <- paste(f, "-AAseq", sep = "")
  if(f == 1){
    write.xlsx(D1, xlsxFS, sheetName = sheetS, col.names = TRUE, row.names = FALSE, append = FALSE)
  }else{
    write.xlsx(D1, xlsxFS, sheetName = sheetS, col.names = TRUE, row.names = FALSE, append = TRUE)
  }
}#TransProtein

#' Create a file with Codon count numbers and a file with codon usage
#'
#' @param infile "DNA sequence file", "Codon Usage Table file", "Result Excel file and "Input file serial number"
#' @return None (Result files with "*_codonSEQ.txt" and "_codonUSTB" and Excel sheets will be created)
#' @export
CodonCount <- function(filenm, CodonPerc, xlsxFS, f) {
  fnm <- filenm
  DD <- read.table(fnm, header = TRUE, sep = "\t")
  N = nrow(DD)
  M <- length(codonList)
  colnames(DD)[1] <- "Gene_ID"
  DD$RanVal <- runif(N)*10
  for(c in c(1:M)){
    DD[c + 3] <- rep(0, N)
    colnames(DD)[c + 3] <- codonList[c]
  }
  DD <- DD[, colnames(DD)[c(1, 3, c(4:(M+3)), 2)]]
  M <- ncol(DD)
  i <- 1
  for(i in c(1:N)){
    x <- toupper(toString(DD[i, M]))
    strVec <- substring(x, seq(1, nchar(x), 3), seq(3, nchar(x), 3))
    ctTab <- as.data.frame(table(strVec))
    for(j in c(3:(M-1))){
      nm <- colnames(DD)[j]
      val <- ctTab[ctTab$strVec==nm, 2]
      if(length(val) > 0) DD[i, j] <- val
    } #for j
    DD[i,]
  } #for i
  fnm0 <- filestem(filenm)
  fnm <- paste(fnm0, "_codonSEQ.txt", sep = "")
  M <- ncol(DD)
  DD <- DD[, c(1, 3:M)]
  write.table(DD, file = fnm, quote = FALSE, sep = "\t", row.names = FALSE)
  sheetS <- paste(f, "-CodonCount", sep = "")
  write.xlsx(DD, xlsxFS, sheetName = sheetS, col.names = TRUE, row.names = FALSE, append = TRUE)
  ######
  RLT <- read.table(CodonPerc, header = TRUE, sep = "\t")
  RLT <- RLT[, 1:4]
  RLT$Perc <- 0
  RLT$Perc_Norm <- 0
  RLT$Count <- 0
  N <- nrow(RLT)
  g <- 1
  for(g in c(1:N)){
    str <- RLT[g, 4]
    col <- which(colnames(DD) == str)
    RLT[g, 7] <- sum(DD[[col]])
  }#g
  grpN <- 21
  g <- 1
  for(g in c(1:grpN)){
    D2a <- RLT[RLT[[2]]==g, ]
    sm <- sum(D2a[[7]])
    D2a[[5]] <- round(D2a[[7]] / sm, digits = 4)
    mx <- max(D2a[[5]])
    D2a[[6]] <- round(D2a[[5]]/mx, digits = 4)
    if(g == 1) RLT2 <- D2a
    else{
      RLT2 <- rbind(RLT2, D2a)
    }
  }#g
  fnm0 <- filestem(filenm)
  fnm <- paste(fnm0, "_codonUSTB.txt", sep = "")
  write.table(RLT2, file = fnm, quote = FALSE, sep = "\t", row.names = TRUE)
  sheetS <- paste(f, "-CodonUsage", sep = "")
  write.xlsx(RLT2, xlsxFS, sheetName = sheetS, col.names = TRUE, row.names = FALSE, append = TRUE)
}#CodonCount

#' Create a file with Motif and Restrict Sites sequence count numbers
#'
#' @param infile "DNA sequence file", "Motif define file", "Result Excel file and "Input file serial number"
#' @return None (Result file with "*_motifCount.txt" and Excel sheet will be created)
#' @export
MotifCount <- function(filenm, motifnm, xlsxFS, f){
  fnm0 <- motifnm
  D1 <- read.table(fnm0, header = TRUE, sep = "\t")
  motifList <- as.vector(D1[, 1])
  L <- length(motifList)
  fnm1 <- filenm
  DD <- read.table(fnm1, header = TRUE, sep = "\t")
  N = nrow(DD)
  M <- ncol(DD)
  l <- 1
  for(l in c(1:L)){
    DD[[l+2]] <- c(rep(0, N))
    colnames(DD)[l+2] <- motifList[l]
  }
  i <- 1
  for(i in c(1:N)){
    seqSS <- toupper(toString(DD[i, M]))
    for(j in c(1:L)){
      motifSS <- motifList[j]
      cnt <- str_count(seqSS, pattern = motifSS)
      DD[i, j+2] <- cnt
    } #j
  } #for i
  M <- ncol(DD)
  DD <- DD[, colnames(DD)[c(1, c(3:(M)), 2)]]
  New <- as.vector(D1[, 3])
  New <- New - 1
  New <- c("Maxium Allowed", New, "...")
  DD <- InsertRow(DD, NewRow = New, RowNum = 1)
  M <- ncol(DD)
  for(j in c(2:(M-1))){
    DD[[j]] <- as.numeric(DD[[j]])
  }
  # DD[, 2:(M-1)] <- as.numeric(DD[, 2:(M-1)])
  fnm0 <- filestem(filenm)
  fnm2 <- paste(fnm0, "_motifCount.txt", sep = "")
  write.table(DD, file = fnm2, quote = FALSE, sep = "\t", row.names = FALSE)
  sheetS <- paste(f, "-MotifCount", sep = "")
  write.xlsx(DD, xlsxFS, sheetName = sheetS, col.names = TRUE, row.names = FALSE, append = TRUE)
}#MotifCount

#' Create a new Codon usage table file with normalized Codon usage values on last column
#'
#' @param infile "Codon Usage Table file"
#' @return None (Result file with "*N.txt"  and Excel sheet will be created)
#' @export
normCodonPerc <- function(CodonPerc){
  D2 <- read.table(CodonPerc, header = TRUE)
  grpN <- 21
  g <- 1
  for(g in c(1:grpN)){
    D2a <- D2[D2[[2]]==g, ]
    mx <- max(D2a[[5]])
    D2a$Perc_Norm <- round(D2a[[5]]/mx, digits = 4)
    if(g == 1) RLT <- D2a
    else{
      RLT <- rbind(RLT, D2a)
    }
  }#g
  # fnm <- "CornCodonTable_V2.txt"
  # fnm <- paste(filestem(CodonPerc), "N.txt", sep = "")
  fnm <- CodonPercN
  write.table(RLT, file = fnm, quote = FALSE, sep = "\t", row.names = FALSE)
}#normCodonPerc

#' Create a file with Rare codon and GC contents information
#'
#' @param infile "DNA sequence file", "New codon Usage Table file", "Result Excel file and "Input file serial number"
#' @return None (Result file with "*_RareCD.txt"  and Excel sheet will be created)
#' @export
RareCodonCalc <- function(filenm, CodonPercN, xlsxFS, f){
  CDFv <- 0.3
  fnm <- filenm
  D1 <- read.table(fnm, header = TRUE, sep = "\t")
  D1N <- nrow(D1)
  D1$CAI <- 0
  D1$CFD <- 0
  colnames(D1)[4] <- paste("CFD < ", CDFv, sep = "")
  D1$GC_Perc <- 0
  D1$GC3_Perc <- 0

  D2 <- read.table(CodonPercN, header = TRUE, sep = "\t")
  s <- 1
  for(s in c(1:D1N)){
    str1 <- toupper(D1[s, 2])
    str2 <- stri_sub(str1, seq(1, stri_length(str1),by=3), length=3)
    strN <- length(str2)
    CAI <- 0.0
    CFD <- 0.0
    GC_Count <- 0.0
    GC3 <- 0.0
    t <- 2
    for(t in c(2:strN)){
      str3 <- str2[t]
      value <- D2[D2[[4]]==str3, 6]
      CAI <- CAI + value
      if(value < CDFv){
        CFD <- CFD + 1
        # print(str3)
      }
      GC_Count <- GC_Count + sum(charToRaw(str3) == charToRaw('G')) #str_count(str3, "G")
      GC_Count <- GC_Count + sum(charToRaw(str3) == charToRaw('C')) #str_count(str3, "C")
      if(substr(str3, 3, 3) == "G") GC3 <- GC3 +1
      if(substr(str3, 3, 3) == "C") GC3 <- GC3 +1
    }#t
    CAI <- round(CAI / (strN-1), digits = 4)
    CFD <- round(CFD / (strN-1), digits = 4)
    GC_Count <- round(GC_Count / ((strN-1)*3), digits = 4)
    GC3 <- round(GC3 / ((strN-1)), digits = 4)

    D1[s, 3] <- CAI
    D1[s, 4] <- CFD
    D1[s, 5] <- GC_Count
    D1[s, 6] <- GC3
  }#s
  D1 <- D1[, c(1, 3, 4, 5, 6, 2)]
  fnm0 <- filestem(filenm)
  fnm <- paste(fnm0, "_RareCD.txt", sep = "")
  write.table(D1, file = fnm, quote = FALSE, sep = "\t", row.names = FALSE)
  sheetS <- paste(f, "-RareCodonInf", sep = "")
  write.xlsx(D1, xlsxFS, sheetName = sheetS, col.names = TRUE, row.names = FALSE, append = TRUE)
}#RareCodonCalc

#' Create one graph file for each DNA sequence related with rare Codon information
#'
#' @param infile "DNA sequence file", "New codon Usage Table file", "Result Excel file and "Input file serial number"
#' @return None (One or more graph PDF file(s) with related DNA sequence ID will be created)
#' @export
RareCodonGraph <- function(filenm, CodonPercN, f, ranges){
  fnm <- filenm
  D1 <- read.table(fnm, header = TRUE, sep = "\t")
  D1N <- nrow(D1)
  D2 <- read.table(CodonPercN, header = TRUE, sep = "\t")
  s <- 1
  if(length(ranges) == 0) ranges <- c(1:D1N)
  for(s in ranges){
    if(s > D1N) break
    str1 <- toupper(D1[s, 2])
    str2 <- stri_sub(str1, seq(1, stri_length(str1),by=3), length=3)
    strN <- length(str2)
    CAI <- 0.0
    RLT1 <- c()
    t <- 2
    for(t in c(2:strN)){
      str3 <- str2[t]
      value <- D2[D2[[4]]==str3, 6]
      RLT1 <- c(RLT1, 1-value)
      CAI <- CAI + value
    }#t
    CAI <- round(CAI / (strN-1), digits = 4)
    mainS <- paste("(", D1[s, 1], ") Codon Adaption Index = ", CAI, sep = "")
    xlabS <- paste("Codon Position (1 - ", strN, ")", sep = "")
    fnm <- paste("CAIGrpah-File", f, "_", s, "-", D1[s, 1], ".pdf", sep = "")
    pdf(file = fnm, pointsize = 8, width = 12)
    barplot(RLT1, ylim = c(0, 1), xlab = xlabS, ylab = "Codon Usage Rate", main = mainS, plot = TRUE, axes = FALSE)#
    axis(2, at = seq(1.0, 0.0, by = -0.1), labels = seq(0.0, 1.0, by = 0.1))
    dev.off()
  }#s
}#RareCodonGraph

#' Create one graph file for each DNA sequence related with GC contents information
#'
#' @param infile "DNA sequence file", "Result Excel file and "Input file serial number"
#' @return None (One or more graph PDF file(s) with related DNA sequence ID will be created)
#' @export
GC_PercGraph <- function(filenm, f, ranges){
  bp1 = 14
  bp2 = 7
  fnm <- filenm
  D1 <- read.table(fnm, header = TRUE, sep = "\t")
  D1N <- nrow(D1)
  s <- 1
  if(length(ranges) == 0) ranges <- c(1:D1N)
  for(s in ranges){
    if(s > D1N) break
    str1 <- toupper(D1[s, 2])
    str2 <- stri_sub(str1, seq(1, stri_length(str1),by=3), length=3)
    strN <- length(str2)
    stepN <- as.integer(strN / bp2) - 2
    RLT <- c()
    w <- 2
    for(w in c(1:stepN)){
      GC_Count <- 0.0
      t1 <- (w-1)*bp2 + 1
      t2 <- t1 + bp1 - 1
      for(t in c(t1:t2)){
        str3 <- str2[t]
        GC_Count <- GC_Count + sum(charToRaw(str3) == charToRaw('G')) #str_count(str3, "G")
        GC_Count <- GC_Count + sum(charToRaw(str3) == charToRaw('C')) #str_count(str3, "C")
      }
      GC_Count <- round(GC_Count / (bp1*3), digits = 4)
      RLT <- c(RLT, GC_Count)
    }#w
    fnm <- paste("GCPercGrpah-File", f, "_", s, "-", D1[s, 1], ".pdf", sep = "")
    mainS <- paste("(", D1[s, 1], ") GC Percentage in Moving Windows", sep = "")
    xlabS <- paste("42 bp Windows with 21 bp Steps (1 - ", stepN, ")", sep = "")
    pdf(file = fnm, pointsize = 8, width = 12)
    # barplot(RLT, ylim = c(0, 0.8), ylab = "GC Percentages", xlab = xlabS, main = mainS)
    plot(RLT, ylim = c(0, 1), ylab = "GC Percentages", xlab = xlabS, main = mainS, type = 'b')
    dev.off()
  }#s
}#GC_PercGraph

#' Create a file with Codon count numbers (adjust by position) and a file with Codon usage
#'
#' @param infile "DNA sequence file", "Codon Usage Table file", "Result Excel file and "Input file serial number"
#' @return None (Result files with "*_codonSEQ.txt" and "_codonUSTB" and Excel sheets will be created)
#' @export
CodonCountPos <- function(filenm, CodonPerc, xlsxFS, f, metds = 1) {
  metds <- 1
  if(metds == 1){ #divided whole seq into segments
    nParts <- 5
    # nWeigt <- c(1.33, 1.0, 0.84, 0.83, 1.0)
    nWeigt <- c(1.8, 1.1, 0.7, 0.5, 0.9)
  }
  fnm <- filenm
  DD <- read.table(fnm, header = TRUE, sep = "\t")
  N = nrow(DD)
  M <- length(codonList)
  colnames(DD)[1] <- "Gene_ID"
  DD$RanVal <- runif(N)*10
  for(c in c(1:M)){
    DD[c + 3] <- rep(0, N)
    colnames(DD)[c + 3] <- codonList[c]
  }
  DD <- DD[, colnames(DD)[c(1, 3, c(4:(M+3)), 2)]]
  M <- ncol(DD)
  i <- 1
  for(i in c(1:N)){
    x <- toupper(toString(DD[i, M]))
    slen1 <- nchar(x) / 3
    slen2 <- as.integer(slen1 / nParts) * 3
    slen3 <- c(slen2, slen2, slen2, slen2, nchar(x)-(slen2*4))
    k <- 1
    for(k in c(1:nParts)){
      if(k == 1){
        xx <- c(substr(x, 1, slen3[k]))
        len <- slen3[k]
      }else{
        xx <- c(xx, substr(x, len+1, len+slen3[k]))
        len <- len + slen3[k]
      }
    }#k
    k <- 2
    for(k in c(1:nParts)){
      strVec <- substring(xx[k], seq(1, nchar(xx[k]), 3), seq(3, nchar(xx[k]), 3))
      ctTab <- as.data.frame(table(strVec))
      # ctTab$Freq <- ctTab$Freq*nWeigt[k]
      if(k == 1){
        j <- 3
        for(j in c(3:(M-1))){
          nm <- colnames(DD)[j]
          val <- ctTab[ctTab$strVec==nm, 2]
          if(length(val) > 0) DD[i, j] <- val*nWeigt[k]
        } #for j
      }else{
        j <- 3
        for(j in c(3:(M-1))){
          nm <- colnames(DD)[j]
          val <- ctTab[ctTab$strVec==nm, 2]
          if(length(val) > 0) DD[i, j] <- DD[i, j] + val*nWeigt[k]
        } #for j
      }
    }#k
  } #for i
  # DD[1:3, 1:15]
  fnm0 <- filestem(filenm)
  fnm <- paste(fnm0, "_codonSEQ_pos.txt", sep = "")
  M <- ncol(DD)
  DD <- DD[, c(1, 3:M)]
  write.table(DD, file = fnm, quote = FALSE, sep = "\t", row.names = FALSE)
  sheetS <- paste(f, "-CodonCount_pos", sep = "")
  write.xlsx(DD, xlsxFS, sheetName = sheetS, col.names = TRUE, row.names = FALSE, append = TRUE)
  ######
  RLT <- read.table(CodonPerc, header = TRUE, sep = "\t")
  RLT <- RLT[, 1:4]
  RLT$Perc <- 0
  RLT$Perc_Norm <- 0
  RLT$Count <- 0
  N <- nrow(RLT)
  g <- 1
  for(g in c(1:N)){
    str <- RLT[g, 4]
    col <- which(colnames(DD) == str)
    RLT[g, 7] <- sum(DD[[col]])
  }#g
  grpN <- 21
  g <- 1
  for(g in c(1:grpN)){
    D2a <- RLT[RLT[[2]]==g, ]
    sm <- sum(D2a[[7]])
    D2a[[5]] <- round(D2a[[7]] / sm, digits = 4)
    mx <- max(D2a[[5]])
    D2a[[6]] <- round(D2a[[5]]/mx, digits = 4)
    if(g == 1) RLT2 <- D2a
    else{
      RLT2 <- rbind(RLT2, D2a)
    }
  }#g
  fnm0 <- filestem(filenm)
  fnm <- paste(fnm0, "_codonUSTB_pos.txt", sep = "")
  write.table(RLT2, file = fnm, quote = FALSE, sep = "\t", row.names = TRUE)
  sheetS <- paste(f, "-CodonUsage_pos", sep = "")
  write.xlsx(RLT2, xlsxFS, sheetName = sheetS, col.names = TRUE, row.names = FALSE, append = TRUE)
}#CodonCountPos
