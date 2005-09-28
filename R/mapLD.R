"mapLD" <-
function(SNPdata, locusID.col, subjectID.col, allele.cols,
                  WhichGene = NA, outgraph = NA){
  
  # Begin all core function definitions.
  getLogLik <- function(h, p, q, N){
    N.AB <- N[1]
    N.Ab <- N[2]
    N.aB <- N[3]
    N.ab <- N[4]
    N.AaBb <- N[5]

    N.AB * log(h) +
      N.Ab * log(p - h) +
        N.aB * log(q - h) +
          N.ab * log(1 - p - q + h) +
            N.AaBb * log(h * (1 - p - q + h) + (p - h) * (q - h))
  }

  getGenotypeInfo <- function(thisSNP, WhichLocus){
    thisSNP[, 3:4] <- t(apply(thisSNP[, 3:4], 1, function(row.i) {
      class(row.i) <- 'character'
      row.i
    })) # To prevent base T from being interpretted as data type of 'logical';
  
    AllAlleleValues <- as.vector(unique(c(thisSNP[, 3], thisSNP[, 4])))
    howMany <- length(AllAlleleValues)
    if (howMany > 2)
      stop(paste('# of allele variants at a locus should be <= 2.',
               howMany, 'variants have been found.'))
    genotypeData <- data.frame(t(apply(thisSNP, 1, function(row.i, AllAlleleValues){
      alleleValues <- row.i[3:4]
      genotype <- NULL
      if (alleleValues[1] != alleleValues[2])
        genotype <- paste(LETTERS[WhichLocus], letters[WhichLocus], sep = '')
      else {
        if (alleleValues[1] == AllAlleleValues[1])
          genotype <- paste(rep(LETTERS[WhichLocus], 2), collapse = '')
        else
          genotype <- paste(rep(letters[WhichLocus], 2), collapse = '')
      }
      c(genotype, row.i[2])
    }, AllAlleleValues = AllAlleleValues)))
    dimnames(genotypeData)[[2]] <- c('Genotype', 'SampleID')
    no.homozygotes <- sum(genotypeData$Genotype ==
                          paste(rep(LETTERS[WhichLocus], 2), collapse = ''))
    no.heterozygotes <- sum(genotypeData$Genotype ==
                            paste(LETTERS[WhichLocus],
                                  letters[WhichLocus], sep = ''))
    list(genotype = genotypeData, locus = WhichLocus,
         fFirst = (2 * no.homozygotes + no.heterozygotes) /
         (2 * nrow(genotypeData)))
  }

  get2LociInfo <- function(SNPdata.ij, locusID.col, subjectID.col, allele.cols){ 
    noneMissing <- SNPdata.ij[!is.na(SNPdata.ij[, allele.cols[1]]) &
                              !is.na(SNPdata.ij[, allele.cols[2]]),
                              c(locusID.col, subjectID.col, allele.cols)]
    AllLociValues <- as.vector(unique(noneMissing[, 1]))

    locus1Info <- getGenotypeInfo(thisSNP = noneMissing[noneMissing[, 1] ==
                                    AllLociValues[1], ],
                                  WhichLocus = 1)
    locus2Info <- getGenotypeInfo(thisSNP = noneMissing[noneMissing[, 1] ==
                                    AllLociValues[2], ],
                                  WhichLocus = 2)
 
    twoLociGenotype <- merge(locus1Info$genotype, locus2Info$genotype,
                             by = 'SampleID', sort = TRUE)
    obsGcounts <- table(paste(twoLociGenotype$Genotype.x,
                              twoLociGenotype$Genotype.y, sep = ''))

    allGcounts <- rep(0, 9)
    names(allGcounts) <- c('AABB', 'AaBB', 'aaBB', 'AABb', 'AAbb',
                           'AaBb', 'Aabb', 'aaBb', 'aabb')
    for (G in 1:length(allGcounts)) {
      for (obsG in 1:length(obsGcounts)){
        if (names(allGcounts)[G] == names(obsGcounts)[obsG]) {
          allGcounts[G] <- obsGcounts[obsG]
          break
        }
      }
    }
    N.AB <- 2 * allGcounts[names(allGcounts) == 'AABB'] +
      allGcounts[names(allGcounts) == 'AABb'] +
        allGcounts[names(allGcounts) == 'AaBB']
    N.Ab <- 2 * allGcounts[names(allGcounts) == 'AAbb'] +
      allGcounts[names(allGcounts) == 'AABb'] +
        allGcounts[names(allGcounts) == 'Aabb']
    N.aB <- 2 * allGcounts[names(allGcounts) == 'aaBB'] +
      allGcounts[names(allGcounts) == 'AaBB'] +
        allGcounts[names(allGcounts) == 'aaBb']
    N.ab <- 2 * allGcounts[names(allGcounts) == 'aabb'] +
      allGcounts[names(allGcounts) == 'Aabb'] +
        allGcounts[names(allGcounts) == 'aaBb']
    N.AaBb <- allGcounts[names(allGcounts) == 'AaBb']

    gameteN <- c(N.AB, N.Ab, N.aB, N.ab, N.AaBb)
    names(gameteN) <- c('AB', 'Ab', 'aB', 'ab', 'AaBb')
    list(N = gameteN, p = locus1Info$fFirst, q = locus2Info$fFirst)
  }

  getDprime <- function(p, q, N){
    h <- optimize(getLogLik, c(max(p + q - 1, 0) + 10^-8, min(p, q) - 10^-8),
                  tol = .Machine$double.eps^0.25,
                  p = p, q = q, N = N, maximum = TRUE)$maximum

    list(Dprime = ifelse(h > p * q, abs(h - p * q) /
           min(q * (1 - p), p * (1 - q)),
           abs(h - p * q) / min(p * q, (1 - p) * (1 - q))),
         HapFreqs = c(AB = h, Ab = p - h, aB = q - h, ab = 1 + h - p - q))
  }

  calNULLh <- function(p, q, Dprime, N){
    possibleH <- c(p * q + Dprime * min(q * (1 - p), p * (1 - q)),
                 p * q - Dprime * min(p * q, (1 - p) * (1 - q)))
    if (possibleH[1] != possibleH[2])
      h <- possibleH[possibleH < min(p, q) & possibleH > max(p + q - 1, 0)]
    else
      h <- possibleH[1]
    if (length(h) > 1) {
      logLik <- getLogLik(h, p, q, N)
      if (logLik[1] == logLik[2])
        h[1]
      else
        h[logLik == max(logLik)]
    }
    else
      h
  }

  calGabrielCI <- function(p, q, N){
    seqK <- c(0, seq(0.001, 0.999, by = 0.001), 0.9999)
    allHs <- sapply(seqK, calNULLh, p = p, q = q, N = N)
    allLogLik <- sapply(allHs, getLogLik, p = p, q = q, N = N)
    scaled.LogLik <- allLogLik - max(allLogLik)
    proportion <- cumsum(exp(scaled.LogLik)) / sum(exp(scaled.LogLik))
    left <- max(proportion[proportion <= 0.05])
    if (is.infinite(left))
      lower <- 0
    else
      lower <- seqK[proportion == left]
    right <- min(proportion[proportion >= 0.95])
    if (is.infinite(right))
      upper <- 1
    else
      upper <- seqK[proportion == right]
    c(lower, upper)
  }

  is19Reached <- function(subData, afterData){
    sCount <- 0
    hCount <- 0
    loci <- unique(subData$LOCUSID2)
    for (id1 in 1:(length(loci) - 1)) {
      for (id2 in (id1 + 1):length(loci)) {
        id12 <- afterData$LDgroup[afterData$LOCUSID1 == loci[id1]
                                  & afterData$LOCUSID2 == loci[id2]]
        if (length(id12) > 0){
          if (id12 == 'S')
            sCount <- sCount + 1
          if (id12 == 'H')
            hCount <- hCount + 1
        }
      }
    }
    sCount <- sCount + sum(subData$LDgroup == 'S')
    hCount <- hCount + sum(subData$LDgroup == 'H')
    if (sCount >= 19 * hCount)
      TRUE
    else
      FALSE
  }

  recursiveBreak <- function(i, sub.i, LDdata){
    tempTails <- sub.i$LOCUSID2[sub.i$LDgroup == 'S']
    tail <- ifelse(length(tempTails) > 0, max(tempTails), i)
    afterData <- LDdata[LDdata$LOCUSID1 <= tail, ]
    subHT <- sub.i[sub.i$LOCUSID2 <= tail, ]
  
    if (nrow(subHT) == 0)
      c(i, i)
    else if (nrow(subHT) == 1){
      if (subHT$LDgroup == 'S')
        c(i, tail)
      else
        c(i, i)
    }
    else if (is19Reached(subHT, afterData))
      c(i, tail)
    else {
      choppedSubH <- subHT[subHT$LOCUSID2 < tail, ]
      recursiveBreak(i, choppedSubH, LDdata)
    }
  }

  genBlock <- function(LDdata){ 
    mList <- unique(LDdata$LOCUSID1)
    blockMap <- NULL
    HT <- NULL
    count <- 0
  
    for (i in mList){
      if (i > count){
        sub.i <- LDdata[LDdata$LOCUSID1 == i, ]
        HT <- matrix(recursiveBreak(i, sub.i, LDdata), nrow = 1)
        blockMap <- rbind(blockMap, HT)
        count <- HT[1, 2]
        HT <- NULL
      }
    }
    last1 <- max(LDdata$LOCUSID2)
    if (max(blockMap[, 2]) < last1)
      blockMap <- rbind(blockMap, matrix(rep(last1, 2), nrow = 1))
    blockMap <- data.frame(blockMap)
    dimnames(blockMap)[[2]] <- c('head', 'tail')
    blockMap
  }


  plotLD <- function(DPrimeData, outgraph, gene, BlockData){
    loci <- unique(c(DPrimeData$LOCUSID1, DPrimeData$LOCUSID2))
    x <- matrix(NA, nrow = length(loci), ncol = length(loci))
    for (i in 1:(length(loci) - 1)) {
      for (j in (i + 1):length(loci)){
        value <- DPrimeData$DPRIME[DPrimeData$LOCUSID1 == loci[i] &
                                   DPrimeData$LOCUSID2 == loci[j]]
        if (length(value) > 0)
          x[i, j] <- value
      }
    }
    for (k in 1:length(loci)){
      x[k, k] <- 1
    }

    plotData <- expand.grid(x = 1:length(loci), y = 1:length(loci))
    plotData$dprime <- as.vector(x)
    
    newBlock <- t(apply(BlockData, 1, function(loci, row.i){
      n <- length(loci)
      start <- (1:n)[loci == row.i[1]]
      end <- (1:n)[loci == row.i[2]]
      c(start, end)
    }, loci = loci))

    library(lattice)
    postscript(outgraph)
    nCol <- length(trellis.par.get("regions")$col)

    newBg <- trellis.par.get('background')
    newBg$col <- 'white'
    trellis.par.set('background', newBg)

    mainGraph <- levelplot(dprime ~ x * y, data = plotData,
                           scales = list(x = list(labels = paste(loci), at = 1:length(loci)),
                             y = list(labels = paste(loci), at = 1:length(loci))),
                           pretty = TRUE, xlab = 'Locus #', ylab = 'Locus #',
                           main = paste('LD Heatmap', ifelse(is.na(gene), '',
                             paste('in Gene', gene))),
                           col.regions = cm.colors(nCol))
  
    print(mainGraph, more = TRUE, split = rep(1, 4))
  
    additional <- levelplot(dprime ~ x * y, data = plotData,
                            scales = list(x = list(labels = paste(loci), at = 1:length(loci)),
                              y = list(labels = paste(loci), at = 1:length(loci))),
                            pretty = TRUE, xlab = 'Locus #', ylab = 'Locus #',
                            main = paste('LD Heatmap', ifelse(is.na(gene), '',
                              paste('in Gene', gene))),
                            col.regions = cm.colors(nCol),
                            panel = function(x, y, z, subscripts, blockMap, ...){
                              for (i in 1:nrow(blockMap)) {
                                row.i <- blockMap[i, ]
                                llines(c(row.i[1] - 0.5, row.i[2] + 0.5),
                                       rep(row.i[2], 2) + 0.5, col = 'blue', lwd = 3)
                                llines(rep(row.i[1], 2) - 0.5,
                                       c(row.i[1] - 0.5, row.i[2] + 0.5), col = 'blue', lwd = 3)
                              }
                            }, blockMap = newBlock)
    print(additional, more = FALSE, split = rep(1, 4))
    dev.off()
  }

  mapColIndex <- function(indicator, dataNames){
    if (is.numeric(indicator)){
      if (ceiling(indicator) == floor(indicator))
        indicator
      else
        stop('Invalid Column Index.') 
    }
    else if (is.character(indicator)){
        getName <- (1:(length(dataNames)))[dataNames == indicator]
        if (!is.na(getName))
          getName
        else
          stop('Invalid Column Index.')
      }
    else
      stop('Invalid Column Index.') 
  }
  # Finish all core function definitions.

  fieldNames <- names(SNPdata)
  locusID.col <- mapColIndex(locusID.col, fieldNames)
  subjectID.col <- mapColIndex(subjectID.col, fieldNames)
  allele.cols <- sapply(allele.cols, mapColIndex, dataNames = fieldNames)
  
  allLOCI <- as.vector(unique(SNPdata[, locusID.col]))
  if (length(allLOCI) > 1){
    LD <- NULL
    for (i in 1:(length(allLOCI) - 1)){
      LocusA <- SNPdata[SNPdata[, locusID.col] == allLOCI[i], ]
      for (j in (i + 1):length(allLOCI)){
        LocusB  <- SNPdata[SNPdata[, locusID.col] == allLOCI[j], ]
        LOCIAB <- rbind(LocusA, LocusB)
        allNeeded <- get2LociInfo(SNPdata.ij = LOCIAB, locusID.col,
                                subjectID.col, allele.cols)
        p <- allNeeded$p
        if (p < 0.1 | p > 0.9)
          break
        else {
          q <- allNeeded$q
          if (q <= 0.9 & q >= 0.1) {
            N <- allNeeded$N
            dprimes <- getDprime(p, q, N)
            coverageCI <- calGabrielCI(p, q, N)
            LD <- rbind(LD, matrix(c(i, j, dprimes$HapFreqs,
                                     dprimes$Dprime, coverageCI), nrow = 1))
          }
        }
      }
    }
    if (!is.null(LD)){
      LD <- data.frame(LD)
      dimnames(LD)[[2]] <- c('LOCUSID1', 'LOCUSID2',
                             'f.AB', 'f.Ab', 'f.aB', 'f.ab',
                             'DPRIME', 'Lower', 'Upper')
      LD$LDgroup <- ifelse(LD$Upper < 0.9, 'H',
                         ifelse(LD$Upper >= 0.98 & LD$Lower >= 0.7, 'S', 'O'))
      LD$LDstrength <- ifelse(LD$Upper < 0.9, 'Historical Evidence of Recombination',
                          ifelse(LD$Upper >= 0.98 & LD$Lower >= 0.7,
                                 'Strong LD', 'Others'))
      LDblock <- genBlock(LD)
      if (is.na(outgraph))
        outgraph <- paste(getwd(), '/LDmap.eps', sep = '')
      plotLD(DPrimeData = LD, outgraph = outgraph,
                   gene = WhichGene, BlockData = LDblock)
      print(paste('LD map has been printed to file', outgraph))
      # shell.exec(outgraph)

      LDblock$contains <- apply(LDblock, 1, function(row.i, fromLDinfo){
			if (row.i[1] == row.i[2])
                          paste(row.i[1])
			else {
                          contains <- c(row.i[1], fromLDinfo$LOCUSID2[
                            fromLDinfo$LOCUSID1 == row.i[1] &
                            fromLDinfo$LOCUSID2 <= row.i[2]])
                          paste(contains, collapse = ', ')
			}
                      }, fromLDinfo = LD[, c('LOCUSID1', 'LOCUSID2')])
      
      list(LDinfo = LD[, names(LD) != 'LDgroup'],
           LDblock = LDblock,
           locusMap = data.frame(LocusName = allLOCI,
             LocusIndex = 1:length(allLOCI)))
    }
  }
  else {    
    print("Less than two loci found. No LD block is generated.")
    return(NULL)
  }
}
