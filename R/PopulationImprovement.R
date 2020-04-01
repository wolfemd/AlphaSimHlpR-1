#' popImprov1Cyc function
#'
#' Function to improve a simulated breeding population by one cycle. This version takes phenotyped individuals and crosses them to create new F1
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @return A records object with a new F1 Pop-class object of progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeSimp(records, bsp, SP)
#' records <- popImprov1(records, bsp, SP)
#' 
#' @export
popImprov1Cyc <- function(records, bsp, SP){
  # Include current year phenotypes for model training?
  trainRec <- records
  if (!bsp$useCurrentPhenoTrain){
    for (stage in 1+1:bsp$nStages){
      trainRec[[stage]] <- trainRec[[stage]][-length(trainRec[[stage]])]
    }
  }
  # Select parents among all individuals
  candidates <- records[[1]]@id
  crit <- bsp$selCritPopImprov(trainRec, candidates, SP)
  if (bsp$useOptContrib){
    progeny <- optContrib(records, bsp, SP, crit)
  } else{
    parents <- records[[1]][candidates[order(crit, decreasing=T)[1:bsp$nParents]]]
    progeny <- randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP)
  }
  records[[1]] <- c(records[[1]], progeny)
  return(records)
}

#' popImprov2Cyc function
#'
#' Function to improve a simulated breeding population by one cycle. This version does two cycles of predicting F1 individuals and making new F1s
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp List of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @return A records object with the F1 Pop-class object updated with new progeny coming out of a population improvement scheme
#' 
#' @details This function uses penotypic records coming out of the product pipeline to choose individuals as parents to initiate the next breeding cycle
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeSimp(records, bsp, SP)
#' records <- popImprov2Cyc(records, bsp, SP)
#' 
#' @export
popImprov2Cyc <- function(records, bsp, SP){
  # Don't include current year (if specified) for the first cycle
  # but do include it for the second cycle
  useCurrentPhenoTrain <- bsp$useCurrentPhenoTrain
  for (cycle in 1:2){
    trainRec <- records
    if (!useCurrentPhenoTrain){
      for (stage in 1+1:bsp$nStages){
        trainRec[[stage]] <- trainRec[[stage]][-length(trainRec[[stage]])]
      }
    }
    candidates <- records[[1]]@id
    crit <- bsp$selCritPopImprov(trainRec, candidates, SP)
    parents <- records[[1]][candidates[order(crit, decreasing=T)[1:bsp$nParents]]]
    progeny <- randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP)
    records[[1]] <- c(records[[1]], progeny)
    useCurrentPhenoTrain <- TRUE
  }
  return(records)
}

#' optContrib function
#'
#' function uses optiSel to identify number of progeny, allocate mates to minimize inbreeding depression, and return progeny
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of product pipeline parameters
#' @param SP The AlphaSimR SimParam object (needed to pull SNPs)
#' @param crit Named vector of selection criterion to be maximized
#' @return Pop class object with the progeny from optimum contribution crosses
#' @details Calculate a grm of individuals with high enough crit, then maximize crit subject to a target increase of relatedness consistent with bsp$targetEffPopSize
#' @examples 
#' crit <- bv(records[[1]]); names(crit) <- records[[1]]@id
#' progeny <- optContrib(records, bsp, SP, crit)
#' @export
optContrib <- function(records, bsp, SP, crit){
  require(optiSel)
  candidates <- names(crit)[order(crit, decreasing=T)[1:bsp$nCandOptCont]]
  grm <- sommer::A.mat(pullSnpGeno(records[[1]][candidates], simParam=SP) - 1)
  grm <- grm[candidates, candidates] # Put it in the right order
  phen <- data.frame(Indiv=candidates, crit=crit[candidates])
  invisible(capture.output(cand <- optiSel::candes(phen, grm=grm, quiet=T)))
  
  Ne <- bsp$targetEffPopSize
  con <- list(
    ub.grm = 1-(1-cand$mean$grm)*(1-1/(2*Ne))
  )
  
  oc <- opticont("max.crit", cand, con, quiet=T, trace=F)$parent[, c("Indiv", "oc")]
  totOffspr <- bsp$nCrosses * bsp$nProgeny
  keep <- oc$oc > 1 / totOffspr / 4
  oc <- oc[keep,]
  grm <- grm[keep, keep]
  oc$nOffspr <- oc$oc * 2 * totOffspr
  # Make sum to 2*totOffspr: very arcane but it works
  curOffspr <- sum(round(oc$nOffspr))
  if (curOffspr != 2*totOffspr){
    nDiff <- 2*totOffspr - curOffspr
    addOrSub <- sign(nDiff)
    decim <- addOrSub * (oc$nOffspr - floor(oc$nOffspr))
    keep <- decim + (addOrSub < 0) < 0.5
    chng <- oc$Indiv[keep]; decim <- decim[keep]
    chng <- chng[order(decim, decreasing=T)[1:abs(nDiff)]]
    oc$nOffspr[oc$Indiv %in% chng] <- oc$nOffspr[oc$Indiv %in% chng] + addOrSub
  }
  oc$nOffspr <- round(oc$nOffspr)
  oc$remOffspr <- oc$nOffspr
  crossPlan <- NULL
  for (curPar in order(oc$nOffspr, decreasing=T)){
    while(oc$remOffspr[curPar] > 0){
      # find which other has minimum relationship with curPar
      mate <- which.min(grm[curPar,])
      if (mate == curPar){ # Happens if last parent somewhat related to all
        redis <- sample(nrow(crossPlan), ceiling(oc$remOffspr[mate]/2))
        redisPar <- c(crossPlan[redis,])
        crossPlan <- rbind(crossPlan[-redis,], cbind(oc$Indiv[mate], redisPar))
        oc$remOffspr[mate] <- 0
      } else{
        nProg <- min(oc$remOffspr[curPar], oc$remOffspr[mate], bsp$nProgeny)
        if (nProg > 0){
          crossPlan <- rbind(crossPlan, matrix(rep(c(oc$Indiv[curPar], oc$Indiv[mate]), each=nProg), nrow=nProg))
          oc$remOffspr[curPar] <- oc$remOffspr[curPar] - nProg
          oc$remOffspr[mate] <- oc$remOffspr[mate] - nProg
        }
        grm[curPar, mate] <- grm[mate, curPar] <- 1e6
      }
    }
  }
  progeny <- makeCross(records[[1]], crossPlan, simParam=SP)
  return(progeny)
}
