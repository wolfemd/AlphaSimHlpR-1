#' prodPipeSimp function
#'
#' Simple function to advance a simulated breeding product pipeline forward by one generation. No selection function, no checks. Selection on phenotype
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#'
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
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
prodPipeSimp <- function(records, bsp, SP){
  toAdd <- list()
  for (stage in 1:bsp$nStages){
    if (stage==1){ # Stage 1: F1 progeny population: take them at random
      nGenoRec <- nInd(records$F1)
      nF1 <- bsp$nCrosses * bsp$nProgeny
      indToAdv <- nGenoRec - nF1 + sort(sample(nF1, bsp$nEntries[1]))
    } else{ # Stage > 1: sort the matrix and make population of best
      sourcePop <- last(records[[stage]])
      indToAdv <- sourcePop$id[order(sourcePop$pheno, decreasing=T)[1:bsp$nEntries[stage]]]
    }
    entries <- records$F1[indToAdv]
    varE <- (bsp$gxeVar + bsp$errVars[stage] / bsp$nReps[stage]) / bsp$nLocs[stage]
    # reps=1 because varE is computed above
    entries <- setPheno(entries, varE=varE, reps=1, simParam=SP)
    phenoRec <- phenoRecFromPop(entries, bsp, stage)
    toAdd <- c(toAdd, list(phenoRec))
  }
  for (stage in 1 + 1:bsp$nStages){
    records[[stage]] <- c(records[[stage]], toAdd[stage-1])
  }

  # Remove old records if needed
  if (length(records[[2]]) > bsp$nCyclesToKeepRecords) records <- removeOldestCyc(records, bsp)

  return(records)
}

#' prodPipeFncChk function
#'
#' function to advance a simulated breeding product pipeline forward by one generation. See Gaynor et al. 2017 for the general idea.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has new records created by advancing by a generation
#'
#' @details The breeding program product pipeline will have been set by initializeFunc. This function moves the breeding program along by one generation and saves all the resulting phenotypes to the records object.
#'
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#' records <- prodPipeFncChk(records, bsp, SP)
#' records <- popImprov1(records, bsp, SP)
#'
#' @export
prodPipeFncChk <- function(records, bsp, SP){
  # Calculate the selection criterion. selCritPipeAdv has to be given in bsp
  candidates <- records$F1@id
  selCrit <- bsp$selCritPipeAdv(records, candidates, bsp, SP)

  # Make summary for the incoming F1s
  year <- max(records$stageOutputs$year)+1 # Add a year relative to last year
  nF1 <- bsp$nCrosses * bsp$nProgeny
  nGenoRec <- nInd(records$F1)
  newF1Idx <- nGenoRec - nF1 + 1:nF1
  id <- records$F1[newF1Idx]@id
  records$stageOutputs <- records$stageOutputs %>%
    bind_rows(stageOutputs(id=id, f1=records$F1, selCrit=selCrit, stage=0, year=year))

  # Will be added to the phenotype records
  toAdd <- list()
  for (stage in 1:bsp$nStages){
    # Make a summary for this stage
    id <- last(records[[stage+1]])$id[1:bsp$nEntries[stage]]
    records$stageOutputs <- records$stageOutputs %>%
      bind_rows(stageOutputs(id=id, f1=records$F1, selCrit=selCrit, stage=stage, year=year))

    if (stage == 1){ # Stage 1 different: no phenotypes but full Pop-class
      # Sample from the most-recent F1s
      indToAdv <- records$F1@id[nGenoRec - nF1 + sort(sample(nF1, bsp$nEntries[stage]))]
    } else{
      # Don't allow checks to be advanced: use 1:bsp$nEntries[stage-1]
      id <- last(records[[stage]])$id[1:bsp$nEntries[stage-1]]
      selCritPop <- selCrit[id]
      indToAdv <- (selCritPop %>% order(decreasing=T))[1:bsp$nEntries[stage]]
      indToAdv <- names(selCritPop)[sort(indToAdv)]
    }
    entries <- records$F1[indToAdv]
    varE <- (bsp$gxeVar + bsp$errVars[stage] / bsp$nReps[stage]) / bsp$nLocs[stage]
    # reps=1 because varE is computed above
    entries <- setPheno(entries, varE=varE, reps=1, simParam=SP)
    phenoRec <- phenoRecFromPop(entries, bsp, stage)
    # If provided, add checks to the population
    if(!is.null(bsp$checks) & bsp$nChks[stage] > 0){
      varE <- (bsp$gxeVar + bsp$errVars[stage] / bsp$chkReps[stage]) / bsp$nLocs[stage]
      chkPheno <- setPheno(bsp$checks[1:bsp$nChks[stage]], varE=varE, reps=1, simParam=SP)
      chkRec <- phenoRecFromPop(chkPheno, bsp, stage, checks=T)
      phenoRec <- bind_rows(phenoRec, chkRec)
    }
    toAdd <- c(toAdd, list(phenoRec))
  }#END 1:nStages
  for (stage in 1 + 1:bsp$nStages){
    records[[stage]] <- c(records[[stage]], toAdd[stage-1])
  }

  # Remove old records if needed
  if (length(records[[2]]) > bsp$nCyclesToKeepRecords) records <- removeOldestCyc(records, bsp)

  return(records)
}

#' stageOutputs function
#'
#' To create a row for a tibble for each cycle and each stage. The tibble will be kept as the last object of the records list
#'
#' @param id Vector of the AlphaSimR ids of the individuals in the stage
#' @param f1 The AlphaSimR pop class with all the individuals
#' @param selCrit Named vector of the selection criterion being used to advance individuals
#' @param stage Integer stage (1 to bsp$nStages) being summarized
#' @param year The current year of the breeding scheme
#' @return A tibble with whatever information from the data you want to store for analysis after simulation is done
#'
#' @details Trying to provide some flexibility in what results AlphaSimHlpR generates from a given simulation.
#'
#' @examples
#' records$stageOutputs <- records$stageOutputs %>% bind_rows(stageOutputs(id, records$F1, selCrit, stage, year))
#'
stageOutputs <- function(id, f1, selCrit, stage, year){
  stageName <- c("F1", bsp$stageNames)[stage+1]
  f1 <- f1[id]
  selCrit <- selCrit[id]
  bestCrit <- which.max(selCrit)
  if (length(bestCrit) == 0){
    gvOfBestCrit <- NA
  } else{
    gvOfBestCrit <- gv(f1[names(selCrit)[bestCrit]])
  }
  highestGV <- max(gv(f1))
  return(tibble(cycle=year-stage, year=year, stage=stageName, first=first(id), last=last(id), genValMean=mean(gv(f1)), genValSD=sd(gv(f1)), evalAtSelMean=mean(selCrit, na.rm=T), evalAtSelSD=sd(selCrit, na.rm=T), accAtSel=cor(gv(f1), selCrit), gvOfBestCrit=gvOfBestCrit, highestGV=highestGV))
}

#' lastCycStgOut function
#'
#' stageOutputs are always made on the previous cycle, so to get the end, you need to run this once
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP the AlphaSimR SimParam object
#' @return A records object that has updated stageOutputs
#'
#' @details Will figure this out later...
#'
#' @examples
#' records <- lastCycStgOut(records, bsp, SP)
lastCycStgOut <- function(records, bsp, SP){
  # Pretty it up by removing the very first row that just had year in it
  records$stageOutputs <- records$stageOutputs[-1,]

  year <- max(records$stageOutputs$year)+1 # Add a year relative to last year

  # Calculate the selection criterion. selCritPipeAdv has to be given in bsp
  candidates <- records$F1@id
  selCrit <- bsp$selCritPipeAdv(records, candidates, bsp, SP)

  # Make summary for the incoming F1s
  nF1 <- bsp$nCrosses * bsp$nProgeny
  nGenoRec <- nInd(records$F1)
  newF1Idx <- nGenoRec - nF1 + 1:nF1
  id <- records$F1[newF1Idx]@id
  records$stageOutputs <- records$stageOutputs %>% bind_rows(stageOutputs(id=id, f1=records$F1, selCrit=selCrit, stage=0, year=year))

  for (stage in 1:bsp$nStages){
    # Make a summary for this stage
    id <- last(records[[stage+1]])$id[1:bsp$nEntries[stage]]
    records$stageOutputs <- records$stageOutputs %>% bind_rows(stageOutputs(id=id, f1=records$F1, selCrit=selCrit, stage=stage, year=year))
  }

  return(records)
}
