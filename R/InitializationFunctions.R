#' initFuncSimp function
#'
#' Function to initialize simulation of a breeding program. A single additive trait is simulated. No checks are used in this scheme
#' 
#' @param bsp A list of breeding scheme parameters.  See \code{specifyPipeline} and \code{specifyPopulation} 
#' @return A list containing: 1. The simulation parameters in \code{SP}; 2. The initial records of the breeding program in \code{records}. See \code{fillPipeline} for details; 3. A completed \code{bsp} object
#' 
#' @details Creates the founders and the initial records at the beginning of the simulation of a breeding program. The initialization function you use needs to be consistent with the functions for advancing the pipeline and improving the population. At a minimum, the SP has to be able to generate a phenotype, so SP$addTraitCCC needs to be invoked. the \code{records} object needs to be "filled" in the sense that there are phenotypes of individuals in all stages of the pipeline that can be advanced to the next stage.
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#'
#' @export
initFuncSimp <- function(bsp){
  # Create haplotypes for founder population of outbred individuals
  founderPop <- runMacs2(nInd=bsp$nFounders, nChr=bsp$nChr, segSites=bsp$segSites, Ne=bsp$effPopSize)
  
  # New global simulation parameters from founder haplotypes
  SP <- SimParam$new(founderPop)
  # Additive and dominance trait architecture
  SP$addTraitA(nQtlPerChr=bsp$nQTL, var=bsp$genVar)
  # Observed SNPs per chromosome
  SP$addSnpChip(bsp$nSNP)
  
  founders <- newPop(founderPop, simParam=SP)

  records <- fillPipeline(founders, bsp, SP)
  
  return(list(SP=SP, records=records, bsp=bsp))
}

#' initFuncADChk function
#'
#' function to initialize simulation of a breeding program. A single additive-dominance trait is simulated. Check are used in this scheme
#' 
#' @param bsp A list of breeding scheme parameters.  See \code{specifyPipeline} and \code{specifyPopulation} 
#' @return A list containing: 1. The simulation parameters in \code{SP}; 2. The initial records of the breeding program in \code{records}. See \code{fillPipeline} for details; 3. A completed \code{bsp} object
#' 
#' @details Creates the founders and the initial records at the beginning of the simulation of a breeding program.
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' initList <- initializeFunc(bsp)
#' SP <- initList$SP
#' bsp <- initList$bsp
#' records <- initList$records
#'
#' @export
initFuncADChk <- function(bsp){
  # Create haplotypes for founder population of outbred individuals
  # Note: default effective population size for runMacs is 100
  founderPop <- runMacs2(nInd=bsp$nFounders, nChr=bsp$nChr, segSites=bsp$segSites, bsp$effPopSize)
  
  # New global simulation parameters from founder haplotypes
  SP <- SimParam$new(founderPop)
  # Additive and dominance trait architecture
  SP$addTraitAD(nQtlPerChr=bsp$nQTL, var=bsp$genVar, meanDD=bsp$meanDD, varDD=bsp$varDD, useVarA=FALSE)
  # Observed SNPs per chromosome
  SP$addSnpChip(bsp$nSNP)
  
  founders <- newPop(founderPop, simParam=SP)
  if (any(bsp$nChks > 0)){
    bsp$checks <- selectInd(founders, nInd=max(bsp$nChks), use="rand", simParam=SP)
  } else bsp$checks <- NULL

  records <- fillPipeline(founders, bsp, SP)
  
  return(list(SP=SP, records=records, bsp=bsp))
}

#' fillPipeline function
#'
#' function to create initial records at the start of a simulation
#'
#' @param founders Pop-class object of the founders of the breeding program
#' @param bsp A list of product pipeline parameters.  It contains nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, stageNames, and nCyclesToKeepRecords.  See \code{runBreedingScheme} for details
#' @return A \code{records} object. A list of lists containing nStages+1 lists. The first list contains one Pop-class of progeny per year of the scheme. The remaining lists contain one matrix per year that has individual id, mother, father, stage, phenotypes, and error variances. The individuals have been phenotyped using \code{setPheno}. The matrix may contain a mix of experimental and check phenotypes with different levels of replication
#' 
#' @details This is a structure for a records object that will be used to simulate breeding schemes
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' founderPop <- runMacs(nInd=bsp$nFounders, nChr=bsp$nChr, segSites=bsp$segSites)
#' SP <- SimParam$new(founderPop)
#' SP$addTraitA(nQtlPerChr=bsp$nQTL, var=bsp$genVar)
#' SP$addSnpChip(bsp$nSNP)
#' founders <- newPop(founderPop, simParam=SP)
#' bsp <- c(bsp, checks=list(NULL))
#' records <- fillPipeline(founders, bsp, SP)
#' 
#' @export
fillPipeline <- function(founders, bsp=NULL, SP){
  records <- list(randCross(founders, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP))
  for (year in 1:bsp$nStages){
    for (stage in 1:year){
      if (stage==1){ # Stage 1: F1 progeny population: random selection use pop
        nGenoRec <- nInd(records[[1]])
        nF1 <- bsp$nCrosses * bsp$nProgeny # Sample from the most-recent F1s
        indToAdv <- nGenoRec - nF1 + sort(sample(nF1, bsp$nEntries[1]))
        entries <- records[[1]][indToAdv]
        # Don't want to bother with phenotypes but want mild selection: use gv
        parents <- selectInd(entries, nInd=nInd(entries)/1.5, use="gv", simParam=SP)
        toAdd <- list(randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP))
      } else{ # Stage > 1: sort the matrix and make population of best
        sourcePop <- last(records[[stage]])
        idBest <- sourcePop$id[order(sourcePop$pheno[1:bsp$nEntries[stage-1]], decreasing=T)[1:bsp$nEntries[stage]]]
        entries <- records[[1]][idBest]
      }
      varE <- (bsp$gxeVar + bsp$errVars[stage] / bsp$nReps[stage]) / bsp$nLocs[stage]
      # reps=1 because varE is computed above
      entries <- setPheno(entries, varE=varE, reps=1, simParam=SP)
      phenoRec <- phenoRecFromPop(entries, bsp, stage)
      if(!is.null(bsp$checks) & bsp$nChks[stage] > 0){
        varE <- (bsp$gxeVar + bsp$errVars[stage] / bsp$chkReps[stage]) / bsp$nLocs[stage]
        chkPheno <- setPheno(bsp$checks[1:bsp$nChks[stage]], varE=varE, reps=1, simParam=SP)
        chkRec <- phenoRecFromPop(chkPheno, bsp, stage, checks=T)
        phenoRec <- bind_rows(phenoRec, chkRec)
      }
      toAdd <- c(toAdd, list(phenoRec))
    }#END stages
    records[[1]] <- c(records[[1]], toAdd[[1]])
    for (i in 2:length(toAdd)){
      if (i > length(records)){
        records <- c(records, list(toAdd[i]))
      } else{
        records[[i]] <- c(records[[i]], toAdd[i])
      }
    }
  }#END years
  names(records) <- c("F1", bsp$stageNames)
  return(records)
}
