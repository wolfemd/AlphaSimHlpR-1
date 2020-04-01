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
  founderHap <- runMacs2(nInd=bsp$nEntries[1], nChr=bsp$nChr, segSites=bsp$segSites, Ne=bsp$effPopSize)
  
  # New global simulation parameters from founder haplotypes
  SP <- SimParam$new(founderHap)
  # Additive and dominance trait architecture
  SP$addTraitA(nQtlPerChr=bsp$nQTL, var=bsp$genVar)
  # Observed SNPs per chromosome
  SP$addSnpChip(bsp$nSNP)
  
  founders <- newPop(founderHap, simParam=SP)

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
  founderHap <- runMacs2(nInd=bsp$nEntries[1], nChr=bsp$nChr, segSites=bsp$segSites, bsp$effPopSize)
  
  # New global simulation parameters from founder haplotypes
  SP <- SimParam$new(founderHap)
  # Additive and dominance trait architecture
  SP$addTraitAD(nQtlPerChr=bsp$nQTL, var=bsp$genVar, meanDD=bsp$meanDD, varDD=bsp$varDD, useVarA=FALSE)
  # Observed SNPs per chromosome
  SP$addSnpChip(bsp$nSNP)
  
  founders <- newPop(founderHap, simParam=SP)
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
#' founderHap <- runMacs(nInd=bsp$nFounders, nChr=bsp$nChr, segSites=bsp$segSites)
#' SP <- SimParam$new(founderHap)
#' SP$addTraitA(nQtlPerChr=bsp$nQTL, var=bsp$genVar)
#' SP$addSnpChip(bsp$nSNP)
#' founders <- newPop(founderHap, simParam=SP)
#' bsp <- c(bsp, checks=list(NULL))
#' records <- fillPipeline(founders, bsp, SP)
#' 
#' @export
fillPipeline <- function(founders, bsp=NULL, SP){
  nF1 <- bsp$nCrosses * bsp$nProgeny
  records <- list(founders)
  for (year in 1 + -(bsp$nStages:1)){
    toAdd <- list()
    for (stage in 1:(year+bsp$nStages)){
      if (stage==1){ # Stage 1: F1 progeny population: random selection use pop
        # Select from the most recent F1s
        indToAdv <- nInd(records[[1]]) - nF1 + sort(sample(nF1, bsp$nEntries[stage]))
      } else{
        # Don't allow checks to be advanced: use 1:bsp$nEntries[stage-1]
        sourcePop <- last(records[[stage]])[1:bsp$nEntries[stage-1],]
        indToAdv <- order(sourcePop$pheno, decreasing=T)[1:bsp$nEntries[stage]]
        indToAdv <- sourcePop$id[sort(indToAdv)]
      }
      entries <- records[[1]][indToAdv]
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
    
    # Make the next F1s with mild selection using gv
    lastGen <- nInd(records[[1]]) - nF1 + 1:nF1
    parents <- selectInd(records[[1]][lastGen], nInd=nF1/1.5, use="gv", simParam=SP)
    toAdd <- c(list(randCross(parents, nCrosses=bsp$nCrosses, nProgeny=bsp$nProgeny, ignoreGender=T, simParam=SP)), toAdd)

    # Actually fill the records
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
  # stageOutputs relies on knowing the year from the previous year
  return(c(records, stageOutputs=list(tibble(year=0))))
}
