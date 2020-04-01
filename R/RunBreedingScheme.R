#' runBreedingScheme function
#'
#' function to run a two-part strategy breeding scheme. See Gaynor et al. 2017 for the general idea.
#'
#' @param replication Integer replication of running the breeding scheme
#' @param nCycles Integer number of cycles to run the breeding scheme
#' @param initializeFunc Function to initialize the breeding program.
#' @param productPipeline Function to advance the product pipeline by one generation
#' @param populationImprovement Function to improve the breeding population and select parents to initiate the next cycle of the breeding scheme
#' @param bsp  A list of breeding scheme parameters. It contains pipeline parameters: nParents, nCrosses, nProgeny, checks, nStages, errVars, nReps, nEntries, nChks, stageNames, and nCyclesToKeepRecords. It contains population parameters: nFounders, nChr, segSites, nQTL, genVar, meanDD, varDD, nSNP
#' @return A \code{records} object containing the phenotypic records retained of the breeding scheme
#' 
#' @details A wrapper to initiate the breeding program then iterate cycles of product pipeline and population improvement
#' 
#' @examples
#' bsp <- specifyPipeline()
#' bsp <- specifyPopulation(bsp)
#' records <- runBreedingScheme(replication=NULL, nCycles=5, initializeFunc=initFuncSimp, productPipeline=prodPipeSimp, populationImprovement=popImprov1, bsp)

#' @export
runBreedingScheme <- function(replication=NULL, nCycles=2, initializeFunc, productPipeline, populationImprovement, bsp){
  cat("******", replication, "\n")
  initList <- initializeFunc(bsp)
  SP <- initList$SP
  bsp <- initList$bsp
  records <- initList$records
  
  for (cycle in 1:nCycles){
    cat(cycle, " ")
    records <- productPipeline(records, bsp, SP)
    records <- populationImprovement(records, bsp, SP)
  }
  cat("\n")
  
  # Finalize the stageOutputs
  records <- lastCycStgOut(records, bsp, SP)
  return(list(records=records, bsp=bsp, SP=SP))
}
