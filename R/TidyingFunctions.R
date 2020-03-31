#' removeOldestCyc function
#'
#' function to remove records of the oldest cycles still in \code{records}. Useful to avoid accumulating too much data which slows simulations down and makes them bulky
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp The breeding scheme parameter list
#' @return A \code{records} object with the first population of each list removed.
#' @details \code{records} is a list of lists. This function deletes the first object of each list, excluding the F1 list.
#'
#' @examples
#' records <- removeOldestCyc(records, bsp)
#'
#' @export
removeOldestCyc <- function(records, bsp){
  nCyclesToKeepRecords <- bsp$nCyclesToKeepRecords
  # Remove the phenotypic records that are older
  for (i in 1 + 1:bsp$nStages){
    nCycStage <- length(records[[i]])
    if (nCycStage > nCyclesToKeepRecords){
      records[[i]] <- records[[i]][-(1:(nCycStage-nCyclesToKeepRecords))]
    }
  }
  # List the id of the individuals remaining
  allID <- NULL
  for (i in 1:length(records[[2]])) allID <- c(allID, records[[2]][[i]]$id)
  for (i in 1 + 2:bsp$nStages) allID <- c(allID, records[[i]][[1]]$id)
  allID <- unique(allID)
  allID <- allID[order(as.integer(allID))]
  records[[1]] <- records[[1]][allID]
  return(records)
}

#' tidysims function
#'
#' function to create a tidy tibble from a list of runBreedingScheme simulations
#'
#' @param sims a list of AlphaSimHlpR sims (result of runBreedingScheme)
#'
#' @return a tidy tibble with one row per simulation, cols: records (tibble of all phenotypic records with indicator added for "year" and "stageName"), simulatedpop (the pop-class object from the sim), bsp and SP for the sim.
#' @export
#'
#' @examples
#' sims <- map(1:2,~runBreedingScheme(replication = .,nCycles = 7,
#'                                 initializeFunc = initFuncADChk,
#'                                 productPipeline = prodPipeFncChk,
#'                                 populationImprovement = popImprov1Cyc,
#'                                 bsp = bsp))
#' sims <- tidysims(sims)
tidysims <- function(sims){
  require(magrittr); require(purrr)
  sims <- tibble(SimRep=1:length(sims),sim=sims) %>%
    dplyr::mutate(sim=map(sim,~tibble(outType=names(.),output=.))) %>%
    tidyr::unnest(sim) %>%
    tidyr::pivot_wider(names_from = "outType",values_from = "output")
  sims %<>%
    # get the pop-class object out of "records"
    dplyr::mutate(simulatedpop=map(records,~.$F1),
                  records=map(records,~.[-1]))
  sims %<>%
    # Get the checks pop-class from "bsp" and add to "simulatedpop"
    dplyr::mutate(simulatedpop=map2(simulatedpop,bsp,~c(.y$checks,.x)),
                  records=map(records,tidyrecords))
  return(sims)
}

#' tidyrecords function
#'
#' function to create a tibble with indicator for year and stage from a "records" object produced by AlphaSimHlpR. Possibly redundant to framePhenoRec() function.
#'
#' @param records a list-class "records" object produced by AlphaSimHlpR
#'
#' @return a tibble with indicator for year and stage
#' @export
#'
#' @examples
#' records <- fillPipeline(founders, bsp, SP)
#' records <- tidyrecords(records)
tidyrecords <- function(records){
  tidyrecs <- tibble(stageName=names(records), recs=records) %>%
    unnest_longer(recs, indices_to = "year") %>%
    unnest(recs)

  return(tidyrecs)
}
