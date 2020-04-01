#' specifyPipeline function
#'
#' function specify the product pipeline. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the pipeline parameters. bsp is short for breeding sheme parameters
#' @param ctrlFileName The name of the text file with parameter values controling the simulation. Must include the path to the file. If NULL a toy example simulation will be set up
#' @return A list containing objects that specify the product pipeline. This list will determine the number of lists in the records object
#'
#' @details Call this function before beginning the simulation
#'
#' @examples
#' bsp <- specifyPipeline()
#'
#' @export
specifyPipeline <- function(bsp=NULL, ctrlFileName=NULL){
  if (is.null(ctrlFileName)){ # NULL control file: make toy example
    nStages <- 4 # Number of stages in the product pipeline
    stageNames <- c("SDN", "CET", "PYT", "AYT")
    nParents <- 20 # Number of parents in the crossing nursery
    nCrosses <- 20 # Number of crosses entering the pipeline
    nProgeny <- 10 # Number of progeny per cross
    # Don't use optimum contributions in simple default. Define other parms in case
    useOptContrib <- FALSE
    nCandOptCont <- 100
    targetEffPopSize <- 30
    # Number of number of entries in each stage
    nEntries <- c(nCrosses*nProgeny, 60, 20, 10)
    nReps <- c(1, 1, 2, 2) # Number of reps used in each stage
    nLocs <- c(1, 2, 2, 3) # Number of locations used in each stage
    # Number of checks used in each stage
    # Checks are replicated the same as experimental entries
    nChks <- c(2, 1, 1, 1)
    entryToChkRatio <- c(20, 20, 20, 10)
    # Error variances estimated from historical data
    # 200 for SDN is a guess
    errVars <- c(200, 146, 82, 40)
    names(nEntries) <- names(nChks) <- names(nReps) <- names(errVars) <- stageNames
    useCurrentPhenoTrain <- FALSE
    nCyclesToKeepRecords <- 5 # How many cycles to keep records
    # Function to advance individuals from one stage to the next
    selCritPipeAdv <- selCritIID
    selCritPopImprov <- selCritIID
    bspNew <- mget(setdiff(ls(), "bspNew"))
    #END no control file
  } else{
    parmNames <- c("nStages", "stageNames", "nParents", "nCrosses", "nProgeny", "useOptContrib", "nCandOptCont", "targetEffPopSize", "nEntries", "nReps", "nLocs", "nChks", "entryToChkRatio", "errVars", "useCurrentPhenoTrain", "nCyclesToKeepRecords", "selCritPipeAdv", "selCritPopImprov")
    bspNew <- readControlFile(ctrlFileName, parmNames)
  }
  bspNew <- calcDerivedParms(bspNew)

  bsp <- c(bsp, bspNew)
  return(bsp)
}

#' specifyPopulation function
#'
#' Function to specify the species and population characteristics. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @param ctrlFileName The name of the text file with parameter values specifying the breeding population. Must include the path to the file. If NULL a toy example simulation will be set up
#' @return A list containing objects that specify the species and population characteristics.
#'
#' @details Call this function before beginning the simulation
#'
#' @examples
#' bsp <- specifyPopulation(bsp)
#'
#' @export
specifyPopulation <- function(bsp=NULL, ctrlFileName=NULL){
  if (is.null(ctrlFileName)){ # NULL control file: make toy example
    # Species characteristics
    nChr <- 2 # Number of chromosomes
    # Population characteristics
    effPopSize <- 100 # Effective size of population generating founders
    nFounders <- 50 # Number of founders
    segSites <- 20 # Number of segregating sites per chromosome
    nQTL <- 5 # Number of QTL per chromosome
    nSNP <- 5 # Number of observed SNP per chromosome
    genVar <- 40 # Initial genetic variance
    gxeVar <- 30 # Initial genetic variance
    meanDD <- 0.8; varDD <- 0.01 # Mean and variance of dominance degree
    bspNew <- mget(setdiff(ls(), "bspNew"))
    #END no control file
  } else{
    parmNames <- c("nChr", "effPopSize", "nFounders", "segSites", "nQTL", "nSNP", "genVar", "gxeVar", "meanDD", "varDD")
    bspNew <- readControlFile(ctrlFileName, parmNames)
  }
  bsp <- c(bsp, bspNew)
  return(bsp)
}

#' specifyCosts function
#'
#' Function to specify the species and population characteristics. This would not need to be a function, but this way all definitions are in one place
#'
#' @param bsp A list of objects to combine with the species and population parameters. bsp is short for breeding sheme parameters
#' @param ctrlFileName The name of the text file with parameter values specifying the breeding costs. Must include the path to the file. If NULL a toy example simulation will be set up
#' @return A list containing objects that specify the species and population characteristics.
#'
#' @details Call this function before beginning the simulation
#'
#' @examples
#' bsp <- specifyCosts(bsp)
#'
#' @export
specifyCosts <- function(bsp=NULL, ctrlFileName=NULL){
  if (is.null(ctrlFileName)){ # NULL control file: make toy example
    # Plot costs
    plotCosts <- c(1, 8, 14, 32)
    # Per location cost
    perLocationCost <- 1000
    # Crossing cost
    crossingCost <- 0.2
    # Genotyping cost
    qcGenoCost <- 1.5
    wholeGenomeCost <- 10
    bspNew <- mget(setdiff(ls(), "bspNew"))
    #END no control file
  } else{
    parmNames <- c("plotCosts", "perLocationCost", "crossingCost", "qcGenoCost", "wholeGenomeCost")
    bspNew <- readControlFile(ctrlFileName, parmNames)
  }
  names(bspNew$plotCosts) <- bsp$stageNames
  bsp <- c(bsp, bspNew)
  # Calculate the program yearly cost
  # Assumptions
  # 1. Every new progeny is both QC and whole-genome genotyped. The number of
  # new progeny is nCrosses*nProgeny, so the cost is
  # nCrosses*nProgeny*(crossingCost + qcGenoCost + wholeGenomeCost)
  # 2. The number of plots in each trial is nEntries*nReps + nChks*chkReps so the cost
  # of the trial is (nEntries*nReps + nChks*chkReps)*plotCost*nLocs
  # NOTE for develCosts not accounting for number of rapid cycles
  develCosts <- bsp$nCrosses * bsp$nProgeny * (bsp$crossingCost + bsp$qcGenoCost + bsp$wholeGenomeCost)
  trialCosts <- ((bsp$nEntries * bsp$nReps + bsp$nChks * bsp$chkReps) * bsp$nLocs) %*% bsp$plotCost
  locationCosts <- max(bsp$nLocs) * bsp$perLocationCost
  totalCosts <- develCosts + trialCosts + locationCosts
  return(c(bsp, c(develCosts=develCosts, trialCosts=trialCosts, totalCosts=totalCosts)))
}

#' specifyBSP function
#'
#' Function to manually specify a breeding scheme parameters (bsp) object in R, rather than using a control file.
#' Currently does not handle costs. For costs (for now), run specifyCosts() on this.
#' Ideally this is useful for programmatically varying breeding schemes.
#'
#' @param schemeDF data.frame columns: stageNames, nReps, nLocs, nChks, nEntries, entryToChkRatio, errVars
#' @param nParents
#' @param nCrosses
#' @param nProgeny
#' @param useOptContrib
#' @param nCandOptCont
#' @param targetEffPopSize
#' @param useCurrentPhenoTrain
#' @param nCyclesToKeepRecords
#' @param selCritPipeAdv
#' @param selCritPopImprov
#' @param nChr
#' @param effPopSize
#' @param nFounders
#' @param segSites
#' @param nQTL
#' @param nSNP
#' @param genVar
#' @param gxeVar
#' @param meanDD
#' @param varDD
#'
#' @return a named list of of the parameters to specify a breeding scheme simulation
#'
#' @details All arguments are exactly as specified in the control files. Main exception is schemeDF, which is just a tibble() or data.frame version of the set of bsp arguments which are vectors (giving values for each breeding stage). Columns must have names exactly as in the corresponding arguments in control file: stageNames, nReps, nLocs, nChks, nEntries, entryToChkRatio, errVars
#'
#' @examples
#' schemeDF <- tibble(stageNames=c("SDN", "CET", "PYT"),
#'                  nReps=c(1, 1, 2),
#'                  nLocs=c(1, 1, 2),
#'                  nChks=c(1, 1, 2),
#'                  nEntries=c(100, 50, 20),
#'                  entryToChkRatio=c(50, 20, 10),
#'                  errVars=c(150,75,40))
#' bsp <- specifyBSP(schemeDF,nParents = 10, nCrosses = 10, nProgeny = 10,
#'                 useOptContrib = FALSE,
#'                 useCurrentPhenoTrain = TRUE,
#'                 nCyclesToKeepRecords = 1,
#'                 selCritPipeAdv = selCritGRM,
#'                 selCritPopImprov = selCritGRM,
#'                 nChr = 2,effPopSize = 50, nFounders = 100,
#'                 segSites = 100, nQTL = 5, nSNP = 10,
#'                 genVar = 50, gxeVar = 0, meanDD = 0.05, varDD = 0.25)
#' test <- runBreedingScheme(replication = 1,nCycles = 1,
#'                         initializeFunc = initFuncADChk,
#'                         productPipeline = prodPipeFncChk,
#'                         populationImprovement = popImprov1Cyc,
#'                         bsp = bsp)
#' @export
specifyBSP <- function(schemeDF,
                     nParents,nCrosses,nProgeny,
                     useOptContrib=FALSE,nCandOptCont=NULL,targetEffPopSize=NULL, # if useOptContrib=TRUE, must specify these args
                     useCurrentPhenoTrain=TRUE,
                     nCyclesToKeepRecords,
                     selCritPipeAdv,selCritPopImprov,
                     nChr,effPopSize,nFounders,
                     segSites,nQTL,nSNP,genVar,gxeVar,meanDD,varDD){
  bspNew <- list()
  bspNew[["nStages"]] <- nrow(schemeDF)
  bspNew[["stageNames"]] <- schemeDF$stageNames
  bspNew[["nReps"]] <- schemeDF$nReps %>% `names <- `(bspNew$stageNames)
  bspNew[["nLocs"]] <- schemeDF$nLocs %>% `names <- `(bspNew$stageNames)
  bspNew[["nChks"]] <- schemeDF$nChks %>% `names <- `(bspNew$stageNames)
  bspNew[["nEntries"]] <- schemeDF$nEntries %>% `names <- `(bspNew$stageNames)
  bspNew[["entryToChkRatio"]] <- schemeDF$entryToChkRatio %>% `names <- `(bspNew$stageNames)
  bspNew[["errVars"]] <- schemeDF$errVars %>% `names <- `(bspNew$stageNames)
  bspNew[["nParents"]] <- nParents
  bspNew[["nCrosses"]] <- nCrosses
  bspNew[["nProgeny"]] <- nProgeny
  bspNew[["useOptContrib"]] <- useOptContrib # if setting this true, there are other arguments that are needed
  bspNew[["useCurrentPhenoTrain"]] <- useCurrentPhenoTrain
  bspNew[["nCyclesToKeepRecords"]] <- nCyclesToKeepRecords
  bspNew[["selCritPipeAdv"]] <- selCritPipeAdv
  bspNew[["selCritPopImprov"]] <- selCritPopImprov
  bspNew[["nChr"]] <- nChr
  bspNew[["effPopSize"]] <- effPopSize
  bspNew[["nFounders"]] <- nFounders
  bspNew[["segSites"]] <- segSites
  bspNew[["nQTL"]] <- nQTL
  bspNew[["nSNP"]] <- nSNP
  bspNew[["genVar"]] <- genVar
  bspNew[["gxeVar"]] <- gxeVar
  bspNew[["meanDD"]] <- meanDD
  bspNew[["varDD"]] <- varDD

  bspNew <- calcDerivedParms(bspNew)
  return(bspNew) 
}

#' calcDerivedParms function
#'
#' Once you have read in parameters from a control file, or set them yourself, there are still a few derived parameters that are needed.  This function calculates them.
#'
#' @param bsp A list. bsp is short for breeding sheme parameters.
#' @return A list bsp that extends the input with a few derived parameters
#'
#' @details This function is only called internally by other functions used to specify the pipeline
#'
calcDerivedParms <- function(bsp){
  # Some parms have to be logical
  bsp$useCurrentPhenoTrain <- as.logical(bsp$useCurrentPhenoTrain)
  bsp$useOptContrib <- as.logical(bsp$useOptContrib)
  
  # In case the function is referred by name, replace with actual function
  if ("character" %in% class(bsp$selCritPipeAdv))
    bsp$selCritPipeAdv <- get(bsp$selCritPipeAdv)
  if ("character" %in% class(bsp$selCritPopImprov))
    bsp$selCritPopImprov <- get(bsp$selCritPopImprov)
  
  # Make sure you keep enough cycles
  bsp$nCyclesToKeepRecords <- max(bsp$nStages+1, bsp$nCyclesToKeepRecords)
  
  # Stop and warn user if not enough crosses specified
  if((bsp$nCrosses * bsp$nProgeny) < bsp$nEntries[1]){
    print("Not enough F1s to fill up Stage 1 trial. [nCrosses * nProgeny >= nEntries for Stage 1] is required")
    stop()
  }
  
  # Figure out how many checks to add to each stage
  pairwiseComp <- function(vec1, vec2, fnc){
    return(apply(cbind(vec1, vec2), 1, fnc))
  }
  nPlots <- bsp$nEntries * bsp$nReps
  nChkPlots <- nPlots / bsp$entryToChkRatio
  nChkPlots <- pairwiseComp(nChkPlots, bsp$nReps, max) # At least one check / rep
  chkReps <- ceiling(nChkPlots / bsp$nChks)
  # Safety if nChks or entryToChkRatio misunderstood
  bsp$nChks <- if_else(bsp$entryToChkRatio == 0, 0, bsp$nChks)
  chkReps <- if_else(is.infinite(chkReps) | is.nan(chkReps) | is.na(chkReps), 0, chkReps)
  
  # Make sure everything has names
  names(bsp$nEntries) <- names(bsp$nChks) <- names(bsp$nReps) <- names(bsp$nLocs) <- names(bsp$errVars) <- names(chkReps) <- bsp$stageNames
  bsp <- c(bsp, list(chkReps=chkReps), list(checks=NULL))

  return(bsp)
}

#' function to read a text control file
#'
#' The text file should be organized as follows
#' 1. Any text after a comment symbol # will be ignored
#' 2. Control parameter names should be on their own line
#' 3. Parameter values should be on the following line. If multiple parameter values are needed they should be separated by white space but on the same line
#' @param fileName The name of the text file to be read. Must include the path to the file
#' @param parmNames A string vector with the names of the control parameters that will be searched in the text file
#' @return A named list of the parameter values read from the control file
#'
#' @details Call this function before beginning the simulation
#'
#' @examples
#' params <- readControlFile("./inputDir/ctrlFile.txt", c("nStages", "nParents", "nCrosses"))
#'
#' @export
readControlFile <- function(fileName, parmNames){
  ctrlLines <- readLines(fileName)
  ctrlLines <- sapply(ctrlLines, function(st) strsplit(st, "#", fixed=T)[[1]][1])
  getParm <- function(parmName){
    parmRow <- grep(parmName, ctrlLines)+1
    parms <- unlist(strsplit(ctrlLines[parmRow], "[[:space:]]"))
    names(parms) <- NULL
    parmsNum <- suppressWarnings(as.numeric(parms))
    if (!any(is.na(parmsNum))) parms <- parmsNum
    return(parms)
  }
  parms <- lapply(parmNames, getParm)
  names(parms) <- parmNames
  return(parms)
}
