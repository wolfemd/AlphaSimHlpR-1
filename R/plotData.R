#' plotRecords function
#'
#' function to make a plot from the data in \code{records}
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @return Makes a plot of the gain over cycles of selection and returns the dataset used to make that plot
#' @details \code{records} is a list of lists of populations and phenotype matrices useful for maintaining the phenotypic observations across years and stages
#' 
#' @examples
#' exptSummary <- plotRecords(records)
#' 
#' @export
plotRecords <- function(replicRecords){
  nCycKept <- bsp$nCyclesToKeepRecords
  nCyc <- bsp$nCyclesToRun
  records <- lapply(replicRecords, mean_records)
  nS2f <- floor(bsp$nStages/2)
  nS2c <- ceiling(bsp$nStages/2)
  stageOrd <- c(2* nS2f:1, 1:nS2c *2-1)
  allRec <- NULL
  for (i in 1:length(records)){
    allRec <- allRec %>% bind_rows(records[[i]] %>% as_tibble %>% dplyr::mutate(repNum=paste0("R", i), timePeriod=1:nCycKept + (nCyc - nCycKept)))
  }
  plotData <- allRec %>% tidyr::pivot_longer(bsp$stageNames, names_to="stage", values_to="meanGenoVal")
  so <- paste0("S", stageOrd); names(so) <- bsp$stageNames
  plotData <- plotData %>% dplyr::mutate(so=so[stage])
  bp <- brkdn.plot(meanGenoVal ~ so + timePeriod, data=plotData, lwd=2, cex=1, col=order(so), md="std.error", stagger=1/(nCycKept+1)/3/bsp$nStages, xlab="Time Period", ylab="Genetic Improvement", main=paste("Annual cost:", bsp$totalCosts))
  legend("topleft", legend=bsp$stageNames, col=1:bsp$nStages, lwd = 2, cex=0.5, horiz = T)
  return(dplyr::select(plotData, -so))
}
