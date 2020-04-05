#' runBreedingSchemeCycOutputs function
#'
#' function to run a two-part strategy breeding scheme. See Gaynor et al. 2017 for the general idea.
#' Adds _per cycle_ output container to records object. Currently that is used to collect the marker effects/varcomps/GEBV and UC selection criterion computed each cycle in the "popImrpovUC" function.
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
#'
runBreedingSchemeCycOutputs <- function(replication=NULL, nCycles=2, initializeFunc, productPipeline, populationImprovement, bsp){
      cat("******", replication, "\n")
      initList <- initializeFunc(bsp)
      SP <- initList$SP
      bsp <- initList$bsp
      records <- initList$records
      records[["cycleOutputs"]]<-tibble()
      for (cycle in 1:nCycles){
            cat(cycle, " ")
            records <- productPipeline(records, bsp, SP)
            records <- populationImprovement(records, bsp, SP, cycle)
      }
      cat("\n")
      # Finalize the stageOutputs
      records <- lastCycStgOut(records, bsp, SP)
      return(list(records=records, bsp=bsp, SP=SP))
}
#' popImprovUC function
#'
#' Version of the popImprov1cyc that runs selCritUCA function to improve a simulated breeding population by one cycle.
#' This version obtains additive Usefulness Criterion (UC) for crosses and creates evenly sized new families from the best crosses.
#' This version also stores the \code{records$cycleOutputs} including e.g. varcomps and SNP effects
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp A list of breeding scheme parameters
#' @param SP The AlphaSimR SimParam object
#' @param cycle An integer for the current breeding cycle. This slot was added so that it could be passed from the runBreedingSchemeCycleOutputs function for tracking purposes.
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
#' records <- popImprovUC(records, bsp, SP)
#'
#' @export
popImprovUC <- function(records, bsp, SP, cycle){
      # Include current year phenotypes for model training?
      trainRec <- records
      if (!bsp$useCurrentPhenoTrain){
            for (stage in 1+1:bsp$nStages){
                  trainRec[[stage]] <- trainRec[[stage]][-length(trainRec[[stage]])]
            }
      }
      # Select parents among all individuals
      candidates<-records$F1@id
      selCrit<-bsp$selCritPopImprov(records = trainRec, candidates, bsp, SP,
                                    nPreselect=bsp$nParents,
                                    ncores=1,propSel=0.1,h=1)
      selectedCrosses<-selCrit$selCrit[[1]] %>%
            arrange(desc(UC)) %>%
            slice(1:bsp$nCrosses) %>%
            dplyr::select(damID,sireID)
      if(bsp$useOptContrib){
            print("Um... can't do optContrib with UC right now."); stop() }
      # makeCross() works if even contributions from each cross...
      progeny<-makeCross(pop = records$F1, crossPlan = as.matrix(selectedCrosses), nProgeny = bsp$nProgeny, simParam=SP)
      records$F1 <- c(records$F1, progeny)
      records$cycleOutputs %<>% bind_rows(selCrit %>% dplyr::mutate(cycle=cycle))
      return(records)
}



#' snpeffPhenoEvalA function
#'
#' function to take a data.frame coming from framePhenoRec and GRM and analyze them with individuals as a random effect with a GRM covariance matrix
#'
#' @param phenoDF A data.frame of phenotypic observations. See \code{framePhenoRec} for details
#' @param grmA The _additive_ genomic relationship matrix
#' @param zA Centered _additive_ marker dosage matrix of use as predictors when backsolving marker effects from GBLUP solutions
#' @return tibble with slots for variance components, marker effects and BLUPs
#' @details Given all the phenotypic records calculate the GEBV and backsolve the SNP effects.
#' _Also_ computes the genomic variance according to Lehermeier et al. 2017a method "M2" (accounts for LD).
#'
#' @examples
#' phenoDF <- framePhenoRec(records, bsp)
#' grm <- makeGRM(records, bsp, SP)
#' grmBLUPs <- grmPhenoEval(phenoDF, grm)
#'
#' @export
#'
snpeffPhenoEvalA <- function(phenoDF, grmA, zA){
      require(sommer); require(predCrossVar)
      phenoDF %<>%
            dplyr::mutate(ga=factor(as.character(id),levels=rownames(grmA)),
                          WT=1/errVar)
      # gblup
      gblup<-mmer(pheno~1,
                  random=~vs(ga,Gu = grmA),
                  weights=WT,
                  data=phenoDF,verbose = F)
      ga<-as.matrix(gblup$U$`u:ga`$pheno,ncol=1)
      # backsolve marker effects
      out<-tibble(addEffects=list(backsolveSNPeff(zA,ga)),
                  GEBV=list(ga))
      # genomic variances (with LD, "M2") in the current TP
      # Not strictly necessary and has a computation cost
      # But I want to see how this estimator behaves under selection
      genoVarCovarMat<-genoVarCovarMatFunc(zA)
      genomicVa_m2<-getM2varcomp(out$addEffects[[1]],genoVarCovarMat,"add")
      # tidy output
      tidyvarcomps<-summary(gblup)$varcomp %>%
            rename(Var=VarComp) %>%
            rownames_to_column(var = "VarComp") %>%
            dplyr::mutate(VarComp=ifelse(grepl("ga",VarComp),"VarA_typeM1","VarE")) %>%
            dplyr::select(VarComp,Var) %>%
            bind_rows(tibble(VarComp="VarA_typeM2",Var=genomicVa_m2)) %>%
            arrange(VarComp) %>%
            tidyr::spread(VarComp,Var)
      out[["varcomps"]]<-list(tidyvarcomps)
      return(out)
}

#' kinshipHlpR function
#'
#' function to make either additive or dominance genomic relationship matrices using predCrossVar::kinship() at either SNPs or QTL.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object. Needed to pull the SNP genotypes
#' @param type string, "add" or "dom". type="add" gives same as rrBLUP::A.mat(), i.e. Van Raden, Method 1. type="dom" gives classical parameterization according to Vitezica et al. 2013.
#' @param snpsORqtl string, "snps" (default) or "qtl". Determines which loci are used to compute the GRM.
#' @return A genomic relationship matrix
#' @details \code{records} maintains the phenotypic and genotypic records across years and stages. For GEBV analysis, you need the GRM of these individuals. \code{makeGRM} assumes the first phenotyping stage (records[[2]]) has all individuals that have been phenotyped. The GRM also includes the unphenotyped new F1 individuals in records[[1]]
#'
#' @examples
#' grm <- makeGRM(records, bsp, SP)
#'
#' @export
kinshipHlpR<-function(records, bsp, SP, type, snpsORqtl="snps"){
      require(predCrossVar)
      allPop <- records[[1]]
      if(!is.null(bsp$checks)){
            putInChks <- setdiff(bsp$checks@id, allPop@id)
            if (length(putInChks > 0)) allPop <- c(allPop, bsp$checks[putInChks])
      }
      if(snpsORqtl=="snps"){ M<-pullSnpGeno(pop=allPop, simParam=SP) }
      if(snpsORqtl=="qtl"){ M<-pullQtlGeno(pop=allPop, simParam=SP) }
      grm<-predCrossVar::kinship(M=M,type=type)
      return(grm)
}

#' predictorsHlpR function
#'
#' function to make centered additive or dominance genomic relationship matrices e.g. for use in whole-genome regressions or backsolving SNP effects from GBLUP.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object. Needed to pull the SNP genotypes
#' @param type string, "add" or "dom". type="add" gives same as rrBLUP::A.mat(), i.e. Van Raden, Method 1. type="dom" gives classical parameterization according to Vitezica et al. 2013.
#' @param snpsORqtl string, "snps" (default) or "qtl". Determines which loci are used to compute the GRM.
#' @return A genomic relationship matrix
#' @details \code{records} maintains the phenotypic and genotypic records across years and stages. For GEBV analysis, you need the GRM of these individuals. \code{makeGRM} assumes the first phenotyping stage (records[[2]]) has all individuals that have been phenotyped. The GRM also includes the unphenotyped new F1 individuals in records[[1]]
#'
#' @examples
#'
#' @export
predictorsHlpR<-function(records, bsp, SP, type, snpsORqtl="snps"){
      require(predCrossVar)
      allPop <- records[[1]]
      if(!is.null(bsp$checks)){
            putInChks <- setdiff(bsp$checks@id, allPop@id)
            if (length(putInChks > 0)) allPop <- c(allPop, bsp$checks[putInChks])
      }
      if(snpsORqtl=="snps"){ M<-pullSnpGeno(pop=allPop, simParam=SP) }
      if(snpsORqtl=="qtl"){ M<-pullQtlGeno(pop=allPop, simParam=SP) }
      if(type=="add"){ Z<-predCrossVar::centerDosage(M) }
      if(type=="dom"){ Z<-predCrossVar::dose2domDev(M) }
      return(Z)
}

#' selCritUCA function
#'
#' function to compute the usefulness cross-selection criteria based on combo of parental GEBV and predicted _additive_ genetic variance of each cross.
#'
#' @param records The breeding program \code{records} object. See \code{fillPipeline} for details
#' @param candidates Character vector of ids of the candidates to be parents
#' @param bsp The breeding scheme parameter list
#' @param SP The AlphaSimR SimParam object (needed to pull SNPs)
#' @param snpsORqtl
#' @param nPreselect
#' @param ncores
#' @param propSel
#' @param h
#' @return tibble with 2 list-columns (selCrit: contains _per cross_ selection criterion, geneticEval: contains varcomps, marker effects and breeding values from GBLUP.)
#' @details Accesses all individuals in \code{records} to pick the highest ones
#' @examples
#' candidates <- records[[1]][[1]]@id
#' parents <- records[[1]][[1]][selCritGRM(records, candidates, bsp, SP)]
#'
#' @export
selCritUCA <- function(records, candidates, bsp, SP,snpsORqtl="snps",
                       nPreselect=bsp$nParents,ncores=1,propSel=0.1,h=1){
      require(predCrossVar)
      phenoDF <- framePhenoRec(records, bsp)
      A<-kinshipHlpR(records, bsp, SP, type="add",snpsORqtl)
      # Center and make dom. devs.
      Za<-predictorsHlpR(records, bsp, SP, type="add", snpsORqtl)
      # Remove individuals with phenotypes but who no longer have geno records
      # I am not sure this can happen but it is a safeguard
      phenoDF<-phenoDF[phenoDF$id %in% rownames(A),]
      geneticEval<-snpeffPhenoEvalA(phenoDF, grmA=A, zA=Za)
      # pre-select parents, default to nPreselect=nParents
      # means analysis will be to
      # find better crosses among same set as truncation
      candParents<-geneticEval$GEBV[[1]] %>%
            .[order(.,decreasing = T),] %>%
            .[1:nPreselect]

      crosses2predict<-crossing(sireID=names(candParents),
                                damID=names(candParents))

      if(snpsORqtl=="snps"){
            haploMat<-pullSnpHaplo(pop=records$F1[names(candParents)], simParam=SP)
            genmap<-getSnpMap(simParam = SP)
            m<-genmap$pos;
            names(m)<-paste0(genmap$chr,"_SNP_",genmap$id)
      }
      if(snpsORqtl=="qtl"){
            haploMat<-pullQtlHaplo(pop=records$F1[names(candParents)], simParam=SP)
            genmap<-getQtlMap(simParam = SP)
            m<-genmap$pos;
            names(m)<-paste0(genmap$chr,"_QTL_",genmap$id)
      }
      # Make the genome-wide recomb. freq. matrix (c1)
      recombFreqMat<-1-(2*genmap2recombfreq(m,nChr = SP$nChr))
      # table(rownames(addEffects)==rownames(recombFreqMat)) # names not matching
      # has to do with AlphaSimR... probably not a problem with real SNP IDs
      rownames(recombFreqMat)<-colnames(recombFreqMat)<-gsub(paste0(1:SP$nChr,"_",collapse = "|"),"",
                                                             rownames(recombFreqMat))
      # table(rownames(addEffects)==rownames(recombFreqMat)) # names matching!

      predictedCrosses<-runCrossVarPredsA(ped=crosses2predict,
                                          addEffects=geneticEval$addEffects[[1]],
                                          haploMat=haploMat,
                                          recombFreqMat=recombFreqMat,
                                          ncores=ncores)

      predictedCrosses %<>%
            dplyr::mutate(Nsegsnps=map_dbl(segsnps,length)) %>%
            dplyr::select(-segsnps) %>%
            left_join(tibble(sireID=names(candParents),sireGEBV=candParents)) %>%
            left_join(tibble(damID=names(candParents),damGEBV=candParents)) %>%
            dplyr::mutate(meanGEBV=(sireGEBV+damGEBV)/2)
      i <- dnorm(qnorm(1-propSel))/propSel
      predictedCrosses %<>%
            dplyr::mutate(UC=as.numeric(meanGEBV+i*h*sqrt(varA))) %>%
            arrange(desc(UC))
      cycleOutputs<-tibble(selCrit=list(predictedCrosses),
                           geneticEval=list(geneticEval))
      return(cycleOutputs)
}


