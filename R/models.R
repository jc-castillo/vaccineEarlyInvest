#' Load data
#'
#' Load summary data with characteristics for vaccines
#'
#' @param par Parameters object with model parameters
#' @param includeVaccines List of experimental vaccines to include in analysis
#'
#' @return Data.table with summary data
#' @export
#' @import data.table
loadData <- function(par, includeVaccines=c()) {
  # We set encoding to BOM UTF to avoid cross-platform issues

  if (par$inputfile=="Default") {
    d <- data.table(read.csv("Data/vaccinesSummary.csv", fileEncoding = "UTF-8-BOM"))
  } else if (par$inputfile=="US") {
    d <- data.table(read.csv("Data/vaccinesSummaryUS.csv", fileEncoding = "UTF-8-BOM"))
  }

  d[is.na(PreClinicalCandidates), PreClinicalCandidates := 0]
  d[is.na(Phase1candidates), Phase1candidates := 0]
  d[is.na(Phase2candidates), Phase2candidates := 0]
  d[is.na(RepurposedCandidates), RepurposedCandidates := 0]

  if ("X" %in% names(d)) {
    d[, X := NULL]
  }

  dchoose <- data.table(read.csv("Data/vaccinesChoose.csv", fileEncoding = "UTF-8-BOM"))
  d <- d[.(dchoose$Platform, dchoose$Subcategory), on=.(Platform, Subcategory)]

  # Delete rows of vaccines to omit in analysis
  if (!"BCG" %in% includeVaccines) {
    d <- d[Subcategory != "Live attenuated bacteria"]
  }
  if (!"Self-assembling" %in% includeVaccines) {
    d <- d[Subcategory != "Self-assembling vaccine"]
  }
  if (!"Dendritic" %in% includeVaccines) {
    d <- d[Subcategory != "Dendritic cells"]
  }
  if (!"AAPC" %in% includeVaccines) {
    d <- d[Subcategory != "Artificial antigen presenting cells"]
  }
  if (!"Unknown" %in% includeVaccines) {
    d <- d[Subcategory != "Unknown"]
  }

  return(d)
}

#' Candidates with fungible costs
#'
#' Based on summary data for candidates, finds the optimal order to select those candidades, as well as some statistics about
#' them. Some of those statistics are computed by a Monte Carlo simulation.
#'
#' @param d Summary data with candidates
#' @param par Parameters object with model parameters
#' @param computeExpComp Whether to compute expected competition for all candidates for all numbers of potential competitors
#'
#' @return List with two data.tables. `dordered` lists all the candidates in the optimal order and with some statistics.
#' `dcanddraws` includes all the random draws and their outcomes.
#' @export
candidatesFung <- function(d, par, computeExpComp=F, seed=10) {

  # Adding platform feasibility probabilities to table
  d[Platform == "DNA", pplat := as.numeric(par$pdna)]
  d[Platform == "RNA", pplat := par$prna]
  d[Platform == "Live attenuated virus", pplat := par$pattenuated]
  d[Platform == "Viral vector", pplat := par$pvector]
  d[Platform == "Protein subunit", pplat := par$psubunit]
  d[Platform == "Inactivated", pplat := par$pinactivated]
  d[Platform == "VLP", pplat := par$pvlp]
  d[Platform == "Dendritic cells", pplat := par$pdendritic]
  d[Platform == "Self-assembling vaccine", pplat := par$psav]
  d[Platform == "Unknown", pplat := par$punknown]
  d[Platform == "Artificial antigen presenting cells", pplat := par$paapc]
  d[Platform == "Live-attenuated bacteria", pplat := par$plivebac]
  d$ptarget<-1.0
  d[Target == "Other", ptarget := par$potherprotein]
  d[Target == "Spike", ptarget := par$pspike]
  d[Target == "Recombinant", ptarget := par$precombinant]

  # Create summary table with candidates by subcategory and phase of trials
  dscPhase <- crossJoin(d, data.table(phase=c("Pre-clinical", "Phase 1", "Phase 2", "Repurposed")))
  dscPhase[phase=="Pre-clinical", Candidates := PreClinicalCandidates]
  dscPhase[phase=="Phase 1", Candidates := Phase1candidates]
  dscPhase[phase=="Phase 2", Candidates := Phase2candidates]
  dscPhase[phase=="Repurposed", Candidates := RepurposedCandidates]

  dscPhase[, PreClinicalCandidates := NULL]
  dscPhase[, Phase1candidates := NULL]
  dscPhase[, Phase2candidates := NULL]
  dscPhase[, RepurposedCandidates := NULL]
  dscPhase$pcand<-1.0
  dscPhase[phase == "Pre-clinical", pcand := par$ppreclinical]
  dscPhase[phase == "Phase 1", pcand := par$pphase1]
  dscPhase[phase == "Phase 2", pcand := par$pphase2]
  dscPhase[phase == "Repurposed", pcand := par$prepurposed]

  # Finding optimal order of vaccines ----

  # Start with no candidates
  dscPhase[, selected := 0]
  dscPhase[, psuccess := 1-(1-pcand)^selected]

  # Vector to store indices of selected vaccines
  selectIndices <- rep(NA, par$maxcand)

  # Vector to store cumulative probabilities
  culProb <- rep(NA, par$maxcand)

  # Vector to store welfare
  welfareVec <- rep(NA, par$maxcand)
  costVec <- rep(NA, par$maxcand)
  benefitVec <- rep(NA, par$maxcand)
  platFungVec <- rep(NA, par$maxcand)
  overallFungVec <- rep(NA, par$maxcand)

  # Only vaccines that have at least one candidate
  dscPhase <- dscPhase[Candidates > 0]

  #Change NAs to 0s
  dscPhase[is.na(dscPhase$pplat),"pplat"]<-0

  # Sorting data and creating index variable
  setkey(dscPhase, Platform, Subcategory, phase)
  setindex(dscPhase, Platform, Subcategory)
  dscPhase[, index := .I]

  # Compute possible target protein success combinations
  targets<- c("Spike","Recombinant","Other")
  probs<- c(as.numeric(par$pspike), as.numeric(par$precombinant), as.numeric(par$potherprotein))
  perm<-permutations(2,length(targets),v=c(0,1),repeats.allowed=TRUE)
  p_perm<- vector(mode="numeric", length=nrow(perm))
  temp<-matrix(0,nrow=nrow(perm),ncol=ncol(perm))
  for (i in 1:length(probs)){
    temp[,i]<- if_else(perm[,i]==1,probs[i], 1-probs[i])
  }
  for (i in 1:length(p_perm)){
    p_perm[i]<- prod(temp[i,])
  }
  perm <- as.data.frame(perm)
  perm <- perm[p_perm!=0,]
  names(perm) <- targets
  perm$p_perm <- p_perm[p_perm!=0]
  perm$perm_index  <-  seq.int(nrow(perm))
  perm <- data.table(perm)

  # Benefits and costs
  # benefit <- par$totmonthben * benefitIntegral(par$capacity * par$TT / 1000 / par$effpop, par) *
  #   par$effpop / par$capacity * 1000

  progben <- benefits(par$totmonthben, list(par$capacity/1000, par$afterCapacity/1000),
                      c(par$TT, par$tau), par)
  noProgBen <- benefits(par$totmonthben, list(1e-10, par$counterCapacity/1000), c(par$TT, par$tau), par)
  #benefit <- benefits(par$totmonthben, list(par$capacity/1000, par$afterCapacity/1000), c(par$TT, par$tau), par)
  benefit <- progben - noProgBen

  basemgcost <- par$c * par$capacity * 12 / 1000

  # Main loop. Sequentially find the next candidate
  for (j in 1:par$maxcand) {

    # get platform and subcategory counts
    dscPhase[, plat_count := sum(selected), by=c("Platform")]
    dscPhase[, subcat_count := sum(selected), by=c("Subcategory")]
    dscPhase[, prot_vector_count := 0]
    dscPhase[Platform %in% c("Viral vector", "Protein subunit"), prot_vector_count := sum(selected)]

    # get candidates
    dcand <- dscPhase[selected < Candidates, .(cand=index, plat_count, subcat_count, prot_vector_count)]

    # stack data frames of increments
    df_stacked <- crossJoin(dcand, dscPhase)
    df_stacked[index==cand, selected:=selected+1]

    # stack data frames of possible permutations
    df_stacked <- crossJoin(perm, df_stacked)

    # drop according to permutation
    df_stacked<- df_stacked[!(Spike==0 & Target=="Spike") &
                              !(Recombinant==0 & Target=="Recombinant") &
                              !(Other==0 & Target=="Other")]

    # get probabilities
    df_stacked[, psuccess := 1-(1-pcand)^selected]
    df_stacked <- df_stacked[, .(psuccess = par$psubcat * (1 - prod(1-psuccess))), by=c("Platform", "Subcategory", "pplat","cand","prot_vector_count","plat_count","subcat_count","perm_index","p_perm")]
    df_stacked <- df_stacked[, .(psuccess = pplat * (1 - prod(1-psuccess))), by=c("Platform", "pplat","cand","prot_vector_count","plat_count","subcat_count","perm_index","p_perm")]
    df_stacked <- df_stacked[, .(psuccess = (1 - prod(1-psuccess))), by=c("prot_vector_count","plat_count","subcat_count","cand","perm_index","p_perm")]
    df_stacked <- df_stacked[, .(psuccess = sum(psuccess*p_perm)), by=c("cand","prot_vector_count","plat_count","subcat_count")]
    # cost-benefit analysis
    #df_stacked$benefit<-df_stacked$psuccess*benefit*par$poverall

    if (j==1){
      # First vaccine doesn't benefit from fungibility
      df_stacked[, ben := df_stacked$psuccess * benefit * par$poverall]
      df_stacked[, mgben := ben]

      df_stacked[, mgcost := basemgcost]
      df_stacked[, cost := basemgcost]
      # df_stacked[, ofung := 0.0]
      # df_stacked[, pfung := 0.0]
    } else{
      # Every other vaccine does
      df_stacked[, ben := psuccess * benefit * par$poverall]
      df_stacked[, mgben := ben - benefitVec[j-1]]

      df_stacked[, mgcost :=
                   (1-par$fracScrap) *
                   if_else(prot_vector_count==0,
                           if_else(plat_count==0,
                                   basemgcost*(1-par$overallFung),
                                   if_else(subcat_count==0,
                                           basemgcost*(1-par$platFung),
                                           basemgcost*(1-par$subcatFung)
                                   )
                           ),
                           if_else(subcat_count==0,
                                   basemgcost*(1-par$protVectorFung),
                                   basemgcost*(1-par$subcatFung)
                           )
      )]

      # TODO: Need to double check this equation
      df_stacked[, mgcost := mgcost + par$poverall * (psuccess - culProb[j-1]) * par$fracScrap * basemgcost]

      df_stacked[, cost := costVec[j-1] + mgcost]

      # df_stacked[, pfung := if_else(df_stacked$plat_count==0,
      #                             0.0,
      #                             as.numeric(par$subcatFung))]
      # df_stacked[, ofung := as.numeric(par$overallFung)]
    }
    df_stacked[, welfareRatio := mgben / mgcost]

    # get best candidate
    stacked_I <- which.max(df_stacked$welfareRatio)
    max_I <- as.numeric(df_stacked$cand[stacked_I])

    # update
    dscPhase[, opt := (max_I == 1:.N)]
    dscPhase[opt==T, selected := selected + 1]
    selectIndices[j] <- max_I
    culProb[j] <- max(df_stacked$psuccess)*par$poverall
    costVec[j] <- df_stacked$cost[stacked_I]
    welfareVec[j] <- df_stacked$welfare[stacked_I]
    benefitVec[j] <- df_stacked$ben[stacked_I]
    # overallFungVec[j] <- df_stacked$ofung[stacked_I]
    # platFungVec[j] <- df_stacked$pfung[stacked_I]
  }

  # Data with candidates in the optimal order
  dordered <- dscPhase[selectIndices]
  dordered[, index := .I]
  dordered$CumulativeProbability<-culProb
  dordered[, MarginalProbability := CumulativeProbability - shift(CumulativeProbability, type="lag", fill=0)]

  dordered[, socialCost := costVec]
  dordered[, socialBenefit := benefitVec]
  dordered[, welfare := welfareVec]
  dordered[, oFung := overallFungVec]
  dordered[, pFung := platFungVec]
  dordered[, margcost := socialCost - shift(socialCost, type="lag", fill=0)]
  dordered[, margbenefit := socialBenefit - shift(socialBenefit, type="lag", fill=0)]

  # Random draws to compute how big competition is for every candidate ----

  # Prepare datasets
  dsubcat <- dscPhase[, c("Platform", "Subcategory")]
  dsubcat <- unique(dsubcat)
  dplat <- dscPhase[, c("Platform", "pplat")]
  dplat <- unique(dplat)
  dtarget <- dscPhase[, c("Target", "ptarget")]
  dtarget <- unique(dtarget)

  # Dataset with all replications. Cartesian product of that dataset with datesets by candidate, subcategory, and platform
  ddraws <- data.table(r=1:par$replications)
  setkey(ddraws, r)
  dplatdraws <- crossJoin(ddraws, dplat)
  setkey(dplatdraws, r, Platform)
  dsubcatdraws <- crossJoin(ddraws, dsubcat)
  setkey(dsubcatdraws, r, Platform, Subcategory)
  dtargetdraws <- crossJoin(ddraws, dtarget)
  setkey(dtargetdraws, r, Target)
  dcanddraws <- crossJoin(ddraws, dordered)
  setkey(dcanddraws, r, Platform, Subcategory)

  # Draw main random variables
  set.seed(seed)
  ddraws[, y := rbernoulli(.N, par$poverall)]
  dplatdraws[, y := rbernoulli(.N, pplat)]
  dsubcatdraws[, y := rbernoulli(.N, par$psubcat)]
  dtargetdraws[, y := rbernoulli(.N, ptarget)]
  dcanddraws[, y := rbernoulli(.N, pcand)]

  # Input random variables into main dataset with all candidates and all draws
  dcanddraws[, ysubcand := dsubcatdraws[.(dcanddraws$r, dcanddraws$Platform, dcanddraws$Subcategory), y]]
  dcanddraws[, yplat := dplatdraws[.(dcanddraws$r, dcanddraws$Platform), y]]
  dcanddraws[, ytarget := dtargetdraws[.(dcanddraws$r, dcanddraws$Target), y]]
  dcanddraws[, yoverall := ddraws[.(dcanddraws$r), y]]
  dcanddraws[, success := yoverall * yplat * ysubcand * ytarget * y]

  setkey(dcanddraws, r, index)

  # Compute expected competition for every candidate when it is the marginal candidate, conditional on success
  for (i in 1:nrow(dordered)) {
    userows <- i

    marginal <- dordered[userows]
    dmarginal <- dcanddraws[.(marginal$Platform, marginal$Subcategory, marginal$index), on=.(Platform, Subcategory, index)]
    setkey(dmarginal, r)
    dcanddraws[, marg_success := dmarginal[.(dcanddraws$r), success]]
    dcanddrawsMargSuccess <- dcanddraws[marg_success == 1 & index <= userows]
    dMargSucessOutcome <- dcanddrawsMargSuccess[, .(successes=sum(success), success=any(success)), by="r"]

    dordered[i, ExpComp := mean(1/dMargSucessOutcome$successes)]
  }

  if (computeExpComp) {
    # dexpComp <- data.table(ind=1:maxcand)

    for (n in 1:maxcand) {
      dcanddraws[, paste0("success_", n) := as.numeric(success==1 & n >= index), by="r"]
      # dcanddraws[, paste0("successes_", n) := sum(as.numeric(success==1 & n >= index)), by="r"]
      dcanddraws[, paste0("successes_", n) := sum(get(paste0("success_", n))), by="r"]

      agg <- dcanddraws[, .(exp_comp = sum(get(paste0("success_", n))/get(paste0("successes_", n)), na.rm=T) /
                              sum(get(paste0("success_", n)))), by="index"]

      dordered[, paste0("exp_comp_", n) := agg$exp_comp]
    }
  }

  dordered[, cost := cumsum(margcost)]

  return(list(dordered=dordered, dcanddraws=dcanddraws))
}

#' Combinations of capacity and candidates
#'
#' Computes social benefits, costs, and program cost for different combinations of the
#' number of candidates and capacity per candidate
#'
#' @param dordered data.table with ordered candidates
#' @param par Parameters object with model parameters
#' @param prodGrid Resolution of the grid of production capacities
#' @param maxprodfac Maximum capacity, as a factor of the capacity constraint
#'
#' @return data.table with all combinations
#' @export
capacityCandidatesCombinations <- function(dordered, par, prodGrid=0.5, maxprodfac=4) {
  margcost <- par$c * par$capacity * 12 / 1000

  combinations <- crossJoin(data.table(productionCapacity=seq(0, maxprodfac * par$capacity, prodGrid)), data.table(candidates=1:par$maxcand))
  combinations[, index := .I]
  setkey(combinations, index)

  # Computing social benefits
  combinations[, prob := dordered[.(candidates), on=.(index), CumulativeProbability]]
  combinations[, totalVaccinated := productionCapacity * par$TT / 1000]
  combinations[totalVaccinated < par$effpop, share := benefitIntegral(productionCapacity * par$TT / 1000 / par$effpop, par) *
                 par$effpop / productionCapacity * 1000 / par$TT]
  combinations[totalVaccinated >= par$effpop, share := benefitIntegral(1, par) *
                 par$effpop / totalVaccinated + (1-par$effpop / totalVaccinated)]
  combinations[, socialBenefit2 := prob * share * par$totmonthben * par$TT]
  combinations[, socialBenefit3 := prob * benefits(par$totmonthben, list(productionCapacity / 1000), par$TT, par)]

  combinations[, progBen := prob *
                 benefits(par$totmonthben, list(productionCapacity/1000, par$afterCapacity/1000),
                          c(par$TT, par$tau), par)]
  combinations[, noProgBen := prob *
                 benefits(par$totmonthben, list(1e-10, par$counterCapacity/1000), c(par$TT, par$tau), par)]
  combinations[, socialBenefit := progBen - noProgBen]

  # Computing social costs and surplus
  combinations[, baseCapacityMarginalCost := dordered[.(candidates), on=.(index), socialCost] / par$capacity]
  combinations[productionCapacity <= par$capkink, socialCost := baseCapacityMarginalCost * productionCapacity]
  combinations[productionCapacity > par$capkink, socialCost := baseCapacityMarginalCost *
                 (productionCapacity^(par$mgcostelast + 1) / par$capkink^(par$mgcostelast) +
                 par$mgcostelast * par$capkink) / (par$mgcostelast + 1)
                 ]
  combinations[, capacityCostPerCandidate := socialCost / candidates]
  combinations[, socialSurplus := socialBenefit - socialCost]

  # Computing funding needed to satisfy firms and funder's surplus
  #combinations[, push := socialCost * (1-par$fcapacity)]
  # browser()
  # combinations[, push0 := socialCost - candidates * par$fcapacity * (1-par$subcatFung) * productionCapacity / par$capacity * margcost]
  combinations[, push := socialCost - candidates * par$fcapacity * (1-par$subcatFung) * capacityCostPerCandidate]
  combinations[, psmarg := dordered[.(candidates), on=.(index), pcand * pplat * par$psubcat * par$poverall * ptarget]]
  combinations[, marketfrac := dordered[.(candidates), on=.(index), ExpComp]] # Need to fix this, since it varies by candidate
  #combinations[, capacityCost := par$c * productionCapacity * 12 / 1000]

  pscands <- dordered[, pcand * pplat * par$psubcat * par$poverall]
  minps <- function(n) {min(pscands[1:n])}
  minpscands <- sapply(1:par$maxcand, minps)

  combinations[, minpull := (productionCapacity * par$TT / 1000) * par$cprod + par$fcapacity * (1-par$subcatFung) * capacityCostPerCandidate / minpscands[candidates] / marketfrac]
  combinations[, pull := minpull + candidates * par$riskPremium * productionCapacity / 500 * par$fcapacity / 0.15]
  combinations[, progCost := push + pull]
  combinations[, funderSurplus := socialBenefit - progCost]

  return(combinations)
}

#' Optimiazation with fungible costs
#'
#' Computes the optimal policy assuming fungible costs. Also computes the optimal policy with a limited budget.
#'
#' @param dordered `data.table` with candidates ordered in the optimal selection order.
#' @param dcanddraws `data.table` with Monte Carlo draws
#' @param par Parameters object with model parameters
#'
#' @return List with several objects related to all the optimization that were performed
#' @export
optimizationsFung <- function(dordered, dcanddraws, par) {
  margcost <- par$c * par$capacity * 12 / 1000

  dcapacities <- data.table(productionCapacity=seq(0, 4 * par$capacity, 2))
  dcapacities[, totalVaccinated := productionCapacity * par$TT / 1000]
  dcapacities[totalVaccinated < par$effpop, share := benefitIntegral(productionCapacity * par$TT / 1000 / par$effpop, par) *
                par$effpop / productionCapacity * 1000 / par$TT]
  dcapacities[totalVaccinated >= par$effpop, share := benefitIntegral(1, par) *
                par$effpop / totalVaccinated + (1-par$effpop / totalVaccinated)]
  dcapacities[productionCapacity == 0, share := 0]

  dcapacities[, progBen := benefits(par$totmonthben, list(productionCapacity/1000, par$afterCapacity/1000),
                          c(par$TT, par$tau), par)]
  dcapacities[, noProgBen := benefits(par$totmonthben, list(1e-10, par$counterCapacity/1000),
                                      c(par$TT, par$tau), par)]
  dcapacities[, socialBenefitIfSuccessful := progBen - noProgBen]

  combinations <- capacityCandidatesCombinations(dordered, par)

  optindU <- combinations[which.max(socialSurplus), index]
  optindC <- combinations[productionCapacity <= par$capacity][which.max(socialSurplus), index]
  optindFU <- combinations[which.max(funderSurplus), index]
  optindFC <- combinations[productionCapacity <= par$capacity][which.max(funderSurplus), index]

  combU <- combinations[.(optindU)]
  combUl <- combinations[.(optindU-1)]
  combUf <- combinations[.(optindU+1)]
  mgcost <- (combUf$socialCost - combUl$socialCost) / (combUf$productionCapacity - combUl$productionCapacity)
  mgcostPcandPyear <- mgcost / combU$candidates / 12 * 1000
  avgcost <- combU$socialCost / combU$productionCapacity
  avgcostPcandPyear <- avgcost / combU$candidates / 12 * 1000

  # Finding optimal policy with limited budget
  dbudget <- data.table(budget=seq(0, 1.1*combinations[.(optindU), progCost], 1))
  dbudgetC <- data.table(budget=seq(0, 1.1*combinations[.(optindU), progCost], 1))
  for (i in 1:nrow(dbudget)) {
    b <- dbudget[i, budget]
    row <- combinations[progCost < b][which.max(combinations[progCost < b, socialSurplus])]
    rowC <- combinations[progCost < b & productionCapacity <= par$capacity][
      which.max(combinations[progCost < b & productionCapacity <= par$capacity, socialSurplus])]

    if (nrow(row) > 0) {
      dbudget[i, c("candidates", "productionCapacity", "prob", "push", "minpull", "pull", "progCost",
                   "socialCost", "socialBenefit", "socialSurplus") :=
                row[, .(candidates, productionCapacity, prob, push, minpull, pull, progCost,
                        socialCost, socialBenefit, socialSurplus)]]
    }
    if (nrow(rowC) > 0) {
      dbudgetC[i, c("candidates", "productionCapacity", "prob", "push", "minpull", "pull", "progCost",
                    "socialCost", "socialBenefit", "socialSurplus") :=
                 rowC[, .(candidates, productionCapacity, prob, push, minpull, pull, progCost,
                          socialCost, socialBenefit, socialSurplus)]]
    }
  }

  dbudget[, type := "Unconstrained"]
  dbudgetC[, type := "Constrained"]

  dbudget <- rbind(dbudget, dbudgetC)

  return(list(dcapacities=dcapacities, combinations=combinations, optindU=optindU, optindC=optindC,
              optindFU=optindFU, optindFC=optindFC, dbudget=dbudget,
              mgcostPcandPyear=mgcostPcandPyear, avgcostPcandPyear=avgcostPcandPyear))
}

#' Funding structure
#'
#' Computes characteristics of the funding structure, including the minimum pull amount that satisfies candidates'
#' participation constraint, every candidate's expected earnings, and the price they would have to expect to prefer to
#' choose an outside market instead of the program.
#'
#' @param data List with output from candidatesFung()
#' @param opt List with output from optimizaitonsFung()
#' @param par Parameters object with model parameters
#'
#' @return
#' @export
fundingStructure <- function(data, opt, par) {
  combinations <- opt$combinations
  dordered <- data$dordered
  dcanddraws <- data$dcanddraws
  margcost <- par$c * par$capacity * 12 / 1000
  optindC <- opt$optindC
  prodCap <- combinations[.(optindC), productionCapacity]

  # Finding details about the funding structure and firms' participation constraint
  capcost <- (1-par$subcatFung) * prodCap / par$capacity * margcost
  optcap <- combinations[.(optindC), productionCapacity]
  optCand <- combinations[.(optindC), candidates]
  push <- combinations[.(optindC), push]
  pushbase <- push - optCand * (1 - par$fcapacity) * (1-par$subcatFung) * prodCap / par$capacity * margcost

  # Computing probability that each candidate is the only succesful candidate and therefore would have an open
  # market to sell vaccines if it doesn't participate
  dcanddraws[, success_nopt := as.numeric(success==1 & optCand >= index), by="r"]
  dcanddraws[, successes_nopt := sum(success_nopt), by="r"]
  dbycand <- dcanddraws[, .(open_outside_market = sum(as.numeric(success_nopt == successes_nopt)), scenarios=.N), by="index"]
  dbycand[, prob_open_market := open_outside_market / scenarios]
  probOpenMarket <- dbycand[, prob_open_market]

  # For all candidates: success probability and expected fraction of revenue obtained
  fcands <- dcanddraws[, .(exp_comp = sum(success_nopt/successes_nopt, na.rm=T) /
                             sum(success_nopt)), by="index"]$exp_comp
  pscands <- dordered[, pcand * pplat * par$psubcat * par$poverall]
  marketfracmg <- dordered[.(optCand), on=.(index), ExpComp]
  #minfps <- function(n) {min(fcands * pscands, na.rm=T)}
  minfps <- min(fcands * pscands[1:par$maxcand], na.rm=T)

  # y intercept and slope for figure with the set of push/pull combinations that satisfy participation constraint
  pullInt <- (optcap * par$TT / 1000) * par$cprod + (capcost + pushbase/optCand) / minfps
  pullSlope <- 1/optCand / minfps
  Smax <- optCand * (minfps * (optcap * par$TT / 1000) * par$cprod + capcost) + pushbase

  # Expected earnings for each candidate
  pullMarg <- combinations[.(optindC), pull]
  pearly <- (pullMarg - par$plate * par$nlate) / par$nearly

  pull <- par$pearly * par$nearly + par$plate * par$nlate
  #candearn <- fcands * pscands * (pull - (par$nearly + par$nlate) * par$cprod) - par$fcapacity * capcost
  candearn <- fcands * pscands * (pull - (prodCap * par$TT / 1000) * par$cprod) - par$fcapacity * capcost

  # Price needed in an outside market so that it's more profitable not to participate in the program
  outsideprices <- 1/par$outquantity * (candearn / (pscands * par$outsidefactor * probOpenMarket * par$outprob) +
                                          par$outcap * par$c) + par$cprod

  return(list(candearn=candearn, outsideprices=outsideprices, pull=pull, push=push,
              pullInt=pullInt, pullSlope=pullSlope, Smax=Smax, pearly=pearly))
}

#' Benefit integral
#'
#' Computes the share of benefits obtained when a fraction `frac` of the population is vaccinated by the end of the analysis period.
#'
#' @param frac Fraction of population that has been vaccinated by the end of the period
#' @param par Parameters object with model parameters
#'
#' @return Share of benefits obtained
#' @export
benefitIntegral <- function(frac, par) {
  if (par$benefitdist == "pnorm") {
    ret <- Re(frac * (1-hypergeo(-1/par$alpha, 1/par$alpha, 1+1/par$alpha, frac^par$alpha)))
  } else if (par$benefitdist == "piecewiseLinear") {
    ret <- ifelse(frac < par$vaccshare,
                  1/2 * frac^2 * par$piecewisepar$slope1,
                  1/2 * par$vaccshare^2 * par$piecewisepar$slope1 + (frac - par$vaccshare) * par$damageshare +
                    1/2 * (frac - par$vaccshare)^2 * par$piecewisepar$slope2
    )
  } else if (par$benefitdist == "piecewiseLinearGen") {
    ret <- rep(0, length(frac))

    damageshares <- c(0, par$piecewisepar$damageshares, 1)
    vaccshares <- c(0, par$piecewisepar$vaccshares, 1)

    for (i in seq_len(length(damageshares)-1)) {
      ds <- damageshares[i]
      vs <- vaccshares[i]
      dsn <- damageshares[i+1]
      vsn <- vaccshares[i+1]

      add <- if_else(frac < vs, 0,
                     if_else(frac > vsn,
                             (vsn - vs) * ds + 1/2 * (vsn - vs) * (dsn - ds),
                             (frac - vs) * ds + 1/2 * (frac - vs)^2 * (dsn - ds) / (vsn - vs),
                             )
                     )

      #browser()

      ret <- ret + add
    }
  }

  return(ret)
}

#' Benefit integral with discounting
#'
#' More general function that allows for linear damage discounting due to introduction of treatment
#'
#' @param frac1 Lower limit of integral
#' @param frac2 Upper limit of integral
#' @param par Parameters object
#' @param share1 Share of damage at time of frac1
#' @param share2 Share of damage at time of frac2
#'
#' @return
#' @export
#'
#' @examples
benefitIntegralDisc <- function(frac1, frac2, par, share1=1, share2=1) {
  slope <- (share2 - share1) / (frac2 - frac1)
  inishare <- share1 - frac1 * slope
  endshare <- share2 + (1-frac2) * slope

  if (par$benefitdist == "pnorm") {
    int1 <- Re(1/2 * frac1 * (2 * inishare - frac1 * (inishare - endshare) -
                                2 * inishare * hypergeo(-1/par$alpha, 1/par$alpha, 1+1/par$alpha, frac1^par$alpha) +
                                (inishare - endshare) * frac1 * hypergeo(-1/par$alpha, 2/par$alpha, 1+2/par$alpha, frac1^par$alpha)))

    int2 <- Re(1/2 * frac2 * (2 * inishare - frac2 * (inishare - endshare) -
                                2 * inishare * hypergeo(-1/par$alpha, 1/par$alpha, 1+1/par$alpha, frac2^par$alpha) +
                                (inishare - endshare) * frac2 * hypergeo(-1/par$alpha, 2/par$alpha, 1+2/par$alpha, frac2^par$alpha)))

    ret <- int2 - int1

  } else if (par$benefitdist == "piecewiseLinear") {
    stop("Not implemented. Use piecewiseLinearGen")
  } else if (par$benefitdist == "piecewiseLinearGen") {
    if (frac1 > 0 | share1 < 1 | share2 < 1) {
      stop("Not implemented yet.")
    }

    ret <- rep(0, length(frac2))

    damageshares <- c(0, par$piecewisepar$damageshares, 1)
    vaccshares <- c(0, par$piecewisepar$vaccshares, 1)

    for (i in seq_len(length(damageshares)-1)) {
      ds <- damageshares[i]
      vs <- vaccshares[i]
      dsn <- damageshares[i+1]
      vsn <- vaccshares[i+1]

      add <- if_else(frac < vs, 0,
                     if_else(frac2 > vsn,
                             (vsn - vs) * ds + 1/2 * (vsn - vs) * (dsn - ds),
                             (frac - vs) * ds + 1/2 * (frac - vs)^2 * (dsn - ds) / (vsn - vs),
                     )
      )

      #browser()

      ret <- ret + add
    }
  }

  return(ret)
}


#' Benefits from vaccination
#'
#' Computes the benefits obtained from vaccination by integrating the share of harm reduction over time.
#'
#' @param monthben Benefits, i.e., total harm (economic + health) due to covid
#' @param capacity Vaccination capacity (billions per month)
#' @param endtime Time at which vaccination ends
#' @param par `Parameters` object with model parameters
#' @param capacity2 Second capacity after endtime (billions per month)
#' @param endtime2 Time at which vaccination with capacity2 ends
#'
#' @return Total benefits from vaccination
#' @export
benefitsOld <- function(monthben, capacity, endtime, par, capacity2=NULL, endtime2=NULL) {
  if ((is.null(capacity2) & !is.null(endtime2)) | (is.null(capacity2) & !is.null(endtime2))) {
    stop("The arguments capacity2 and endtime2 must either both of them or none of them be null.")
  }

  # Compute benefits from first stage of vaccination
  fracImmunized <- capacity * endtime / par$effpop
  timeVaccAll <- par$effpop / capacity

  benefits <- if_else(fracImmunized < 1, # Different computation when capacity is enough to vaccinate everyone before endtime
                      benefitIntegral(fracImmunized, par) * timeVaccAll * monthben,
                      (benefitIntegral(1, par) * timeVaccAll + (endtime - timeVaccAll)) * monthben
  )

  # Compute benefits from second stage of vaccination (if there is a second stage)
  if (!is.null(capacity2) & !is.null(endtime2)) {
    fracImmunizedNew <- fracImmunized + capacity2 * (endtime2 - endtime) / par$effpop
    timeVaccAllNew <- par$effpop / capacity2
    timeFinish <- endtime + (1-fracImmunized) * par$effpop / capacity2

    # Different computation when everyone was already vaccinated in the fist stage
    benefitsNew <- if_else(fracImmunized >= 1,
                           (endtime2 - endtime) * monthben,
                           if_else(fracImmunizedNew < 1, # Different computation when capacity is enough to vaccinate everyone before endtime2
                                   (benefitIntegral(fracImmunizedNew, par) - benefitIntegral(fracImmunized, par)) * timeVaccAllNew * monthben,
                                   ((benefitIntegral(1, par) - benefitIntegral(fracImmunized, par)) * timeVaccAllNew + (endtime2 - timeFinish)) * monthben,
                           )
    )

    benefits <- benefits + benefitsNew
  }

  return(benefits)
}

#' Benefits from vaccination
#'
#' Computes the benefits obtained from vaccination by integrating the share of harm reduction over time.
#'
#' @param monthben Benefits, i.e., total harm (economic + health) due to covid
#' @param capacities List of vectors with capacities. Each vector represents a period. Elements
#' in the vector represent the capacities for that period in each scenario.
#' @param endtimes Vector of endtimes of all the periods
#' @param par `Parameters` object with model parameters
#'
#' @return Total benefits from vaccination
#' @export
benefits <- function(monthben, capacities, endtimes, par) {

  if (length(capacities) != length(endtimes)) {
    stop("The capacities and endtimes vectors should be of the same length")
  }

  benefits <- rep(0, length(capacities[[1]]))
  fracImmunized <- rep(0, length(capacities[[1]]))
  begintime <- rep(0, length(capacities[[1]]))

  for (i in seq_along(capacities)) {
    endtime <- endtimes[i]
    capacity <- capacities[[i]]

    fracImmunizedNew <- fracImmunized + capacity * (endtime - begintime) / (par$effpop)
    timeVaccAll <- par$effpop / capacity
    timeFinish <- begintime + (1-fracImmunized) * (par$effpop) / capacity

    # Different computation when everyone was already vaccinated in the fist stage
    benefitsNew <- if_else(fracImmunized >= 1,
                           (endtime - begintime) * monthben,
                           if_else(fracImmunizedNew < 1, # Different computation when capacity is enough to vaccinate everyone before endtime2
                                   (benefitIntegral(fracImmunizedNew, par) - benefitIntegral(fracImmunized, par)) * timeVaccAll * monthben,
                                   if_else(fracImmunized < 1,
                                           ((benefitIntegral(1, par) - benefitIntegral(fracImmunized, par)) * timeVaccAll + (endtime - timeFinish)) * monthben,
                                           (endtime - timeFinish) * monthben
                                   )

                           )
                    )

    benefits <- benefits + benefitsNew
    fracImmunized <- fracImmunizedNew
    begintime <- endtime
  }

  return(benefits)
}

#' Cross join
#'
#' Cartesian product of two data.tables
#'
#' @param db1 First dataset
#' @param db2 Second dataset
#'
#' @return Cartesian product of both datasets
#' @export
crossJoin <- function(db1, db2) {
  db1[, ..dummy.. := 1]
  db2[, ..dummy.. := 1]

  db <- db1[db2, on="..dummy..", allow.cartesian=T]

  db[, ..dummy.. := NULL]
  db1[, ..dummy.. := NULL]
  db2[, ..dummy.. := NULL]

  return(db)
}

# Old functions -----

#' Optimization with fungible costs
#'
#' Computes the optimal policy assuming fungible costs. Also computes the optimal policy with a limited budget.
#'
#' @param dordered `data.table` with candidates ordered in the optimal selection order.
#' @param dcanddraws `data.table` with Monte Carlo draws
#' @param par Parameters object with model parameters
#'
#' @return List with several objects related to all the optimization that were performed
#' @export
optimizationsFungSudden <- function(dordered, dcanddraws, par) {
  margcost <- par$c * par$capacity * 12 / 1000

  dcapacities <- data.table(productionCapacity=seq(0, 4 * par$capacity, 2))
  dcapacities[, totalVaccinated := productionCapacity * par$TT / 1000]
  dcapacities[totalVaccinated < par$effpop, share := benefitIntegral(productionCapacity * par$TT / 1000 / par$effpop, par) *
                par$effpop / productionCapacity * 1000 / par$TT]
  dcapacities[totalVaccinated >= par$effpop, share := benefitIntegral(1, par) *
                par$effpop / totalVaccinated + (1-par$effpop / totalVaccinated)]
  dcapacities[productionCapacity == 0, share := 0]
  dcapacities[, socialBenefitIfSuccessful := share * par$totmonthben * par$TT]

  combinations <- crossJoin(data.table(productionCapacity=seq(0, 4 * par$capacity, 0.5)), data.table(candidates=1:par$maxcand))
  combinations[, index := .I]
  setkey(combinations, index)

  combinations[, prob := dordered[.(candidates), on=.(index), CumulativeProbability]]
  combinations[, totalVaccinated := productionCapacity * par$TT / 1000]
  combinations[totalVaccinated < par$effpop, share := benefitIntegral(productionCapacity * par$TT / 1000 / par$effpop, par) *
                 par$effpop / productionCapacity * 1000 / par$TT]
  combinations[totalVaccinated >= par$effpop, share := benefitIntegral(1, par) *
                 par$effpop / totalVaccinated + (1-par$effpop / totalVaccinated)]
  combinations[, socialBenefit := prob * share * par$totmonthben * par$TT]

  combinations[, socialCost := productionCapacity / par$capacity * dordered[.(candidates), on=.(index), socialCost]]
  combinations[, socialSurplus := socialBenefit - socialCost]

  #combinations[, push := socialCost * (1-par$fcapacity)]
  combinations[, push := socialCost - candidates * par$fcapacity * (1-par$subcatFung) * productionCapacity / par$capacity * margcost]
  combinations[, psmarg := dordered[.(candidates), on=.(index), pcand * pplat * par$psubcat * par$poverall * ptarget]]
  combinations[, marketfrac := dordered[.(candidates), on=.(index), ExpComp]] # Need to fix this, since it varies by candidate
  combinations[, capacityCost := par$c * productionCapacity * 12 / 1000]

  pscands <- dordered[, pcand * pplat * par$psubcat * par$poverall]
  minps <- function(n) {min(pscands[1:n])}
  minpscands <- sapply(1:par$maxcand, minps)

  combinations[, minpull := (productionCapacity * par$TT / 1000) * par$cprod + par$fcapacity * (1-par$subcatFung) * capacityCost / minpscands[candidates] / marketfrac]
  combinations[, pull := minpull + candidates * par$riskPremium]
  combinations[, progCost := push + pull]
  combinations[, funderSurplus := socialBenefit - progCost]

  optindU <- combinations[which.max(socialSurplus), index]
  optindC <- combinations[productionCapacity <= par$capacity][which.max(socialSurplus), index]
  optindFU <- combinations[which.max(funderSurplus), index]
  optindFC <- combinations[productionCapacity <= par$capacity][which.max(funderSurplus), index]

  # Finding optimal policy with limited budget
  dbudget <- data.table(budget=seq(0, 1.1*combinations[.(optindU), progCost], 1))
  dbudgetC <- data.table(budget=seq(0, 1.1*combinations[.(optindU), progCost], 1))
  for (i in 1:nrow(dbudget)) {
    b <- dbudget[i, budget]
    row <- combinations[progCost < b][which.max(combinations[progCost < b, socialSurplus])]
    rowC <- combinations[progCost < b & productionCapacity <= par$capacity][
      which.max(combinations[progCost < b & productionCapacity <= par$capacity, socialSurplus])]

    if (nrow(row) > 0) {
      dbudget[i, c("candidates", "productionCapacity", "prob", "push", "minpull", "pull", "progCost",
                   "socialCost", "socialBenefit", "socialSurplus") :=
                row[, .(candidates, productionCapacity, prob, push, minpull, pull, progCost,
                        socialCost, socialBenefit, socialSurplus)]]
    }
    if (nrow(rowC) > 0) {
      dbudgetC[i, c("candidates", "productionCapacity", "prob", "push", "minpull", "pull", "progCost",
                    "socialCost", "socialBenefit", "socialSurplus") :=
                 rowC[, .(candidates, productionCapacity, prob, push, minpull, pull, progCost,
                          socialCost, socialBenefit, socialSurplus)]]
    }
  }

  dbudget[, type := "Unconstrained"]
  dbudgetC[, type := "Constrained"]

  dbudget <- rbind(dbudget, dbudgetC)

  return(list(dcapacities=dcapacities, combinations=combinations, optindU=optindU, optindC=optindC,
              optindFU=optindFU, optindFC=optindFC, dbudget=dbudget))
}

#' Compute optimal candidate table
#'
#' Computes the optimal ordering of candidates to maximize the probability of getting at least one vaccine.
#' Also computes a few additional statistics that are necessary for rest of the model.
#'
#' @param d data.table with summarized data
#' @param par Parameters object with all parameters for the model
#' @param computeExpComp
#'
#' @return data.table with the optimal candidates ordered by the order in which they should be chosen
#' @export
computeCandidateTable <- function(d, par, computeExpComp=F) {
  if (par$pspike==1 & par$precombinant==1 & par$potherprotein==1 ){
    return(computeCandidateTable_nested(d, par, computeExpComp=F))
  }

  # Adding platform feasibility probabilities to table
  d[Platform == "DNA", pplat := as.numeric(par$pdna)]
  d[Platform == "RNA", pplat := par$prna]
  d[Platform == "Live attenuated virus", pplat := par$pattenuated]
  d[Platform == "Viral vector", pplat := par$pvector]
  d[Platform == "Protein subunit", pplat := par$psubunit]
  d[Platform == "Inactivated", pplat := par$pinactivated]
  d[Platform == "VLP", pplat := par$pvlp]
  d[Platform == "Dendritic cells", pplat := par$pdendritic]
  d[Platform == "Self-assembling vaccine", pplat := par$psav]
  d[Platform == "Unknown", pplat := par$punknown]
  d[Platform == "Artificial antigen presenting cells", pplat := par$paapc]
  d[Platform == "Live-attenuated bacteria", pplat := par$plivebac]
  d$ptarget<-1.0
  d[Target == "Other", ptarget := par$potherprotein]
  d[Target == "Spike", ptarget := par$pspike]
  d[Target == "Recombinant", ptarget := par$precombinant]

  # Create summary table with candidates by subcategory and phase of trials
  dscPhase <- crossJoin(d, data.table(phase=c("Pre-clinical", "Phase 1", "Phase 2", "Repurposed")))
  dscPhase[phase=="Pre-clinical", Candidates := PreClinicalCandidates]
  dscPhase[phase=="Phase 1", Candidates := Phase1candidates]
  dscPhase[phase=="Phase 2", Candidates := Phase2candidates]
  dscPhase[phase=="Repurposed", Candidates := RepurposedCandidates]

  dscPhase[, PreClinicalCandidates := NULL]
  dscPhase[, Phase1candidates := NULL]
  dscPhase[, Phase2candidates := NULL]
  dscPhase[, RepurposedCandidates := NULL]
  dscPhase$pcand<-1.0
  dscPhase[phase == "Pre-clinical", pcand := par$ppreclinical]
  dscPhase[phase == "Phase 1", pcand := par$pphase1]
  dscPhase[phase == "Phase 2", pcand := par$pphase2]
  dscPhase[phase == "Repurposed", pcand := par$prepurposed]

  # Finding optimal order of vaccines

  # Start with no candidates
  dscPhase[, selected := 0]
  dscPhase[, psuccess := 1-(1-pcand)^selected]

  # Vector to store indices of selected vaccines
  selectIndices <- rep(NA, par$maxcand)

  # Vector to store cumulative probabilities
  culProb <- rep(NA, par$maxcand)

  # Only vaccines that have at least one candidate
  dscPhase <- dscPhase[Candidates > 0]

  #Change NAs to 0s
  dscPhase[is.na(dscPhase$pplat), pplat := 0]

  # Sorting data and creating index variable
  setkey(dscPhase, Platform, Subcategory, phase)
  setindex(dscPhase, Platform, Subcategory)
  dscPhase[, index := .I]

  # Compute possible target protein success combinations
  targets<- c("Spike","Recombinant","Other")
  probs<- c(as.numeric(par$pspike), as.numeric(par$precombinant), as.numeric(par$potherprotein))
  perm<-permutations(2,length(targets),v=c(0,1),repeats.allowed=TRUE)
  p_perm<- vector(mode="numeric", length=nrow(perm))
  temp<-matrix(0,nrow=nrow(perm),ncol=ncol(perm))
  for (i in 1:length(probs)){
    temp[,i]<- if_else(perm[,i]==1,probs[i], 1-probs[i])
  }
  for (i in 1:length(p_perm)){
    p_perm[i]<- prod(temp[i,])
  }
  perm<-as.data.frame(perm)
  perm<-perm[p_perm!=0,]
  names(perm)<-targets
  perm$p_perm<-p_perm[p_perm!=0]
  perm$perm_index <- seq.int(nrow(perm))
  perm<-data.table(perm)


  # Main loop. Sequentially find the next candidate

  for (j in 1:par$maxcand) {
    # get candidates
    dcand <- data.table(cand=dscPhase$index[dscPhase$selected<dscPhase$Candidates])

    # stack data frames of increments
    df_stacked <- crossJoin(dcand, dscPhase)
    df_stacked[index==cand, selected:=selected+1]

    # stack data frames of possible permutations
    df_stacked <- crossJoin(perm, df_stacked)

    # drop according to permutation
    df_stacked<- df_stacked[!(Spike==0 & Target=="Spike")]
    df_stacked<- df_stacked[!(Recombinant==0 & Target=="Recombinant")]
    df_stacked<- df_stacked[!(Other==0 & Target=="Other")]

    # get probabilities
    df_stacked[, psuccess := 1-(1-pcand)^selected]
    df_stacked <- df_stacked[, .(psuccess = par$psubcat * (1 - prod(1-psuccess))), by=c("Platform", "Subcategory", "pplat","cand","perm_index","p_perm")]
    df_stacked <- df_stacked[, .(psuccess = pplat * (1 - prod(1-psuccess))), by=c("Platform", "pplat","cand","perm_index","p_perm")]
    df_stacked<-df_stacked[, .(psuccess = (1 - prod(1-psuccess))), by=c("cand","perm_index","p_perm")]
    df_stacked<-df_stacked[, .(psuccess = sum(psuccess*p_perm)), by=c("cand")]

    # get best candidate
    max_I<-as.numeric(df_stacked$cand[which.max(df_stacked$psuccess)])

    # update
    dscPhase[, opt := (max_I == 1:.N)]
    dscPhase[opt==T,"selected"] <- dscPhase[opt==T,"selected"] + 1
    selectIndices[j] <- max_I
    culProb[j]<-max(df_stacked$psuccess)*par$poverall

    dscPhase[,plat_count:= sum(selected),by=c("Platform")]
    dscPhase[plat_count==1 & opt==T,]

  }

  # Data with candidates in the optimal order
  dordered <- dscPhase[selectIndices]
  dordered[, index := .I]
  dordered[, CumulativeProbability := culProb]
  dordered[, MarginalProbability := CumulativeProbability - shift(CumulativeProbability, type="lag", fill=0)]

  # Compute variables about progress given optimal order

  dordered[, ind_in_platform := 1:.N, by="Platform"]
  dordered[, cum_candidates := cumsum(rep(1,.N)), by=c("Platform", "Subcategory")]

  # Random draws to compute how big competition is for every candidate

  # Prepare datasets
  dsubcat <- dscPhase[, c("Platform", "Subcategory")]
  dsubcat <- unique(dsubcat)
  dplat <- dscPhase[, c("Platform", "pplat")]
  dplat <- unique(dplat)
  dtarget <- dscPhase[, c("Target", "ptarget")]
  dtarget <- unique(dtarget)

  # Dataset with all replications. Cartesian product of that dataset with datesets by candidate, subcategory, and platform
  ddraws <- data.table(r=1:par$replications)
  setkey(ddraws, r)
  dplatdraws <- crossJoin(ddraws, dplat)
  setkey(dplatdraws, r, Platform)
  dsubcatdraws <- crossJoin(ddraws, dsubcat)
  setkey(dsubcatdraws, r, Platform, Subcategory)
  dtargetdraws <- crossJoin(ddraws, dtarget)
  setkey(dtargetdraws, r, Target)
  dcanddraws <- crossJoin(ddraws, dordered)
  setkey(dcanddraws, r, Platform, Subcategory)

  # Draw main random variables
  ddraws[, y := rbernoulli(.N, par$poverall)]
  dplatdraws[, y := rbernoulli(.N, pplat)]
  dsubcatdraws[, y := rbernoulli(.N, par$psubcat)]
  dtargetdraws[, y := rbernoulli(.N, ptarget)]
  dcanddraws[, y := rbernoulli(.N, pcand)]

  # Input random variables into main dataset with all candidates and all draws
  dcanddraws[, ysubcand := dsubcatdraws[.(dcanddraws$r, dcanddraws$Platform, dcanddraws$Subcategory), y]]
  dcanddraws[, yplat := dplatdraws[.(dcanddraws$r, dcanddraws$Platform), y]]
  dcanddraws[, ytarget := dtargetdraws[.(dcanddraws$r, dcanddraws$Target), y]]
  dcanddraws[, yoverall := ddraws[.(dcanddraws$r), y]]
  dcanddraws[, success := yoverall * yplat * ysubcand * ytarget * y]

  setkey(dcanddraws, r, index)

  # Compute expected competition for every candidate when it is the marginal candidate, conditional on success
  for (i in 1:nrow(dordered)) {
    userows <- i

    marginal <- dordered[userows]
    dmarginal <- dcanddraws[.(marginal$Platform, marginal$Subcategory, marginal$index), on=.(Platform, Subcategory, index)]
    setkey(dmarginal, r)
    dcanddraws[, marg_success := dmarginal[.(dcanddraws$r), success]]
    dcanddrawsMargSuccess <- dcanddraws[marg_success == 1 & index <= userows]
    dMargSucessOutcome <- dcanddrawsMargSuccess[, .(successes=sum(success), success=any(success)), by="r"]

    dordered[i, ExpComp := mean(1/dMargSucessOutcome$successes)]
  }

  if (computeExpComp) {
    # dexpComp <- data.table(ind=1:maxcand)

    for (n in 1:maxcand) {
      dcanddraws[, paste0("success_", n) := as.numeric(success==1 & n >= index), by="r"]
      # dcanddraws[, paste0("successes_", n) := sum(as.numeric(success==1 & n >= index)), by="r"]
      dcanddraws[, paste0("successes_", n) := sum(get(paste0("success_", n))), by="r"]

      agg <- dcanddraws[, .(exp_comp = sum(get(paste0("success_", n))/get(paste0("successes_", n)), na.rm=T) /
                              sum(get(paste0("success_", n)))), by="index"]

      dordered[, paste0("exp_comp_", n) := agg$exp_comp]
    }
  }

  # Total benefit
  benefit <- par$damageint * par$totmonthben * par$TT

  # Marginal cost of developing a vaccine
  mgcost <- par$c * par$capacity * 12 / 1000

  dordered[, margcost := mgcost * ifelse(ind_in_platform==1, 1, 1-par$fcostred)]
  dordered[, margbenefit := MarginalProbability * benefit]
  dordered[, cost := cumsum(margcost)]

  return(list(dordered=dordered, dcanddraws=dcanddraws))
}


computeCandidateTable_nested <- function(d, par, computeExpComp=F, seed=1) {

  # Adding platform feasibility probabilities to table
  d[Platform == "DNA", pplat := as.numeric(par$pdna)]
  d[Platform == "RNA", pplat := par$prna]
  d[Platform == "Live attenuated virus", pplat := par$pattenuated]
  d[Platform == "Viral vector", pplat := par$pvector]
  d[Platform == "Protein subunit", pplat := par$psubunit]
  d[Platform == "Inactivated", pplat := par$pinactivated]
  d[Platform == "VLP", pplat := par$pvlp]
  d[Platform == "Dendritic cells", pplat := par$pdendritic]
  d[Platform == "Self-assembling vaccine", pplat := par$psav]
  d[Platform == "Unknown", pplat := par$punknown]
  d[Platform == "Artificial antigen presenting cells", pplat := par$paapc]
  d[Platform == "Live-attenuated bacteria", pplat := par$plivebac]

  d$ptarget<-1.0
  d[Target == "Other", ptarget := par$potherprotein]
  d[Target == "Spike", ptarget := par$pspike]
  d[Target == "Recombinant", ptarget := par$precombinant]

  # Create summary table with candidates by subcategory and phase of trials
  dscPhase <- crossJoin(d, data.table(phase=c("Pre-clinical", "Phase 1", "Phase 2", "Repurposed")))
  dscPhase[phase=="Pre-clinical", Candidates := PreClinicalCandidates]
  dscPhase[phase=="Phase 1", Candidates := Phase1candidates]
  dscPhase[phase=="Phase 2", Candidates := Phase2candidates]
  dscPhase[phase=="Repurposed", Candidates := RepurposedCandidates]

  dscPhase[, PreClinicalCandidates := NULL]
  dscPhase[, Phase1candidates := NULL]
  dscPhase[, Phase2candidates := NULL]
  dscPhase[, RepurposedCandidates := NULL]

  dscPhase[phase == "Pre-clinical", pcand := par$ppreclinical]
  dscPhase[phase == "Phase 1", pcand := par$pphase1]
  dscPhase[phase == "Phase 2", pcand := par$pphase2]
  dscPhase[phase == "Repurposed", pcand := par$prepurposed]

  # Finding optimal order of vaccines ----

  # Start with no candidates
  dscPhase[, selected := 0]
  dscPhase[, psuccess := 1-(1-pcand)^selected]

  # Vector to store indices of selected vaccines
  selectIndices <- rep(NA, par$maxcand)

  # Only vaccines that have at least one candidate
  dscPhase <- dscPhase[Candidates > 0]

  # Sorting data and creating index variable
  setkey(dscPhase, Platform, Subcategory, phase)
  setindex(dscPhase, Platform, Subcategory)
  dscPhase[, index := .I]

  # Main loop. Sequentially find the next candidate
  for (j in 1:par$maxcand) {
    # Find the optimal candidate within each subcategory
    dscPhase[selected < Candidates, log1mp := -log(1 - pcand)]
    dscPhase[selected == Candidates, log1mp := 0]
    dscPhase[, opt := (which.max(log1mp) == 1:.N), by=c("Platform", "Subcategory")]
    dscPhase[selected == Candidates, opt := F]
    dscPhase[, Nselected := selected]
    dscPhase[opt==T, Nselected := selected + 1]
    dscPhase[, Npsuccess := 1-(1-pcand)^Nselected]

    # Find the optimal candidate within each platform
    dsubcat <- dscPhase[, .(psuccess = par$psubcat * (1 - prod(1-psuccess)), Npsuccess = par$psubcat * (1 - prod(1-Npsuccess))),
                        by=c("Platform", "Subcategory", "pplat")]
    dsubcat[, opt := (which.max((1-psuccess)/(1-Npsuccess)) == 1:.N), by=c("Platform")]
    dsubcat[opt==F, Npsuccess := psuccess]
    setkey(dsubcat, Platform, Subcategory)

    # Find the optimal platform
    dplat <- dsubcat[, .(psuccess = pplat * (1 - prod(1-psuccess)), Npsuccess = pplat * (1 - prod(1-Npsuccess))),
                     by=c("Platform", "pplat")]
    dplat[, opt := (which.max((1-psuccess)/(1-Npsuccess)) == 1:.N)]
    setkey(dplat, Platform)

    # Recover optimal platform, subcategory, and phase
    opt_plat <- dplat[opt == T, Platform]
    opt_subcat <- dsubcat[.(opt_plat)][opt == T, Subcategory]
    opt_phase <- dscPhase[.(opt_plat, opt_subcat), on=.(Platform, Subcategory)][opt == T, phase]

    # Increase by one the number of selected candidates for the optimal candidate
    dscPhase[.(opt_plat, opt_subcat, opt_phase), selected := selected + 1]
    dscPhase[.(opt_plat, opt_subcat, opt_phase), psuccess := 1-(1-pcand)^selected]

    selectIndices[j] <- dscPhase[.(opt_plat, opt_subcat, opt_phase), index]

    # print(paste0(dscPhase[.(opt_plat, opt_subcat, opt_phase), Subcategory]))
  }

  # Data with candidates in the optimal order
  dordered <- dscPhase[selectIndices]
  dordered[, index := .I]

  # Compute variables about progress given optimal order ----

  dordered[, ind_in_platform := 1:.N, by="Platform"]
  dordered[, cum_candidates := cumsum(rep(1,.N)), by=c("Platform", "Subcategory")]

  # Cumulative probabilities within subcategory
  dordered[, cum_subcat_prob := par$psubcat * (1 - exp(cumsum(log(1-pcand)))), by=c("Platform", "Subcategory")]
  dordered[, prev_cum_subcat_prob := shift(cum_subcat_prob, type="lag", fill=0), by=c("Platform", "Subcategory")]

  # Main loop through the rows of the table
  for (i in 1:nrow(dordered)) {
    # Cumulative probabilities within platform
    if (dordered[i, ind_in_platform] == 1) {
      dordered[i, prev_plat_cumprob := 0]
      dordered[i, plat_cumprob := pplat * cum_subcat_prob ]
    } else {
      dordered[i, prev_plat_cumprob :=
                 {
                   pl <- Platform;
                   iip <- ind_in_platform - 1;
                   ret <- dordered[.(pl, iip), on=.(Platform, ind_in_platform), plat_cumprob];
                   ret
                 }
               ]
      dordered[i, plat_comp_cumprob := (1-prev_plat_cumprob/pplat) * (1-cum_subcat_prob) / (1-prev_cum_subcat_prob) ]
      dordered[i, plat_cumprob := pplat * (1 - plat_comp_cumprob) ]
    }
  }

  # Overall cumulative probability
  dordered[1, cumprob := par$poverall * plat_cumprob]
  for (i in 2:nrow(dordered)) {
    dordered[i, prev_cumprob := dordered[i-1, cumprob]]
    dordered[i, comp_cumprob := (1-prev_cumprob/par$poverall) * (1-plat_cumprob) / (1-prev_plat_cumprob)]
    dordered[i, cumprob := par$poverall * (1-comp_cumprob)]
  }

  dordered[, c("cum_candidates", "cum_subcat_prob", "prev_cum_subcat_prob", "prev_plat_cumprob",
               "plat_cumprob", "plat_comp_cumprob", "prev_cumprob", "comp_cumprob") := NULL]

  setnames(dordered, "cumprob", "CumulativeProbability")

  # Computing marginal probabilities
  dordered[, MarginalProbability := CumulativeProbability - shift(CumulativeProbability, type="lag", fill=0)]

  # Random draws to compute how big competition is for every candidate ----

  # Dataset with all replications. Cartesian product of that dataset with datesets by candidate, subcategory, and platform
  ddraws <- data.table(r=1:par$replications)
  setkey(ddraws, r)
  dplatdraws <- crossJoin(ddraws, dplat)
  setkey(dplatdraws, r, Platform)
  dsubcatdraws <- crossJoin(ddraws, dsubcat)
  setkey(dsubcatdraws, r, Platform, Subcategory)
  dcanddraws <- crossJoin(ddraws, dordered)
  setkey(dcanddraws, r, Platform, Subcategory)

  # Draw main random variables
  set.seed(seed)
  ddraws[, y := rbernoulli(.N, par$poverall)]
  dplatdraws[, y := rbernoulli(.N, pplat)]
  dsubcatdraws[, y := rbernoulli(.N, par$psubcat)]
  dcanddraws[, y := rbernoulli(.N, pcand)]

  # Input random variables into main dataset with all candidates and all draws
  dcanddraws[, ysubcand := dsubcatdraws[.(dcanddraws$r, dcanddraws$Platform, dcanddraws$Subcategory), y]]
  dcanddraws[, yplat := dplatdraws[.(dcanddraws$r, dcanddraws$Platform), y]]
  dcanddraws[, yoverall := ddraws[.(dcanddraws$r), y]]
  dcanddraws[, success := yoverall * yplat * ysubcand * y]

  setkey(dcanddraws, r, index)

  # Compute expected competition for every candidate when it is the marginal candidate, conditional on success
  for (i in 1:nrow(dordered)) {
    userows <- i

    marginal <- dordered[userows]
    dmarginal <- dcanddraws[.(marginal$Platform, marginal$Subcategory, marginal$index), on=.(Platform, Subcategory, index)]
    setkey(dmarginal, r)
    dcanddraws[, marg_success := dmarginal[.(dcanddraws$r), success]]
    dcanddrawsMargSuccess <- dcanddraws[marg_success == 1 & index <= userows]
    dMargSucessOutcome <- dcanddrawsMargSuccess[, .(successes=sum(success), success=any(success)), by="r"]

    dordered[i, ExpComp := mean(1/dMargSucessOutcome$successes)]
  }

  if (computeExpComp) {
    # dexpComp <- data.table(ind=1:maxcand)

    for (n in 1:maxcand) {
      dcanddraws[, paste0("success_", n) := as.numeric(success==1 & n >= index), by="r"]
      # dcanddraws[, paste0("successes_", n) := sum(as.numeric(success==1 & n >= index)), by="r"]
      dcanddraws[, paste0("successes_", n) := sum(get(paste0("success_", n))), by="r"]

      agg <- dcanddraws[, .(exp_comp = sum(get(paste0("success_", n))/get(paste0("successes_", n)), na.rm=T) /
                              sum(get(paste0("success_", n)))), by="index"]

      dordered[, paste0("exp_comp_", n) := agg$exp_comp]
    }
  }

  # Total benefit
  benefit <- par$damageint * par$totmonthben * par$TT

  # Marginal cost of developing a vaccine
  mgcost <- par$c * par$capacity * 12 / 1000

  dordered[, margcost := mgcost * ifelse(ind_in_platform==1, 1, 1-par$fcostred)]
  dordered[, margbenefit := MarginalProbability * benefit]
  dordered[, cost := cumsum(margcost)]

  return(list(dordered=dordered, dcanddraws=dcanddraws))
}


#' Optimal number of candidates
#'
#' Based on a dataset of candidates ordered optimally, compute optimal number of candidates as well as a few other numbers.
#'
#' @param data data.table with ordered candidates
#' @param par Parameters object with all parameters for the model
#'
#' @return List with the optimal number of candidates, the marginal probability at the optimum, the marginal
#' cost at the optimum, and the benefit of finding a vaccine
#' @export
optCandidates <- function(dordered, par) {
  # Total benefit
  #benefit <- par$damageint * par$totmonthben * par$TT

  # Marginal cost of developing a vaccine
  mgcost <- par$c * par$capacity * 12 / 1000

  # Value of the optmial marginal probability
  #optMargProb <- mgcost / benefit

  # Find the row with the optimal number of candidates
  # dordered[, NextMarginalProbability := shift(MarginalProbability, type="lead", fill=0)]
  # relrow <- dordered[MarginalProbability >= optMargProb & optMargProb > NextMarginalProbability]
  dordered[, nextMargBenefit := shift(margbenefit, type="lead", fill=0)]
  relrow <- dordered[margbenefit >= margcost & margcost > nextMargBenefit]

  optCand <- relrow[, index]
  optProb <- relrow[, CumulativeProbability]
  socialBenefit <- relrow[, socialBenefit]
  socialCost <- relrow[, cost]

  return(list(optCand=optCand, margcost=dordered$margcost,
              margbenefit=dordered$margbenefit, optProb=optProb, mgcost=mgcost,
              socialBenefit=socialBenefit, socialCost=socialCost))
}

optCandidatesB <- function(dordered, par) {
  # Total benefit
  benefit <- par$damageint * par$totmonthben * par$TT

  # Marginal cost of developing a vaccine
  mgcost <- par$c * par$capacity * 12 / 1000

  # Value of the optmial marginal probability
  #optMargProb <- mgcost / benefit

  # Find the row with the optimal number of candidates
  # dordered[, NextMarginalProbability := shift(MarginalProbability, type="lead", fill=0)]
  # relrow <- dordered[MarginalProbability >= optMargProb & optMargProb > NextMarginalProbability]
  dordered[, nextMargBenefit := shift(margbenefit, type="lead", fill=0)]
  relrow <- dordered[margbenefit >= margcost & margcost > nextMargBenefit]

  optCand <- relrow[, index]
  optProb <- relrow[, CumulativeProbability]
  socialBenefit <- optProb * benefit
  socialCost <- relrow[, cost]

  return(list(optCand=optCand, margcost=dordered$margcost,
              margbenefit=dordered$margbenefit, optProb=optProb, mgcost=mgcost,
              socialBenefit=socialBenefit, oFung=dordered$oFung, pFung=dordered$pFung))
}


#' Program structure
#'
#' Computes some characteristics of the optimal program
#'
#' @param dordered data.table with candidates with the optimal order
#' @param dcanddraws data.table with draws from Monte Carlo simulation
#' @param optCand Optimal number of candidates
#' @param par Parameters object with model parameters
#'
#' @return List with pull and push amounts, cost of capacity, intercept of the pull quantity, slope of the pull quantity,
#' optimal number of candidates, and maximum push amount
#' @export
progStructure <- function(dordered, dcanddraws, optCand, par) {

  avgCapcost <- par$capcost * (optCand - par$fcostred * (optCand - length(unique(dordered[1:optCand, Platform])))) / optCand

  # Total push quantity the program must pay
  push <- (optCand - par$fcostred * (optCand - length(unique(dordered[1:optCand, Platform])))) * (1-par$fcapacity) * par$capcost

  # For the marginal candidate: success probability and expected fraction of revenue obtained
  psmarg <- dordered[optCand, pcand * pplat * par$psubcat * par$poverall ]
  fmarg <- dordered[optCand, ExpComp]

  # Amount of pull funding needed for the marginal candidate to break even, and to justify the risky investment
  minpull <- (par$capacity * par$TT / 1000) * par$cprod + par$fcapacity * avgCapcost / psmarg / fmarg
  #pull <- minpull + optCand * par$riskPremium
  #pearly <- (pull - par$plate * par$nlate) / par$nearly
  pull <- par$pearly * par$nearly + par$plate * par$nlate
  pearly <- par$pearly

  # y intercept and slope for figure with the set of push/pull combinations that satisfy participation constraint
  pullInt <- (par$capacity * par$TT / 1000) * par$cprod + avgCapcost / psmarg / fmarg
  pullSlope <- 1/optCand / psmarg / fmarg
  Smax <- optCand * (psmarg * fmarg * (par$capacity * par$TT / 1000) * par$cprod + avgCapcost)

  # Computing probability that each candidate is the only succesful candidate and therefore would have an open
  # market to sell vaccines if it doesn't participate
  dcanddraws[, success_nopt := as.numeric(success==1 & optCand >= index), by="r"]
  dcanddraws[, successes_nopt := sum(success_nopt), by="r"]
  dbycand <- dcanddraws[, .(open_outside_market = sum(as.numeric(success_nopt == successes_nopt)), scenarios=.N), by="index"]
  dbycand[, prob_open_market := open_outside_market / scenarios]
  probOpenMarket <- dbycand[, prob_open_market]

  # For all candidates: success probability and expected fraction of revenue obtained
  fcands <- dcanddraws[, .(exp_comp = sum(success_nopt/successes_nopt, na.rm=T) /
                             sum(success_nopt)), by="index"]$exp_comp
  pscands <- dordered[, pcand * pplat * par$psubcat * par$poverall]

  # Expected earnings for each candidate
  candearn <- fcands * pscands * (pull - (par$nearly + par$nlate) * par$cprod) - par$fcapacity * par$capcost

  # Price needed in an outside market so that it's more profitable not to participate in the program
  outcap <- 0.1 # Billions
  outprob <- 0.2
  outquantity <- 0.5 # Billions
  outsidefactor <- 2
  outsideprices <- 1/outquantity * (candearn / (pscands * outsidefactor * probOpenMarket * outprob) + outcap * par$c) + par$cprod

  return(list(pull=pull, push=push, pearly=pearly, minpull=minpull, pullInt=pullInt, pullSlope=pullSlope,
              optCand=optCand, Smax=Smax, candearn=candearn, outsideprices=outsideprices))
}

#' Optimal number of candidates
#'
#' Based on a dataset of candidates ordered optimally, compute optimal number of candidates as well as a few other numbers.
#'
#' @param data data.table with ordered candidates
#' @param par Parameters object with model parameters
#'
#' @return List with the optimal number of candidates, the marginal probability at the optimum, the marginal
#' cost at the optimum, and the benefit of finding a vaccine
#' @export
optCandidatesCapacity <- function(data, par) {

  # damageint=0.5, TT=6, monthben=375, mortality=200000, statLife=1.18, lifeExp=71, yearsLost=5, pdrug=0.5, sharm=0.5,
  #                         c=1, capacity=500, cprod=1) {
  # Health benefit of a vaccine per month
  hben <- par$mortality * par$statLife * par$yearsLost / par$lifeExp / 1000

  # Marginal cost of developing a vaccine
  margcost <- par$c * par$capacity * 12 / 1000

  maxT <- 6

  combinations <- crossJoin(data.table(TT=seq(0, maxT, 0.01)), data.table(candidates=1:30))
  combinations[, prob := data[.(combinations$candidates), on=.(index), CumulativeProbability]]
  combinations[, socialBenefit := prob * maxT * par$totmonthben - (1-par$damageint) * prob * TT * par$totmonthben]
  combinations[, socialCost := candidates * margcost / TT * 6]
  combinations[, socialSurplus :=  socialBenefit - socialCost]

  optrow <- combinations[which.max(combinations$socialSurplus)]

  optCand <- optrow$candidates
  optTT <- optrow$TT
  optCapacity <- par$capacity * maxT / optTT
  optCost <- optrow$socialCost
  optBenefit <- optrow$socialBenefit
  optSurplus <- optrow$socialSurplus

  candOptimization <- combinations[TT==optTT]
  setkey(candOptimization, candidates)
  candOptimization[, prevBenefit := shift(socialBenefit, type="lag", fill=0)]
  candOptimization[, margBenefit := socialBenefit - prevBenefit]
  candOptimization[, margCost := margcost / TT * 6]

  TToptimization <- combinations[candidates==optCand]
  setkey(TToptimization, TT)
  TToptimization[, prevBenefit := shift(socialBenefit, type="lag", fill=0)]
  TToptimization[, margBenefit := socialBenefit - prevBenefit]
  TToptimization[, prevCost := shift(socialCost, type="lag", fill=NA)]
  TToptimization[, margCost := socialCost - prevCost]
  TToptimization[, productionCapacity := par$capacity * maxT / TT]
  TToptimization[, prevProductionCapacity := shift(productionCapacity, type="lag", fill=0)]
  TToptimization[, margBenefitCap := margBenefit / (productionCapacity - prevProductionCapacity) ]
  TToptimization[, margCostCap := margCost / (productionCapacity - prevProductionCapacity) ]
  TToptimization[, margCostCapB := par$c * 12 / 1000 * candidates]

  return(list(optCand=optCand, optTT=optTT, optCost=optCost, optBenefit=optBenefit, optSurplus=optSurplus,
              candOptimization=candOptimization,
              TToptimization=TToptimization, optCapacity=optCapacity, combinations=combinations))
}

optBudgetB <- function(dordered, dcanddraws, par) {

  # Marginal cost of developing a vaccine
  #margcost <- par$c * par$capacity * 12 / 1000

  maxT <- 6

  combinations <- crossJoin(data.table(productionCapacity=seq(par$capacity, 5 * par$capacity, 2)), data.table(candidates=1:30))
  combinations[, TT := par$capacity * maxT / productionCapacity]
  combinations[, prob := dordered[.(candidates), on=.(index), CumulativeProbability]]
  combinations[, socialBenefit := prob * maxT * par$totmonthben - (1-par$damageint) * prob * TT * par$totmonthben]
  combinations[, socialCost := candidates * margcost / TT * 6]
  combinations[, socialSurplus :=  socialBenefit - socialCost]

  combinations[, push := socialCost * (1-par$fcapacity)]
  combinations[, psmarg := dordered[.(candidates), on=.(index), pcand * pplat * par$psubcat * par$poverall *ptarget]]
  combinations[, fmarg := dordered[.(candidates), on=.(index), ExpComp]]
  combinations[, minpull := (par$capacity * par$TT / 1000) * par$cprod + par$fcapacity * par$capcost *oFung*pFung / psmarg / fmarg]
  combinations[, pull := minpull + candidates * par$riskPremium]
  combinations[, progCost := push + pull]

  dbudget <- data.table(budget=seq(0, 200, 1))
  dbudgetC <- data.table(budget=seq(0, 200, 1))
  for (i in 1:nrow(dbudget)) {
    b <- dbudget[i, budget]
    row <- combinations[progCost < b][which.max(combinations[progCost < b, socialSurplus])]
    rowC <- combinations[progCost < b & productionCapacity <= par$capacity][
      which.max(combinations[progCost < b & productionCapacity <= par$capacity, socialSurplus])]

    if (nrow(row) > 0) {
      dbudget[i, c("TT", "candidates", "productionCapacity", "prob", "push", "minpull", "pull", "progCost",
                   "socialCost", "socialBenefit", "socialSurplus") :=
                row[, .(TT, candidates, productionCapacity, prob, push, minpull, pull, progCost,
                        socialCost, socialBenefit, socialSurplus)]]
    }
    if (nrow(rowC) > 0) {
      dbudgetC[i, c("TT", "candidates", "productionCapacity", "prob", "push", "minpull", "pull", "progCost",
                    "socialCost", "socialBenefit", "socialSurplus") :=
                 rowC[, .(TT, candidates, productionCapacity, prob, push, minpull, pull, progCost,
                          socialCost, socialBenefit, socialSurplus)]]
    }
  }

  dbudget[, type := "Unconstrained"]
  dbudgetC[, type := "Constrained"]

  dbudget <- rbind(dbudget, dbudgetC)

  return(list(combinations=combinations, dbudget=dbudget))
}


#' Optimal policy under a budget constraint
#'
#' Compute the optimal policy when there is a budget constraint. It computes the optima both with and without
#' a capacity constraint.
#'
#' @param dordered Ordered dataset with candidate vaccines
#' @param dcanddraws Monte Carlo draws of vaccines
#' @param par Parameter object with model parametesr
#'
#' @return A list with all the potential combinations of candidates and capacity, and a list with the optimal policiy for
#' different levels of the budget
#' @export
optBudget <- function(dordered, dcanddraws, par) {
  # Marginal cost of developing a vaccine
  margcost <- par$c * par$capacity * 12 / 1000

  maxT <- 6

  combinations <- crossJoin(data.table(productionCapacity=seq(par$capacity, 5 * par$capacity, 2)), data.table(candidates=1:30))
  combinations[, TT := par$capacity * maxT / productionCapacity]
  combinations[, prob := dordered[.(candidates), on=.(index), CumulativeProbability]]
  combinations[, socialBenefit := prob * maxT * par$totmonthben - (1-par$damageint) * prob * TT * par$totmonthben]
  combinations[, socialCost := candidates * margcost / TT * 6]
  combinations[, socialSurplus :=  socialBenefit - socialCost]

  combinations[, push := socialCost * (1-par$fcapacity)]
  combinations[, psmarg := dordered[.(candidates), on=.(index), pcand * pplat * par$psubcat * par$poverall *ptarget]]
  combinations[, fmarg := dordered[.(candidates), on=.(index), ExpComp]]
  combinations[, minpull := (par$capacity * par$TT / 1000) * par$cprod + par$fcapacity * par$capcost / psmarg / fmarg]
  combinations[, pull := minpull + candidates * par$riskPremium]
  combinations[, progCost := push + pull]

  dbudget <- data.table(budget=seq(0, 200, 1))
  dbudgetC <- data.table(budget=seq(0, 200, 1))
  for (i in 1:nrow(dbudget)) {
    b <- dbudget[i, budget]
    row <- combinations[progCost < b][which.max(combinations[progCost < b, socialSurplus])]
    rowC <- combinations[progCost < b & productionCapacity <= par$capacity][
      which.max(combinations[progCost < b & productionCapacity <= par$capacity, socialSurplus])]

    if (nrow(row) > 0) {
      dbudget[i, c("TT", "candidates", "productionCapacity", "prob", "push", "minpull", "pull", "progCost",
                   "socialCost", "socialBenefit", "socialSurplus") :=
                row[, .(TT, candidates, productionCapacity, prob, push, minpull, pull, progCost,
                        socialCost, socialBenefit, socialSurplus)]]
    }
    if (nrow(rowC) > 0) {
      dbudgetC[i, c("TT", "candidates", "productionCapacity", "prob", "push", "minpull", "pull", "progCost",
                    "socialCost", "socialBenefit", "socialSurplus") :=
                 rowC[, .(TT, candidates, productionCapacity, prob, push, minpull, pull, progCost,
                          socialCost, socialBenefit, socialSurplus)]]
    }
  }

  dbudget[, type := "Unconstrained"]
  dbudgetC[, type := "Constrained"]

  dbudget <- rbind(dbudget, dbudgetC)

  return(list(combinations=combinations, dbudget=dbudget))
}
