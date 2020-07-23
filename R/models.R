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
loadData <- function(par, candidateFile=NULL, includeVaccines=c()) {
  # We set encoding to BOM UTF to avoid cross-platform issues

  if (is.null(candidateFile)) {
    if (par$inputfile=="Default") {
      d <- data.table(read.csv("Data/vaccinesSummary.csv", fileEncoding = "UTF-8-BOM"))
    } else if (par$inputfile=="US") {
      d <- data.table(read.csv("Data/vaccinesSummaryUS.csv", fileEncoding = "UTF-8-BOM"))
    }
  } else {
    d <- data.table(read.csv(candidateFile, fileEncoding = "UTF-8-BOM"))
  }

  d[is.na(PreClinicalCandidates), PreClinicalCandidates := 0]
  d[is.na(Phase1candidates), Phase1candidates := 0]
  d[is.na(Phase2candidates), Phase2candidates := 0]
  d[is.na(RepurposedCandidates), RepurposedCandidates := 0]

  if ("X" %in% names(d)) {
    d[, X := NULL]
  }

  # dchoose <- data.table(read.csv("Data/vaccinesChoose.csv", fileEncoding = "UTF-8-BOM"))
  # d <- d[.(dchoose$Platform, dchoose$Subcategory), on=.(Platform, Subcategory)]

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
#' @param seed Seed for random generator
#'
#' @return List with two data.tables. `dordered` lists all the candidates in the optimal order and with some statistics.
#' `dcanddraws` includes all the random draws and their outcomes.
#' @export
#' @import gtools
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
  perm <- permutations(2,length(targets),v=c(0,1),repeats.allowed=TRUE)
  p_perm <- vector(mode="numeric", length=nrow(perm))
  temp <- matrix(0,nrow=nrow(perm),ncol=ncol(perm))
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

#' Benefit integral
#'
#' Computes the share of benefits obtained when a fraction `frac` of the population is vaccinated by the end of the analysis period.
#'
#' @param frac Fraction of population that has been vaccinated by the end of the period
#' @param par Parameters object with model parameters
#'
#' @return Share of benefits obtained
#' @import hypergeo
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
#' @return Values of benefits integral
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
