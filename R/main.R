#' Optimal portfolio for a price taker
#'
#' Computes the optimal portfolio for a price-taking country.
#' This is done by calling optimisation function [optimizeGrid()]
#' by optimising using decreasing step sizes (see `step`)
#' on benefits function [countryNetBenefits()]
#'
#' @param parameters `Parameters` object with model parameters
#' @param population Country population (in millions)
#' @param gdp_pc Country GDP per capita (in thousand $)
#' @param frac_high_risk Fraction of population that is high risk
#' @param loss2yr Cumulative percent of GDP lost because of pandemic over two years
#' @param price Market price for manufacturing capacity ($ per course / year)
#' @param steps Steps to optimize over. By default three step-sizes are used in succession:
#'              10, 1 and 0.1.
#' @param candidateFile File with candidate data
#' @param return_benefit_args if `TRUE`, then an extra object is returned, that can then be used
#'                            as input to [countryNetBenefits()]
#' @param lambda Lagrange multiplier. Should be one, unless trying to solve
#'               problem with constrained budget.
#'
#' @return List with information on optimal portfolio. Includes investment
#'  in every candidate, benefits, cost, and total capacity.
#' @export
#' @examples
#' \dontrun{
#' population <- 31.99
#' gdp_pc <- 6.71
#' frac_high_risk <- 0.131
#' loss2yr <- 0.269
#' price <- 4
#' portfolio <- portfolioPriceTaker(population=population, gdp_pc=gdp_pc,
#'                    frac_high_risk=frac_high_risk, loss2yr=loss2yr, price=price)
#' }
portfolioPriceTaker <- function(parameters=NULL, population, gdp_pc, frac_high_risk, loss2yr,
                                price, steps=c(10,1,0.1), candidateFile=NULL, lambda=1,
                                return_benefit_args = FALSE) {
  . <- Platform <- pplat <- prob <- capacity <- NULL

  if (is.null(parameters)) {
    # Create parameters from function arguments
    par <- Parameters$new(population=population, gdp_pc=gdp_pc,
                          frac_high_risk=frac_high_risk, loss2yr=loss2yr)
  } else {
    # Use parameters object passed as an argument
    par <- parameters
  }

  d <- loadData(par, candidateFile)

  d$Target <- "Other"
  d$Target[1:5]<-"Spike"
  d$Target[10:15]<-"Recombinant"

  dordered <- candidatesFung(d, par)$dordered

  dordered <- dordered[,1:11]
  dplatforms <- unique(dordered[, .(Platform, pplat)])
  setkey(dplatforms, Platform, pplat)
  dcandidate <- copy(dordered)

  # Getting information about permutations
  targets <- c("Spike","Recombinant","Other")
  probs <- c(as.numeric(par$pspike),
             as.numeric(par$precombinant),
             as.numeric(par$potherprotein))
  targetPermutations <- getTargetPermutations(targets, probs)

  # Creating objective function
  objectiveFun <- function(capacities, grid, price, par) {
    countryNetBenefits(capacities, dcandidate, targetPermutations,
                       dplatforms, grid, price, par, lambda=lambda)
  }

  # Initialize at zero
  capacities <- rep(0, nrow(dcandidate))

  # Optimize over grids with decreasing step sizes
  for (s in steps) {
    capacities <- optimizeGrid(capacities, objectiveFun, step=s, price=price, par=par,
                               name = "country net benefits")
  }

  if(return_benefit_args)
    benefit_args <- list(capacities=capacities, dcandidate=dcandidate, targetPermutations=targetPermutations,
                         dplatforms=dplatforms, grid=min(steps), price=price, par=par, lambda=lambda)

  # Computing other values
  totCapacity <- sum(capacities)
  cost <- priceTakerCost(capacities, price)
  expBenefits <-
    countryExpectedBenefits(capacities, dcandidate,
                            targetPermutations, dplatforms, par, grid=s)
  distribution <-
    countryDistribution(capacities, dcandidate,
                        targetPermutations, dplatforms, par, grid=s)
  expCapacity <- sum(distribution[, prob * capacity])

  out <- list(capacities=capacities, totCapacity=totCapacity, cost=cost,
              expBenefits=expBenefits, distribution=distribution, expCapacity=expCapacity,
              dordered=dordered)

  if(return_benefit_args)
    out$benefit_args <-  benefit_args

  return(out)
}


#' Demand curve for a price-taking country
#'
#' Builds the demand curve for a price-taking country
#'
#' @param parameters `Parameters` object with model parameters
#' @param population Country population (in millions)
#' @param gdp_pc Country GDP per capita (in thousand $)
#' @param frac_high_risk Fraction of population that is high risk
#' @param loss2yr Cumulative percent of GDP lost because of pandemic over two years
#' @param prices Vector of market price for manufacturing capacity ($ per course / year)
#' @param candidateFile File with candidate data
#' @param inisteps Step sizes to optimize over for the initial point
#' @param mainstep Step size for main optimization
#' @param candidateFile File with candidate data
#' @param verbose How much output to produce. 0 means no output, 1 means limited output, 2 means full output.
#'
#' @return List with information on demand curve. Includes a `data.table` with total demand, benefits, and cost,
#' and a matrix with demand for individual candidates at every price
#' @export
#' @examples
#' \dontrun{
#' population <- 31.99
#' gdp_pc <- 6.71
#' frac_high_risk <- 0.131
#' loss2yr <- 0.269
#' par <- Parameters$new(population=population, gdp_pc=gdp_pc,
#'                    frac_high_risk=frac_high_risk,
#'                    loss2yr=loss2yr)
#'                    demand <- demandPriceTaker(parameters=par)
#' }
demandPriceTaker <- function(parameters=NULL, population, gdp_pc, frac_high_risk, loss2yr,
                             prices=seq(100,1,-1), inisteps=c(10,1), mainstep=0.1, candidateFile=NULL, verbose=0) {

  . <- Platform <- pplat <- ind <- expBenefit <- cost <- expCapacity <- prob <- capacity <-
    successProb <- totalCapacity <- NULL

  if (is.null(parameters)) {
    # Create parameters from function arguments
    par <- Parameters$new(global=F, inputfile="Default", maxcand=30, monthben=500, popshare=population/7800,
                          gdpshare=population*gdp_pc/1e3/87.3, fracHighRisk=frac_high_risk, afterCapacity=population/7800*500,
                          counterCapacity=population/7800*500, loss2yr=loss2yr)
  } else {
    # Copy parameters object passed as an argument
    par <- parameters
  }

  d <- loadData(par, candidateFile)

  d$Target <- "Other"
  d$Target[1:5]<-"Spike"
  d$Target[10:15]<-"Recombinant"

  dordered <- candidatesFung(d, par)$dordered

  dordered <- dordered[,1:11]
  dplatforms <- unique(dordered[, .(Platform, pplat)])
  setkey(dplatforms, Platform, pplat)
  dcandidate <- copy(dordered)

  # Getting information about permutations
  targets <- c("Spike","Recombinant","Other")
  probs <- c(as.numeric(par$pspike), as.numeric(par$precombinant), as.numeric(par$potherprotein))
  targetPermutations <- getTargetPermutations(targets, probs)

  # Creating objective function
  objectiveFun <- function(capacities, grid, price, par) {
    countryNetBenefits(capacities, dcandidate, targetPermutations, dplatforms, grid, price, par, lambda=1)
  }

  # Initialize capacities at zero, initialize price
  capacities <- rep(0, nrow(dcandidate))
  price <- prices[1]

  # Optimize over grids with decreasing step sizes to get initial point
  for (s in inisteps) {
    capacities <- optimizeGrid(capacities, objectiveFun, verbose=verbose, step=s, price=price, par=par)
  }

  # Setting up objects to save infomation
  optimizations <- data.table(price=prices)
  optimizations[, ind := .I]
  setkey(optimizations, ind)
  allCapacities <- matrix(0, length(prices), length(capacities))

  # Loop over prices to get demand levels
  for (p in prices[1:length(prices)]) {
    capacities <- optimizeGrid(capacities, objectiveFun, verbose=verbose, step=mainstep, price=p, par=par)

    i <- optimizations[price==p, ind]
    allCapacities[i, ] <- capacities
    optimizations[.(i), expBenefit :=
                    countryExpectedBenefits(capacities, dcandidate, targetPermutations, dplatforms, par, grid=mainstep)]
    optimizations[.(i), cost := priceTakerCost(capacities, p)]

    distribution <- countryDistribution(capacities, dcandidate, targetPermutations, dplatforms, par, grid=mainstep)

    optimizations[.(i), expCapacity := sum(distribution[, prob * capacity])]
    optimizations[.(i), successProb := sum(distribution[capacity > 0.011, prob])]
    optimizations[.(i), totalCapacity := sum(capacities)]
  }

  return(list(optimizations=optimizations, allCapacities=allCapacities, dordered=dordered))
}
