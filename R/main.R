#' Optimal portfolio for a price taker
#'
#' Computes the optimal portfolio for a price-taking country
#'
#' @param population Country population (in millions)
#' @param gdp_pc Country GDP per capita (in thousand $)
#' @param frac_high_risk Fraction of population that is high risk
#' @param loss2yr Cumulative percent of GDP lost because of pandemic over two years
#' @param price Market price for capacity ($ per course / year)
#' @param steps Steps to optimize over
#' @param candidateFile File with candidate data
#'
#' @return Vector with optimal capacities
#' @export
portfolioPriceTaker <- function(population, gdp_pc, frac_high_risk, loss2yr, price, steps=c(10,1,0.1), candidateFile=NULL) {
  par <- Parameters$new(global=F, inputfile="Default", maxcand=30, monthben=500, popshare=population/7800,
                         gdpshare=population*gdp_pc/1e3/87.3, fracHighRisk=frac_high_risk, afterCapacity=population/7800*500,
                         counterCapacity=population/7800*500, econlossratio=loss2yr/0.138)

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

  # Initialize at zero
  capacities <- rep(0, nrow(dcandidate))

  # Optimize over grids with decreasing step sizes
  for (s in steps) {
    capacities <- optimizeGrid(capacities, objectiveFun, step=s, price=price, par=par)
  }

  # Computing other values
  totCapacity <- sum(capacities)
  cost <- priceTakerCost(capacities, price)
  expBenefits <- countryExpectedBenefits(capacities, dcandidate, targetPermutations, dplatforms, par, grid=s)

  return(list(capacities=capacities, totCapacity=totCapacity, cost=cost, expBenefits=expBenefits))
}
